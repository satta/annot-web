class HardWorker
  include Sidekiq::Worker
  sidekiq_options :retry => 1
  include Sidekiq::Status::Worker

  def expiration
    @expiration ||= 60*60*24*30*12*15 # keep jobs around for a long time
  end

  def with_environment(variables={})
    if block_given?
      old_values = variables.map{ |k,v| [k,ENV[k]] }
      begin
         variables.each{ |k,v| ENV[k] = v }
         result = yield
      ensure
        old_values.each{ |k,v| ENV[k] = v }
      end
      result
    else
      variables.each{ |k,v| ENV[k] = v }
    end
  end

  def add_result_file(job, file)
    file_path = "#{job.job_directory}/#{file}"
    if File.exist?(file_path) then
      rf = ResultFile.new
      rf.file = File.new(file_path)
      rf.md5 = Digest::MD5.file(file_path).hexdigest
      rf.save!
      job.result_files << rf
    else
      raise "required file not present: #{job.job_directory}/#{file}"
    end
  end

  def perform(id)
    # wait a bit to minimize timing issues
    Kernel.sleep(5)
    job = Job.find(id)
    store name: job[:name]
    store job_id: @jid

    read, write = IO.pipe
    # we need to run this in a separate process to be able to kill it
    fork_pid = fork do
      read.close

      my_stdout = nil
      my_stderr = nil
      status = nil

      # start job
      job[:started_at] = DateTime.now
      job[:job_id] = @jid
      job.save!

      # send job start notification email
      if job[:email] and job[:email].length > 0 then
        JobMailer.start_job_email(job).deliver_later
      end

      begin
        # make necessary run directories
        job.make_directories!

        # prepare starting Sidekiq job
        cf = JobsHelper::ConfigFactory.new
        uf = job.sequence_file
        cf.use_target_seq(uf)
        r = Reference.find(job[:reference_id])
        cf.select_reference(r)
        cf.use_prefix(job[:prefix])
        cf.do_contiguation(job[:do_contiguate] && r.has_chromosomes?)
        if job[:use_transcriptome_data] then
          tf = job.transcript_file
          if tf then
            cf.use_transcript_file(tf)
          end
        end
        cf.make_embl()
        cf.use_reference()
        cf.do_circos()
        cf.do_exonerate(job[:do_exonerate])
        cf.run_ratt(job[:do_ratt])
        cf.do_pseudo(job[:do_pseudo])
        job[:config_file] = cf.get_file(job).path
        job.save!

        # Nextflow run
        run = "#{CONFIG['nextflowpath']}/nextflow -c " + \
              "#{CONFIG['locationconfig']} -c " + \
              "#{job[:config_file]} run " + \
              "#{CONFIG['nextflowscript']} #{CONFIG['dockerconf']} " + \
              "#{'-resume' unless job[:no_resume]} " + \
              "--dist_dir #{job.job_directory}"
        Rails.logger.info run
        nf_pid, stdin, stdout, stderr = nil, nil, nil, nil
        with_environment("ROOTDIR" => "#{CONFIG['rootdir']}",
                         "NXF_WORK" => job.work_directory,
                         "NXF_TEMP" => job.temp_directory) do
          nf_pid, stdin, stdout, stderr = Open4::popen4(run)
        end

        # once we're here, we can set up the code to interrupt a running job
        Signal.trap("INT") do
          if nf_pid then
            # if Nextflow is still running, stop it
            Process.kill("INT", nf_pid)
            # collect output from pipes
            stdout.readlines
            stderr.readlines
            stdout.close
            stdin.close
            # collect exit status
            Process::waitpid(nf_pid)
          end
          exit
        end

        loop do
          Kernel.sleep(5)
          if !Process::waitpid(nf_pid, Process::WNOHANG).nil? then
            puts "Nextflow finished with exit code #{$?} and output #{read.read}!"
            my_stdout = stdout.readlines.join
            my_stderr = stderr.readlines.join
            break
          end
        end

        Rails.logger.info "STDOUT:"
        Rails.logger.info my_stdout
        Rails.logger.info "STDERR:"
        Rails.logger.info my_stderr

        job[:stderr] = my_stderr
        job[:stdout] = my_stdout
        job.save!

        # job finished successfully but no annotations were created
        if status == 0 and !File.exist?("#{job.job_directory}/pseudo.out.gff3") then
          gstat = GenomeStat.new
          gstat[:nof_genes] = 0
          gstat.save!
          job.genome_stat = gstat
          job.save!
        else
          # zip up EMBL files
          Kernel.system("cd #{job.job_directory}; tar -czf #{job.job_directory}/embl.tar.gz *.embl")
          if File.exist?("#{job.job_directory}/embl.tar.gz") then
            Kernel.system("rm -f #{job.job_directory}/*.embl")
          end

          add_result_file(job, "pseudochr.fasta.gz")
          add_result_file(job, "pseudo.out.gff3")
          add_result_file(job, "pseudo.pseudochr.agp")
          add_result_file(job, "scafs.fasta.gz")
          add_result_file(job, "scaffold.out.gff3")
          add_result_file(job, "pseudo.scafs.agp")
          add_result_file(job, "embl.tar.gz") if File.exist?("#{job.job_directory}/embl.tar.gz")
          add_result_file(job, "out.gaf")
          add_result_file(job, "proteins.fasta")
          job.save!

          if (not status) or status.exitstatus != 0 then
            raise "run failed at nextflow stage"
          end

          # stats gathering
          if File.exist?("#{job.job_directory}/stats.txt") then
            gstat = GenomeStat.new
            File.open("#{job.job_directory}/stats.txt").read.each_line do |l|
              l.chomp!
              m = /^([^:]+):\s+(.+)$/.match(l)
              if m and gstat.has_attribute?(m[1].to_sym)
                gstat[m[1].to_sym] = m[2]
              end
            end
            gstat.save!
            job.genome_stat = gstat
            job.save!
          end

          # store circos images
          Dir.glob("#{job.job_directory}/chr*.png") do |f|
            m = f.match(/chr(.+)\.png$/)
            if not m.nil? then
              chrname = m[1]
            else
              chrname = "unnamed chromosome"
            end
            img = CircosImage.new
            img.file = File.new(f)
            img.chromosome = chrname
            img.save!
            job.circos_images << img
            File.unlink(f)
          end

          # import genes
          if File.exist?("#{job.job_directory}/genelist.csv") then
            genes = []
            File.open("#{job.job_directory}/genelist.csv").read.each_line do |l|
              l.chomp!
              id, type, product, seqid, start, stop, strand = l.split("\t")
              g = Gene.new(:gene_id => id, :product => product, :loc_start => start,
                           :loc_end => stop, :strand => strand, :job => job,
                           :seqid => seqid, :gtype => type, :species => job[:prefix],
                           :section => r[:section])
              genes << g
            end
            Gene.import(genes)
          end

          # import clusters
          if File.exist?("#{job.job_directory}/orthomcl_out") then
            File.open("#{job.job_directory}/orthomcl_out").each_line do |l|
              m = l.match(/^(ORTHOMCL[0-9]+)[^:]+:\s+(.+)/)
              next unless m
              c = Cluster.new(:cluster_id => m[1], :job => job)
              r = m[2].scan(/([^ ()]+)\([^)]+\)/)
              r.each do |memb|
                # HACK! needs to be done correctly for all possible transcript namings!
                memb_id = memb[0].gsub(/(:.+$|\.\d+$|\.mRNA$)/,"")
                g = Gene.where(["gene_id LIKE ? AND (job_id = #{job[:id]} OR job_id IS NULL)", "#{memb_id}%"]).take
                if g then
                  c.genes << g
                else
                  Rails.logger.info("#{memb[0]}: #{memb_id} (with job ID #{job[:id]}) not found!")
                end
              end
              c.save!
            end
          end

          # store tree files
          if File.exist?("#{job.job_directory}/tree_selection.genes") and
            File.exist?("#{job.job_directory}/tree.aln") then
            genes = []
            t = Tree.new(job: job)
            t.seq = File.new("#{job.job_directory}/tree.aln")
            t.save!
            File.open("#{job.job_directory}/tree_selection.genes").each_line do |l|
              l.chomp!
              l.split(/\s+/).each do |memb|
                # HACK! needs to be done correctly for all possible transcript namings!
                memb_id = memb.gsub(/(:[^:]+$|\.\d+$)/,"")
                g = Gene.where(["gene_id LIKE ? AND (job_id = #{job[:id]} OR job_id IS NULL)", "#{memb_id}%"]).take
                if g then
                  t.genes << g
                else
                  Rails.logger.info("#{memb}: #{memb_id} (with job ID #{job[:id]}) not found!")
                end
              end
              t.save!
            end
          end
        end

        job[:finished_at] = DateTime.now
        job.save!

        # remove temporary directories
        job.cleanup_work_directories!

        # send finish notification email
        if job[:email] and job[:email].length > 0 then
          JobMailer.finish_success_job_email(job).deliver_later
        end

        write.write(job[:stdout])
      rescue => e
        job[:finished_at] = DateTime.now
        errstr = "#{e.backtrace.first}: #{e.message} (#{e.class})\n"
        errstr += e.backtrace.drop(1).map{|s| "\t#{s}"}.join("\n")
        if job[:stderr] then
          job[:stderr] = job[:stderr] + "\n" + errstr
        else
          job[:stderr] = errstr
        end

        write.write(job[:stderr])
        job.save!
        # send error notification email
        if job[:email] and job[:email].length > 0 then
          JobMailer.finish_failure_job_email(job).deliver_later
        end
        JobMailer.finish_failure_job_email_notify_dev(job).deliver_later
        raise e
      end
    end

    write.close
    puts "actual worker PID: #{fork_pid}"
    loop do
      # wait for another round
      Kernel.sleep(5)

      # if the job has finished, just return
      if !Process::waitpid(fork_pid, Process::WNOHANG).nil? then
        Rails.logger.info("finished, exit code #{$?}:\n#{read.read}!")
        break
      end

      # kill Nextflow process and return if the job was cancelled
      if self.cancelled? then
        puts "killing forked worker with PID #{fork_pid}"
        Process.kill("INT", fork_pid)
        Process::waitpid(fork_pid)
        job[:finished_at] = DateTime.now
        job.save!
        # remove temporary directories
        job.cleanup_work_directories!
        puts "finished with exit code #{$?} and output #{read.read}!"
        break
      end
    end
  end

  def cancelled?
    Sidekiq.redis {|c| c.exists("cancelled-#{@jid}") }
  end

  def self.cancel!(jid)
    Sidekiq.redis {|c| c.setex("cancelled-#{jid}", 86400, 1) }
  end
end