class JobsController < ApplicationController
  def new
    if @closed then
      flash[:info] = "New job creation is temporarily closed for technical " + \
                     "reasons."
      redirect_to :welcome
    else
      # set defaults
      @job = Job.new
      @job[:do_contiguate] = true
      @job[:do_exonerate] = false
      @job[:do_ratt] = true
      @job[:do_pseudo] = true
      @job[:use_transcriptome_data] = false
      @job[:no_resume] = false
      @job[:max_gene_length] = 50000
      @job[:abacas_match_size] = 500
      @job[:abacas_match_sim] = 85
      @job[:augustus_score_threshold] = 0.8
      @job[:taxon_id] = 5653
      @job[:db_id ] = "Companion"
      @job[:ratt_transfer_type] = 'Species'
      @job.build_sequence_file
      @job.build_transcript_file
    end
  end

  def index
    if !logged_in? then
      flash[:info] = "You do not have permission to view all jobs in the " + \
                     "queue. Please access your job using the URL you were " + \
                     "given or by searching for your job ID."
      redirect_to :welcome
    else
      jobs = Job.all
      @outjobs = []
      @running = 0
      if jobs then
        jobs.each do |job|
          jobname = Sidekiq::Status::get(job[:job_id], :name)
          if not jobname then
            jobname = job[:name]
          end
          @outjobs <<  {:job_id => job[:job_id],
                        :id => job[:id],
                        :created_at => job[:created_at],
                        :started_at => job[:started_at],
                        :finished_at => job[:finished_at],
                        :name => jobname,
                        :email => job[:email],
                        :queued => Sidekiq::Status::queued?(job[:job_id]),
                        :working => Sidekiq::Status::working?(job[:job_id]),
                        :failed => Sidekiq::Status::failed?(job[:job_id]),
                        :complete => Sidekiq::Status::complete?(job[:job_id])}
          @running = @outjobs.reduce(0) do |a,job|
            if job[:working] then
              a + 1
            else
              a
            end
          end
        end
      end
    end
  end

  def create
    @job = Job.new(jobs_params(params))
    if simple_captcha_valid? && @job.save then
      Rails.logger.info @job.inspect
      file = @job.sequence_file
      file.file.canonicalize_seq!
      # start job and record Sidekiq ID
      jobid = HardWorker.perform_async(@job[:id])
      @job[:job_id] = jobid
      @job.save!
      url = Rails.application.routes.url_helpers.job_url(id: @job[:job_id], :host => request.host_with_port)
      flash[:success] = "Your job with ID <b>#{jobid}</b> has just been successfully
           created and enqueued. You can go back to this page and check the
           status of your job using the following URL:
           <a href=\"#{url}\">#{url}</a>."
      redirect_to job_path(id: @job[:job_id])
    else
      flash[:info] = "Invalid human verification code. Please try again."
      render 'new'
    end
  end

  def cancel
    if !logged_in? then
      flash[:info] = "You do not have permission to cancel jobs in the " + \
                     "queue."
      redirect_to :welcome
    else
      thisjob = Job.find_by(:job_id => params[:id])
      if thisjob then
        HardWorker.cancel!(thisjob[:job_id])
        queue = Sidekiq::Queue.new
        queue.each do |job|
          job.delete if job.jid == params[:id]
        end
        if File.exist?("#{thisjob.job_directory}") then
          FileUtils.rm_rf("#{thisjob.job_directory}")
        end
        flash[:info] = "Job '#{thisjob[:name]}' was cancelled."
        if logged_in? then
          redirect_to :jobs
        else
          redirect_to "welcome/index"
        end
      else
        render plain: "job #{params[:id]} not found", status: 404
      end
    end
  end

  def destroy
    if !logged_in? then
      flash[:info] = "You do not have permission to delete jobs in the " + \
                     "queue."
      redirect_to :welcome
    elsif CONFIG['example_job_id'] == params[:id] then
      flash[:info] = "Can't delete the example job. Remove the job from " + \
                      "the configuration to allow deletion."
      redirect_to :jobs
    else
      thisjob = Job.find_by(:job_id => params[:id])
      if thisjob then
        if thisjob.sequence_file then
          thisjob.sequence_file.destroy
        end
        if thisjob.transcript_file then
          thisjob.transcript_file.destroy
        end
        flash[:info] = "Job '#{thisjob[:name]}' was deleted."
        queue = Sidekiq::Queue.new
        queue.each do |job|
          job.delete if job.jid == params[:id]
        end
        if File.exist?("#{thisjob.job_directory}") then
          FileUtils.rm_rf("#{thisjob.job_directory}")
        end
      end
      thisjob.destroy
      if logged_in? then
        redirect_to :jobs
      else
        redirect_to "welcome/index"
      end
    end
  end

  def destroy_by_int_id
    if !logged_in? then
      flash[:info] = "You do not have permission to delete jobs in the " + \
                     "queue."
      redirect_to :welcome
    else
      thisjob = Job.find_by(:job_id => params[:id])
      if thisjob then
        flash.now[:info] = "Deleted job #{thisjob[:name]}"
        thisjob.destroy
        queue = Sidekiq::Queue.new
        queue.each do |job|
          job.delete if job.jid == params[:id]
        end
        if File.exist?("#{thisjob.job_directory}") then
          FileUtils.rm_rf("#{thisjob.job_directory}")
        end
      end
      render "welcome/index"
    end
  end

  def show
    @job = Job.find_by(:job_id => params[:id])
    if not @job then
      flash.now[:danger] = "The job with ID '#{params[:id]}' could not be found."
      render "welcome/index" , status: 404
    else
      @job_hash = {:queued => Sidekiq::Status::queued?(@job[:job_id]),
                   :working => Sidekiq::Status::working?(@job[:job_id]),
                   :failed => Sidekiq::Status::failed?(@job[:job_id]),
                   :complete => Sidekiq::Status::complete?(@job[:job_id])}
      @ref = Reference.find(@job[:reference_id])
      if @job_hash[:failed] then
        render 'jobs/show_failed'
      elsif @job_hash[:complete] and (!@job.genome_stat or @job.genome_stat[:nof_genes] == 0) then
        render 'jobs/show_empty'
      else
        render 'jobs/show'
      end
    end
  end

  def orths
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if not job then
      render plain: "job #{params[:id]} not found", status: 404
    else
      ref = Reference.find(job[:reference_id])
      cl_a = job.clusters.joins(:genes).where(:genes => {:species => job[:prefix], :gtype => 'coding'}).distinct.collect {|c| c[:cluster_id]}
      #cl_a = Cluster.find_by_sql("select * from clusters c join clusters_genes cg on cg.cluster_id = c.id join genes g on g.id = cg.gene_id where c.job_id = '#{job[:id]}' and g.species = '#{job[:prefix]}'")
      cl_b = job.clusters.joins(:genes).where(:genes => {:species => ref[:abbr], :gtype => 'coding'}).distinct.collect {|c| c[:cluster_id]}
      #cl_b = Cluster.find_by_sql("select * from clusters c join clusters_genes cg on cg.cluster_id = c.id join genes g on g.id = cg.gene_id where c.job_id = ''#{job[:id]}'' and g.species = '#{ref[:abbr]}'")
      a = cl_a - cl_b
      b = cl_b - cl_a
      ab = cl_a & cl_b
      orths = [{:name => {:A => job[:prefix], :B => ref[:abbr]},
                :data => {:A => a, :B => b, :AB => ab}}]
      render json: orths
    end
  end

  def orths_for_cluster
    clusters = params[:cluster]
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if not job then
      render plain: "job #{params[:id]} not found", status: 404
    else
      ref = Reference.find(job[:reference_id])
      # the code below is horribly hacky IMHO and should eventually be
      # replaced by proper SQL or iterative set operations
      cl_a = job.clusters.joins(:genes).where(:genes => {:species => job[:prefix], :gtype => 'coding'}).distinct
      cl_b = job.clusters.joins(:genes).where(:genes => {:species => ref[:abbr], :gtype => 'coding'}).distinct
      v = {}
      v[job[:prefix]] = cl_a - cl_b
      v[ref[:abbr]] = cl_b - cl_a
      v["#{job[:prefix]} #{ref[:abbr]}"] = cl_a & cl_b
      results = []
      if not v[clusters] then
        render plain: "clusters #{params[:cluster]} not valid", status: 500
      else
        v[clusters].each do |cl|
          cl.genes.each do |g|
            results << {id: g[:gene_id], product: g[:product], cluster: cl[:cluster_id]}
          end
        end
        respond_to do |format|
          format.html do
            out = []
            results.each do |r|
              out << "#{r[:id]}\t#{r[:product]}\t#{r[:cluster]}"
            end
            render plain: out.join("\n")
          end
          format.json do
            render json: results
          end
        end
      end
    end
  end

  def get_clusters
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if job and File.exist?("#{job.job_directory}/orthomcl_out") then
      render file: "#{job.job_directory}/orthomcl_out", layout: false, \
        content_type: 'text/plain'
    else
      render plain: "job #{params[:id]} not found", status: 404
    end
  end

  def get_singletons
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if job then
      ref = Reference.find(job[:reference_id])
      this_s = job.genes.includes(:clusters).where(:clusters => {id: nil}, :gtype => 'coding')
      #ref_s = Gene.where(job: nil, species: ref[:abbr]).includes(:clusters).where(:clusters => { id: nil})
      ref_s = Gene.where(species: ref[:abbr], :gtype => 'coding') -
       Gene.joins('JOIN clusters_genes ON genes.id = clusters_genes.gene_id ' + \
                  'JOIN clusters ON clusters_genes.cluster_id = clusters.id')
           .where('genes.species = ? and clusters.job_id = ?', ref[:abbr], job[:id])
      respond_to do |format|
        format.html do
          out = []
          results.each do |r|
            out << "#{r[:id]}\t#{r[:product]}\t#{r[:cluster]}"
          end
          render plain: out.join("\n")
        end
        format.json do
          render json: {:ref => ref_s, :this => this_s}
        end
      end
    else
      render plain: "job #{params[:id]} not found", status: 404
    end
  end

  def get_tree
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if job and File.exist?("#{job.job_directory}/tree.out") then
      data = File.open("#{job.job_directory}/tree.out").read
      send_data data, :filename => "#{job[:job_id]}.nwk"
    else
      render plain: "job #{params[:id]} not found or completed", status: 404
    end
  end

  def get_report
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if job then
      if File.exist?("#{job.job_directory}/pseudo.report.html") then
        data = File.open("#{job.job_directory}/pseudo.report.html").read
        send_data data, :filename => "#{job[:job_id]}.pseudo.report.html"
      else
        render html: "No report found for this job."
      end
    else
      render plain: "job #{params[:id]} not found or completed", status: 404
    end
  end

  def get_tree_genes
    job = Job.find_by(:job_id => params[:id])
    expires_in 1.month, :public => true
    if job then
      render json: job.tree.genes.order(:product)
    else
      render plain: "job #{params[:id]} not found", status: 404
    end
  end

  def get_all_synteny_images
    job = Job.find_by(:job_id => params[:id])
    if job then
      if job.circos_images and job.circos_images.size > 0 then
        t = Tempfile.new(job.job_id)
        Zip::OutputStream.open(t.path) do |zos|
          job.circos_images.each do |image|
            zos.put_next_entry("#{image.chromosome}.png")
            zos.print IO.read(image.file.path)
          end
        end
        send_file t.path, :type => 'application/zip',
                          :disposition => 'attachment',
                          :filename => "#{job.name.parameterize}.zip"
        t.close
      else
        render html: "No synteny plots found for this job."
      end
    else
      render plain: "job #{params[:id]} not found or completed", status: 404
    end
  end

  private

  def jobs_params(params)
    params.require(:job).permit(:name, :sequence_file, :transcript_file_id, \
                                :reference_id, :prefix, :do_pseudo, \
                                :do_contiguate, :do_exonerate, :do_ratt, \
                                :use_transcriptome_data, :organism, \
                                :max_gene_length, :augustus_score_threshold , \
                                :abacas_match_size, :abacas_match_sim, \
                                :taxon_id, :db_id, :ratt_transfer_type, \
                                :no_resume, :email, sequence_file_attributes: [:id, :file],
                                transcript_file_attributes: [:id, :file])

  end
end
