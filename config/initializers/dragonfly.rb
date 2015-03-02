require 'dragonfly'

class GFF3Analyser
  def call(content)
    path = content.path
    rval = system("gt gff3 -retainids -tidy -show no #{path}")
#    puts "rval #{rval}"
    rval
  end
end

class SequenceValidAnalyser
  def call(content)
    path = content.path
    rval = system("gt convertseq -noseq #{path}")
#    puts "rval #{rval}"
    rval
  end
end

class SequenceNumberAnalyser
  def call(content)
    path = content.path
    rval = `gt seqstat -contigs #{path}`
    m = rval.match(/number of contigs:\s+(\d+)/)
    nof_contigs = 0
    if m then
      nof_contigs = m[1].to_i
    end
#     puts "nof_contigs #{nof_contigs}"
    return nof_contigs
  end
end

Dragonfly.app.configure do
  secret "b2909acc0dbe34c4a88d89c8c465f0f45b4a7b1d4e0a957b3181e6ac1063e2d3"

  url_format "/media/:job/:name"

  datastore :file,
    root_path: Rails.root.join('public/system/dragonfly', Rails.env),
    server_root: Rails.root.join('public')

# do not forget dos2unix!

  processor :canonicalize_gff3 do |content|
    content.shell_update do |old_path, new_path|
      "gt gff3 -retainids -fixregionboundaries -tidy -retainids #{old_path} > #{new_path}"
    end
  end

  analyser :valid_gff3, GFF3Analyser.new
  analyser :valid_sequence, SequenceValidAnalyser.new
  analyser :number_of_sequences, SequenceNumberAnalyser.new
end

# Logger
Dragonfly.logger = Rails.logger

# Mount as middleware
Rails.application.middleware.use Dragonfly::Middleware

# Add model functionality
if defined?(ActiveRecord::Base)
  ActiveRecord::Base.extend Dragonfly::Model
  ActiveRecord::Base.extend Dragonfly::Model::Validations
end