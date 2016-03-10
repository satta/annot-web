class Job < ActiveRecord::Base
    belongs_to :user
    has_one :sequence_file, dependent: :destroy
    has_one :transcript_file, dependent: :destroy
    has_one :reference
    has_one :genome_stat, dependent: :destroy
    has_many :circos_images, dependent: :destroy
    has_many :result_files, dependent: :destroy
    has_many :genes, dependent: :delete_all
    has_many :clusters, dependent: :delete_all
    has_one :tree, dependent: :destroy
    validates_format_of :email, :with => /\A\z|\A([^@\s]+)@((?:[-a-z0-9]+\.)+[a-z]{2,})\z/i
    validates_presence_of :sequence_file
    accepts_nested_attributes_for :sequence_file, allow_destroy: true
    accepts_nested_attributes_for :transcript_file, allow_destroy: true

    def job_directory
      if not self[:job_id] then
        raise "no job ID yet"
      end
      if self[:job_id].length > 1 then
        Rails.root.join('public', 'jobs', self[:job_id])
      end
    end

    def temp_directory
      if not self[:job_id] then
        raise "no job ID yet"
      end
      if not CONFIG['tmpdir'] then
        raise "CONFIG['tmpdir'] not set"
      end
      "#{CONFIG['tmpdir']}/#{self[:job_id]}"
    end

    def work_directory
      if not self[:job_id] then
        raise "no job ID yet"
      end
      if not CONFIG['workdir'] then
        raise "CONFIG['workdir'] not set"
      end
      "#{CONFIG['workdir']}/#{self[:job_id]}"
    end

    def make_directories!
      if not Dir.exist?(self.job_directory) then
        Dir.mkdir(self.job_directory)
      end
      if not Dir.exist?(self.temp_directory) then
        Dir.mkdir(self.temp_directory)
      end
      if not Dir.exist?(self.work_directory) then
        Dir.mkdir(self.work_directory)
      end
    end

    def cleanup_work_directories!
      if not CONFIG['keep_work_directories'] then
        if File.exist?(self.temp_directory) then
          FileUtils.rm_rf(self.temp_directory)
        end
        if File.exist?(self.work_directory) then
          FileUtils.rm_rf(self.work_directory)
        end
      end
    end

    def cleanup_job_directory!
      if File.exist?(self.job_directory) then
        FileUtils.rm_rf(self.job_directory)
      end
    end
end
