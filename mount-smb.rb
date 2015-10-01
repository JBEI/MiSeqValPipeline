require 'fileutils'
require 'optparse'

MISEQ_SHARE = "//smb.jbei.org/miseq"

def create_mount_point(directory)
  FileUtils.mkdir_p(directory)
end

def mount(username, password, share, directory)
  "mount -t cifs -o username=#{username},password=#{password} \
   #{share} #{directory}"
end

def shell(command)
  puts command
  puts `#{command}`
end

options = {}

optparser = OptionParser.new do |opts|
  opts.on("-u", "--username USERNAME") do |u|
    options[:username] = u
  end

  opts.on("-p", "--password PASSWORD") do |p|
    options[:password] = p
  end
end

optparser.parse!
options[:directory] = ARGV.pop

create_mount_point(options[:directory])
shell mount(options[:username],
            options[:password], 
            MISEQ_SHARE,
            options[:directory])
