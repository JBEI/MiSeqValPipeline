require 'optparse'

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: example.rb [options]"

  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    options[:verbose] = v
  end

  opts.on("-h", "--host H", "Server to send request to") do |h|
    options[:host] = h
  end

  opts.on("-e", "--entryid E", "ID of the entry to upload to") do |e|
    options[:entryid] = e
  end


  opts.on("-s", "--sessionid S", "Login Session Id") do |s|
    options[:sessionid] = s
  end
  opts.on("-f", "--file F", "File to upload") do |f|
    options[:file] = f
  end
end.parse!

def shell(command)
  puts command
  `#{command}`
end

shell "curl --form \"file=@#{options[:file]}\" \
  -kv -H 'X-ICE-Authentication-SessionId: #{options[:sessionid]}' \
  #{options[:host]}/rest/parts/#{options[:entryid]}/shotgunsequences"
