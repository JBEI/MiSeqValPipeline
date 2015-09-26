def shell(command)
  puts command
  system(command)
end

task :default do
  cmd = 'docker exec -it `docker ps -q` bash'
  puts cmd
  system(cmd)
end

task :build do
  shell 'docker build --tag seqval .'
end

task :run do
  cmd = <<-CMD
  ([ `docker ps -aq | wc -l` -ne 0 ] && docker stop $(docker ps -aq)) &&
     docker run --tty=false \
                   --rm -p 8003:80 \
                   -v /home/oge/code/MiSeqValPipeline:/src \
                   -v /home/oge/code/MiSeqValPipeline/home:/root \
                   --env SMB_USERNAME="${SMB_USERNAME}" \
                   --env SMB_PASSWORD="${SMB_PASSWORD}" \
                   --privileged seqval
CMD
  puts cmd
  system(cmd)
end

task :clearoldvms do
  cmd = "docker ps -a | grep 'weeks ago' | awk '{print $1}' | xargs --no-run-if-empty docker rm"
  puts cmd
  system(cmd)
end

task :s => :setup
