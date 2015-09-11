task :default do
  cmd = 'docker exec -it `docker ps -ql` bash'
  puts cmd
  system(cmd)
end

task :setup do
  cmd = "docker stop $(docker ps -aq) && docker build --tag seqval . && docker run --tty=false --rm -p 8003:80 -v /home/oge/code/MiSeqValPipeline:/src -v /home/oge/code/MiSeqValPipeline/home:/root  --privileged seqval"
  puts cmd
  system(cmd)
end

task :clearoldvms do
  cmd = "docker ps -a | grep 'weeks ago' | awk '{print $1}' | xargs --no-run-if-empty docker rm"
  puts cmd
  system(cmd)
end

task :s => :setup
