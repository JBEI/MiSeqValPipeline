task :default do
  cmd = 'docker exec -it `docker ps -ql` bash'
  puts cmd
  system(cmd)
end

task :setup do
  cmd = "docker build --tag seqval . && docker run --tty=false --rm --privileged seqval"
  puts cmd
  system(cmd)
end

task :clearoldvms do
  cmd = "docker ps -a | grep 'weeks ago' | awk '{print $1}' | xargs --no-run-if-empty docker rm"
  puts cmd
  system(cmd)
end
