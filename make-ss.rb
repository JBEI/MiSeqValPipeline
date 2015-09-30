require 'json'

def svelt_cmd(dir, clone, pool)
  "python ../MiSeqValPipeline/svelt.py MS #{dir} #{clone} #{pool} #{dir}/analysis"
end

def call_cmd(dir, call)
  "echo #{call} > #{dir}/call"
end

def score_cmd(dir, score)
  "echo #{score} > #{dir}/score"
end

def call(index, clone, pool)
  index[clone][pool]['call']
end

def score(index, clone, pool)
  index[clone][pool]['display']
end

dir = ARGV[0]
clone = ARGV[1]
pool = ARGV[2]
index = JSON.parse(File.read(ARGV[3]))

puts "#{svelt_cmd(dir, clone, pool)} && \
  #{call_cmd(dir, call(index, clone, pool))} && \
  #{score_cmd(dir, score(index, clone, pool))}"
