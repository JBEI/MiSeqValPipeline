def svelt_cmd(dir, clone, pool)
  "python ../MiSeqValPipeline/svelt.py MS #{dir} #{clone} #{pool} #{dir}/analysis"
end

def call(dir, call)
  "echo #{call} > #{dir}/call"
end

def score(dir, score)
  "echo #{score} > #{dir}/score"
end

dir = ARGV[0]
clone = ARGV[1]
pool = ARGV[2]

puts "#{svelt_cmd(ARGV[0], ARGV[1], ARGV[2])} && \
  #{call(dir, 'bae')} && \
  #{score(dir, 419)}"
