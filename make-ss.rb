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

index = JSON.parse(File.read(ARGV[0]))

index.keys.each do |clone|
  index[clone].keys.each do |pool|
    puts "#{svelt_cmd(pool, clone, pool)} && \
      #{call_cmd(pool, call(index, clone, pool))} && \
      #{score_cmd(pool, score(index, clone, pool))}"
  end
end
