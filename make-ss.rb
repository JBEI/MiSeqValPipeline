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

def zip(dir, clone, pool)
  intermediates = "analysis.bam analysis.bam.bai \
      analysis.vcf analysis.bed analysis.fasta call score"

  "cd #{dir} && \
    zip #{clone}.ss.zip #{intermediates} && \
    rm #{intermediates} && \
    cd .."
end

index = JSON.parse(File.read(ARGV[0]))

index.keys.each do |clone|
  index[clone].keys.each do |pool|
    cmd = "#{svelt_cmd(pool, clone, pool)} && \
      #{call_cmd(pool, call(index, clone, pool))} && \
      #{score_cmd(pool, score(index, clone, pool))} && \
      #{zip(pool, clone, pool)}"
    puts cmd
    `#{cmd}`
  end
end
