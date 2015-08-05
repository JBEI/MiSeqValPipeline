GP_HEADER = '''#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=8:00:00
#$ -P gentech-rnd.p
#$ -pe pe_slots 8
#$ -w e
#$ -j y
#$ -R y

export _JAVA_OPTIONS="-Xmx12G"
umask 0002
'''

XFER_HEADER = '''#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=1:00:00
#$ -P gentech-rnd.p
#$ -w e
#$ -j y
#$ -R y
#$ -q xfer.q

umask 0002

module load python
module use /global/dna/projectdirs/RD/synbio/Modules
module load jgisb_analysis

SBA_migrate --destination %(archive_path)s --name %(analysis_name)s < config.xml
SBA_summarize --fspath %(archive_path)s --urlpath %(url_path)s < %(archive_path)s/%(analysis_name)s/config.xml

echo -e "Results from your analysis are ready to be viewed:\\n%(url_path)s/%(analysis_name)s/' mail -c mlbendall@lbl.gov -s "[JGISB_Analysis] Analysis results for %(analysis_name)s" %(notify)s
'''

GP_HEADER_MENDEL = '''#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=4:00:00
#$ -l exclusive.c
#$ -P gentech-rnd.p
#$ -w e
#$ -j y
#$ -R y

export _JAVA_OPTIONS="-Xmx12G"
umask 0002
'''

SMRT_BLOCK = '''
module load %(smrtmodule)s
smrtpipe.py %(options)s --params=%(settings_xml)s --output=%(outdir)s xml:%(input_xml)s &> smrtpipe.err
module unload %(smrtmodule)s
'''

GATK_BLOCK = '''
# A bug in GATK 2.8-1 produces this error: Null alleles are not supported
# The bug is supposedly fixed in later releases, not currently available on Genepool (3/3/2014)
# module load GATK/2.8-1
module load %(gatkmodule)s
module load picard/1.92
module load samtools/20130508

mkdir -p GATK

pbbam=data/aligned_reads.bam
ref=../../references/seqs/sequence/seqs.fasta
unsorted=GATK/unsorted.bam
bamfile=GATK/aligned_reads.bam

# Adjust base quality scores to work with GATK
# CCS and ROI can generate very high base qualities, sometimes > 80. GATK detects the high
# scores and assumes the encoding is incorrect.
# To work around this, simply change any quality score >41 to 41. (Any character >"K" -> "J")
samtools view -h $pbbam | perl -n -e '@c=split(/\\t/); $c[10]=~s/[K-~]/J/g unless /^@/; print join("\\t",@c)' > $unsorted.sam

# Reorder BAM file so it matches reference file
picard ReorderSam I=$unsorted.sam R=$ref O=$unsorted

# Sort and index BAM
picard SortSam I=$unsorted O=$bamfile SORT_ORDER=coordinate
picard BuildBamIndex I=$bamfile O=$bamfile.bai
rm $unsorted.sam $unsorted

[ ! -e "$bamfile" ] && echo "ERROR: bamfile $bamfile does not exist" && exit 1
[ ! -e "$ref" ] && echo "ERROR: reference file $ref does not exist" && exit 1

# FindCoveredIntervals
# Output intervals that are "uncovered" with depth < %(mincov)s
GenomeAnalysisTK -T FindCoveredIntervals -R $ref -I $bamfile -cov %(mincov)s -u -o GATK/uncovered.txt

# Callable loci
GenomeAnalysisTK -T CallableLoci -R $ref -I $bamfile -summary GATK/callable.summary -o GATK/callable.bed

# Generate coverage at every base
GenomeAnalysisTK -T DepthOfCoverage -R $ref -I $bamfile -o GATK/covdepth
### python $SBPATH/analyze_coverage.py --covfile GATK/covdepth --reffile $ref > GATK/coverage_summary.txt

# Call variants using Unified Genotyper
GenomeAnalysisTK -T UnifiedGenotyper -R $ref -I $bamfile -glm BOTH --max_deletion_fraction 0.55 -o GATK/snps.gatk.vcf -nt 8
### python $SBPATH/analyze_variants.py --vcffile GATK/snps.gatk.vcf --reffile $ref > GATK/variants_summary.txt

# Call variants using Haplotype caller
# This is not multi-threaded, so runs slowly
# GenomeAnalysisTK -T HaplotypeCaller -R $ref -I $bamfile -o GATK/snps.gatkHC.vcf
'''

GATK_POST_BLOCK = '''
# Postanalysis for GATK
module load dnassemble/1.0
module use /global/dna/projectdirs/RD/synbio/Modules
module load jgisb_analysis

ref=../../references/seqs/sequence/seqs.fasta
SBA_calls_gatk --reffile $ref --vcffile GATK/snps.gatk.vcf --covfile GATK/covdepth > call_summary.txt
'''

SUB_POST_BLOCK = '''
# Postanalysis for PacBio consensus calls
module load samtools/20130508

pbbam=data/aligned_reads.bam
ref=../../references/seqs/sequence/seqs.fasta

# Calculate depth using samtools
samtools depth $pbbam > data/covdepth

# Postanalysis for gencons
module load dnassemble/1.0
module use /global/dna/projectdirs/RD/synbio/Modules
module load jgisb_analysis
SBA_calls_gencons --reffile $ref --gfffile data/variants.gff.gz --covfile data/covdepth > call_summary.txt
'''

std_opts = {'smrtmodule':'smrtanalysis/2.1.1',
            'options':'',
            'input_xml':'../input.xml',
            'settings_xml': 'settings.xml',
            'outdir': '.',
            'jobname':'test',
            'gatkmodule':'GATK/2.7-2',
            'mincov':'5'
           }


def generate_jobsh(protocol,**kwargs):
  opts = dict(std_opts)
  for k,v in kwargs.iteritems():
    opts[k] = v
  shstr =  GP_HEADER % opts
  shstr += SMRT_BLOCK % opts
  if protocol in ['ccs','roi']:
    shstr += GATK_BLOCK % opts
    shstr += GATK_POST_BLOCK
  elif protocol == 'sub':
    shstr += SUB_POST_BLOCK
  return shstr

sum_opts = {'url_path':'http://synbio.jgi-psf.org:8077/analysis',
           }
#            'archive_path':'',
#            'analysis_name':'',
#            'notify':'',            
           
def generate_summarysh(**kwargs):
  opts = dict(sum_opts)
  for k,v in kwargs.iteritems():
    opts[k] = v

  shstr = XFER_HEADER % opts
  return shstr
