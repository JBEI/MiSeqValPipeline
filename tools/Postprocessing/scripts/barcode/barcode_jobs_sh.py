
ALIGN_CCS = """#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=8:00:00
#$ -P gentech-rnd.p
#$ -l exclusive.c
#$ -w e
#$ -j y
#$ -R y

##########################################################################################
# Aligns barcoded reads to reference and calls variants
##########################################################################################
#--- Run options
# PacBio sample name, i.e. HAR33306
[[ -z ${samp} ]] && echo "No sample!" && exit 1
# List of barcode numbers. This should be a colon-delimited list; had difficulty passing 
# a comma-delimited list through qsub -v
[[ -z ${bclist} ]] && echo "No barcodes!" && exit 1
# Temporary directory
[[ -z ${dest} ]] && dest=`mktemp -d`
# Final destination (on DnA)
dnadest="/global/dna/projectdirs/RD/synbio/Runs/$samp"



##########################################################################################
# Initialize run directory
##########################################################################################
umask 0002
export _JAVA_OPTIONS="-Xmx12G"

module load smrtanalysis/2.2.0

RUNDIR=/global/dna/projectdirs/RD/synbio/Runs
SAMP=HAR33306
BC=65

mkdir -p ccs
ccsreads=$RUNDIR/$SAMP/ccs/bins/B$(printf "%03d" $BC).fasta
ccsfastq=$RUNDIR/$SAMP/ccs/bins/B$(printf "%03d" $BC).fastq
## ccspulse=$RUNDIR/$SAMP/ccs.fofn
ref=../ref.fasta
alnreads=ccs/aligned_reads

##########################################################################################
# Trim reads to remove barcodes
##########################################################################################
[[ -z ${bctrim} ]] && bctrim=20

SBA_trimmer --trim $bctrim < $ccsfastq > trimmed.fastq
SBA_fastq_to_fasta < trimmed.fastq > trimmed.fasta

##########################################################################################
# Align CCS reads using pbalign.py
##########################################################################################

#--- Do not load pulse data into cmpH5
pbalign.py \
  --nproc 8 \
  --minAccuracy 0.75 \
  --minLength 500 \
  --maxHits 1 \
  --hitPolicy randombest \
  trimmed.fasta \
  $ref \
  ${alnreads}.orig.sam

SBA_add_qual_to_sam $ccsfastq ${alnreads}.orig.sam | \
  perl -n -e '@c=split(/\t/); $c[10]=~s/[K-~]/J/g unless /^@/; print join("\t",@c)' | \
  perl -n -e '@c=split(/\t/); $c[8]="0" unless /^@/; print join("\t",@c)' > ${alnreads}.qual.sam

module load picard/1.92

# Reorder BAM file so it matches reference file
picard ReorderSam I=${alnreads}.qual.sam R=$ref O=${alnreads}.unsorted.bam

# Sort and index BAM
picard SortSam I=${alnreads}.unsorted.bam O=${alnreads}.bam SORT_ORDER=coordinate
picard BuildBamIndex I=${alnreads}.bam O=${alnreads}.bam.bai

rm ${alnreads}.qual.* ${alnreads}.unsorted.*
bamfile=${alnreads}.bam


##########################################################################################
# Run GATK
##########################################################################################
module load GATK/2.7-2

# FindCoveredIntervals
# Output intervals that are "uncovered" with depth < 5
GenomeAnalysisTK -T FindCoveredIntervals -R $ref -I $bamfile -cov 5 -u -o ccs/uncovered.txt

# Callable loci
GenomeAnalysisTK -T CallableLoci -R $ref -I $bamfile -summary ccs/callable.summary -o ccs/callable.bed

# Generate coverage at every base
GenomeAnalysisTK -T DepthOfCoverage -R $ref -I $bamfile -o ccs/covdepth

# Call variants using Unified Genotyper
GenomeAnalysisTK -T UnifiedGenotyper -R $ref -I $bamfile -glm BOTH --max_deletion_fraction 0.55 -o ccs/snps.gatk.vcf -nt 8

"""