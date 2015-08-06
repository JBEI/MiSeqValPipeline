#!/bin/sh

usage() {
    echo "Usage: postprocessing.sh SEQVAL_FOLDER MAINLIBRARY_FOLDER REF_FASTA GAKT_DIR PYTHON_DIR
            SEQVAL_FOLDER=Path to Seq Validation folder
            MAINLIBRARY_FOLDER=Path to mainLibrary folder
            REF_FASTA=Path to references.fasta file
            GATK_DIR=Path to GATK
            PYTHON_DIR=Path to Python Post Processing Scripts"
}

if [ "$#" -ne 2 ]; then
    usage
    exit
fi

SEQUENCEVALIDATIONPIPELINE_FOLDER=$1

MAINLIBRARY_FOLDER=$2
REF_FASTA=$3

GATK_DIR=$4
PYTHON_DIR=$5

if ! [ -f ${REF_FASTA} ]; then
    echo "${REF_FASTA} does not exist."
    exit
fi

for POOL_NAME in `ls -d *_*`;
do
    if [ -f ${POOL_NAME}/bwa_dir/${POOL_NAME}.bam ]; then
        mv ${POOL_NAME}/bwa_dir/${POOL_NAME}.bam ${POOL_NAME}/bwa_dir/aligned_reads.bam
    fi    

    if [ -f ${POOL_NAME}/bwa_dir/${POOL_NAME}.bam.bai ]; then
        mv ${POOL_NAME}/bwa_dir/${POOL_NAME}.bam.bai ${POOL_NAME}/bwa_dir/aligned_reads.bam.bai
    fi

    echo "Running GATK Depth of Coverage..."
    # DEPTH OF COVERAGE
    java -jar ${GATK_DIR}/GenomeAnalysisTK.jar -T DepthOfCoverage --includeRefNSites -R ${REF_FASTA} -I ${POOL_NAME}/bwa_dir/aligned_reads.bam -o ${POOL_NAME}/bwa_dir/covdepth > /dev/null 2>&1
    echo "covdepth file generated"

    echo ""

    ## VCF
    echo "Running GATK Unified Genotyper..."
    java -jar ${GATK_DIR}/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${REF_FASTA} -I ${POOL_NAME}/bwa_dir/aligned_reads.bam -glm BOTH --max_deletion_fraction 0.55 -o ${POOL_NAME}/bwa_dir/snps.gatk.vcf -nt 8 > /dev/null 2>&1
    echo "snps.gatk.vcf file generated"

    echo ""

    ## CALLABLE
    echo "Running GATK Callable Loci..."
    java -jar ${GATK_DIR}/GenomeAnalysisTK.jar -T CallableLoci -R ${REF_FASTA} -I ${POOL_NAME}/bwa_dir/aligned_reads.bam -summary ${POOL_NAME}/call_summary.txt -o ${POOL_NAME}/bwa_dir/callable.bed > /dev/null 2>&1
    echo "callable.bed file generated"

    echo ""

    # CALL SUMMARY
    echo "Running make_calls_gatk.py script..."
    python ${PYTHON_DIR}/make_calls_gatk.py --reffile ${REF_FASTA} --vcffile ${POOL_NAME}/bwa_dir/snps.gatk.vcf --covfile ${POOL_NAME}/bwa_dir/covdepth > ${POOL_NAME}/bwa_dir/call_summary.txt
    echo "call_summary.txt file generated"

    echo ""
done


## generate HTML
 # Dependencies:
 #        lxml
 #        openpyxl 
 #        scipy 
 #        numpy
python ${PYTHON_DIR}/summarize_analysis.py --fspath ${SEQUENCEVALIDATIONPIPELINE_FOLDER} < ${MAINLIBRARY_FOLDER}/config.xml

##
## WEB
##

## copy files to web-documents
# WEB_DIR=${SEQUENCEVALIDATIONPIPELINE_FOLDER}/website

# if [ -d ${WEB_DIR}/${LIBRARY_NAME} ]; then
#     rm -Rf ${WEB_DIR}/${LIBRARY_NAME}
# fi
# mkdir ${WEB_DIR}/${LIBRARY_NAME}

#cp ${LIBRARY_DIR}/${LIBRARY_NAME}/index.html ${WEB_DIR}/${LIBRARY_NAME}/index.html 
#cp -r ${LIBRARY_DIR}/${LIBRARY_NAME}/ref ${WEB_DIR}/${LIBRARY_NAME}/ref
#cp -r ${LIBRARY_DIR}/${LIBRARY_NAME}/results ${WEB_DIR}/${LIBRARY_NAME}/results
#cp -r ${LIBRARY_DIR}/${LIBRARY_NAME}/config.xml ${WEB_DIR}/${LIBRARY_NAME}/config.xml


## for every library, copy the aligned_reads.bam, aligned_reads.bam.bai, call_summary.txt, callable.bed, snps.gatk.vcf, and snps.gatk.vcf.idx 
## into the $WEB_DIR

# for POOL_NAME in `ls -d *_*`;
# do
#     mkdir ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}
#     mkdir ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir

#     cp ${LIBRARY_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/aligned_reads.bam ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/aligned_reads.bam  
#     cp ${LIBRARY_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/aligned_reads.bam.bai ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/aligned_reads.bam.bai 

#     cp ${LIBRARY_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/call_summary.txt ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/call_summary.txt  
#     cp ${LIBRARY_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/callable.bed ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/callable.bed 
#     cp ${LIBRARY_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/snps.gatk.vcf ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/snps.gatk.vcf 
#     cp ${LIBRARY_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/snps.gatk.vcf.idx ${WEB_DIR}/${LIBRARY_NAME}/${POOL_NAME}/bwa_dir/snps.gatk.vcf.idx 
# done