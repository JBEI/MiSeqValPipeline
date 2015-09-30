__author__ = 'simi'
#!/usr/bin/env python
import csv
import glob
import sys
import os
from Bio import SeqIO


def create_fasta(type, path, out_name, constructs):
    print("creating new fasta file...")
    # parse out the sequences needed for the refrerence fasta
    # create new fasta file
    new_fasta = open(out_name+".fasta", "w+")
    old_fasta_path = get_fasta_path(type, path)
    # get sequences from old fasta
    handle = open(old_fasta_path, "rU")
    sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    # loop through constructs
    for con in constructs:
        clone = con['clone']
        print("clone: "+clone)
        s = sequences[clone]
        seq = ''.join(s.seq)
        # write to new fasta
        # print("seq: "+seq)
        new_fasta.write(">"+clone+"\n"+seq+"\n")
    new_fasta.close()
    print("Done creating new reference fasta file")


def create_bam(type, path, out_name, constructs):
    print("creating new bam and bai files...")
    # create the shell script that will create the mini bam files,
    # merge them and then create the bai file
    count=0
    merge_names =[]
    # loop through constructs, and make mini bam files - one for each
    for con in constructs:
        count+=1
        pool = con['pool']
        clone = con['clone']
        print("clone: "+clone+", pool: "+pool)
        bam = get_bam_path(type, path, pool)
        mini_bam = "in"+str(count)+".bam"
        # add mini_bam name to the merge_names for the merge command
        merge_names.append(mini_bam)
        # create command string for creating the mini_bam
        cmd = 'samtools view -h -b '+bam+' "'+clone+'" > '+mini_bam
        # run the command
        os.system(cmd)


    if len(constructs) == 1:
       merge_cmd = "mv in1.bam "+out_name+".bam"
    else:
        merge_cmd = "samtools merge "+out_name+".bam "+" ".join(merge_names)

    print("merge command: "+merge_cmd)
    os.system(merge_cmd)
    # now create the bai - index
    index_cmd = "samtools index -b "+out_name+".bam"
    print("index command: "+index_cmd)
    os.system(index_cmd)
    print("Done creating new bam and bai files")


def create_vcf(type, path, outname, constructs):
    print("creating new vcf file...")
    # create a new vcf
    new_vcf = open(outname+".vcf", "w+")
    # so keep track of the pool files that are open - only works if clones in csv file are in same order as the analysis
    open_pool_files = {}
    #the first one opened will be the one we copy the header from, and get the contig ID lines for each construct
    first = True
    for con in constructs:
        pool = con['pool']
        clone = con['clone']
        print("clone: "+clone+", pool: "+pool)
        # see if pool file is already open
        if not first and  (pool in open_pool_files):
            pool_file = open_pool_files[pool]
        else:
            # get pool file
            pool_file = read_vcf(type, path, pool)
            # add pool file to open pool files dict
            open_pool_files[pool]= pool_file
        for line in iter(pool_file):
            if first and line.startswith("##"):
                if line.startswith("##contig"):
                    first = False
                else:
                    new_vcf.write(line)
            if line.startswith("##contig=<ID="+clone):
                new_vcf.write(line)
                break
    # now we will add our ##reference line and #Chrom line -
    # TODO make sure that these are relatively generic for our uses
    new_vcf.write("##reference=file://"+path+"/"+outname+".fasta\n")
    new_vcf.write("#CHROM   POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  c100843992550000001823187012311500\n")
    # now loop through the constructs list again, this time adding any line that starts with the clone name
    # also, replace any white space in the line with a tab
    print("second loop through constructs...")
    for con in constructs:
        pool = con['pool']
        clone = con['clone']
        print("clone: "+clone+", pool: "+pool)
        # the pool file will already be open from the first step, but will need to reset the cursor
        pool_file = open_pool_files[pool]
        #reset read start
        pool_file.seek(0)
        for line in iter(pool_file):
            if line.startswith(clone):
                # replace white spaces with tabs
                '\t'.join(line.split())
                new_vcf.write(line)
                print("writing: "+line)
    # now close all the files
    # close all the open vcf files
    for f in open_pool_files.values():
        f.close()
    #close new bed
    new_vcf.close()
    print("Done creating new vcf file")


def create_bed(type, path, outname, constructs):
    # parse the lines needed from each bed, for the correct pool for each construct
    # need to create the new file for writing
    print("creating new bed file...")
    new_bed = open(outname+".bed", "w+")
    # for each construct in constructs, open pool bed file, and copy the needed line
    # it's more efficient to copy all constructs from each pool...
    # so keep track of the pool files that are open - only works if clones in csv file are in same order as the analysis
    open_pool_files = {}
    for con in constructs:
        pool = con['pool']
        clone = con['clone']
        print("clone: "+clone+", pool: "+pool)
        # see if pool file is already open
        if (pool in open_pool_files):
            pool_file = open_pool_files[pool]
        else:
            # get pool file
            pool_file = read_bed(type, path, pool);
            # add pool file to open pool files dict
            open_pool_files[pool]= pool_file
        # find our clone line in bed
        for line in iter(pool_file):
            if line.startswith(clone):
                print("writing: "+line)
                # printlin in new bed
                new_bed.write(line)
                break
    # close all the open bed files
    for f in open_pool_files.values():
        f.close()
    # close new bed
    new_bed.close()
    print("Done creating new bed file")


##### utility functions ######

def get_fasta_path(type, path):
    # for PB, fasta is in: path/references/seqs/sequence/seqs.fasta
    print("getting old fasta path")
    if (type == "PB"):
         f_path = path+"/references/seqs/sequence/seqs.fasta"
         return f_path
    elif (type == "MS"):
        f_path = glob.glob("ref/*.fasta")[0]
        return f_path
    else:
         print("Type: "+type+" is not currently supported.")


def get_bam_path(type, path, pool):
     print("getting bam path for pool: "+pool)
     if (type == "PB"):
         # for pac bio analyses, this file <path>/<pool>/roi/aligned_reads.bam
         bam_path = path+"/"+pool+"/roi/aligned_reads.bam"
         return bam_path
     elif (type == "MS"):
        bam_path = path+"/bwa_dir/aligned_reads.bam"
        return bam_path
     else:
         print("Type: "+type+" is not currently supported.")


def read_vcf(type, path, pool):
    '''
    Opens the vcf file - similar to red_bed
    :param type:
    :param path:
    :param pool:
    :return:
    '''
    print("getting vcf file for pool: "+pool)
    if (type == "PB"):
        # for pac bio analyses, this file <path>/roi/<pool>/callable.bed
        vcf_file = open(path+"/"+pool+"/roi/snps.gatk.vcf", "rU")
        return vcf_file
    elif (type == "MS"):
        vcf_file = open(path+"/bwa_dir/snps.gatk.vcf", "rU")
        return vcf_file
    else:
        print("Type: "+type+" is not currently supported.")

def read_bed(type, path, pool):
    '''
    Opens the bed file for reading and returns the file handel
    :param type: PB \ MS, for pac bio and miseq data may have different file structures
    :param path: directory where files are kept for this analysis
    :return: the bed file handel for reading
    '''
    print("getting bed file for pool: "+pool)
    if (type == "PB"):
        # for pac bio analyses, this file <path>/roi/<pool>/callable.bed
        bed_file = open(path+"/"+pool+"/roi/callable.bed", "rU")
        return bed_file
    elif (type == "MS"):
        bed_file = open(path+"/bwa_dir/callable.bed", "rU")
        return bed_file
    else :
         print("Type: "+type+" is not currently supported.")



def get_constructs(filename):
    print("parsing csv file...")
    # parses csv and returns a dict, keys are the construct names, value is the pool
    constructs = []
    # open csv file filename
    with open(filename, 'rU') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            clone = row[0]
            pool = row[1]
            # ignore header row if exists
            if (pool != "pool"):
                # put in dict and list
                # TODO make this able to acept a 3rd column - for users name,
                # so we can swap out our name for the users name
                constructs.append({'clone':clone, 'pool':pool})
                print("put in list: clone: "+clone+", pool: "+pool)
    return constructs


def set_up(type, path, filename, out_name):
    # parse the vcf and bed files, then create the shell script to create the bam and bai files

    print("Setting up svelt for "+out_name)
    constructs = get_constructs(filename)
    if len(constructs) >0:
        create_bed(type, path, out_name, constructs)
        create_vcf(type, path, out_name, constructs)
        create_bam(type, path, out_name, constructs)
        create_fasta(type, path, out_name, constructs)
        print("Done running svelt for "+out_name)
    else:
        print("There were no clones specified in your .csv file")


def main():
    '''
    upload directory path for analysis, type [ MS | PB], file.csv name, and the output name for the files
    Assumes you are running this in the same directory as the .csv file, and al the new files will be created in the same directory
    Right now only supports PB (the MiSeq directory structure may be different
    It would be nice to modify this to substitute the users name for constructs with our anominized name
    need to use modules for biopython and smatools
    '''
    if len(sys.argv) < 5:
        print("Usage: python svelt.py PB path/to/analysis/directory file.csv output_name")
    else:
        type = sys.argv[1]
        path = sys.argv[2]
        name = sys.argv[3]
        out_name = sys.argv[4]
        if (type != "PB") and (type != "MS"):
            print("Only PacBio (PB) or MiSeq (MS) is currently supported")
        else:
            set_up(type, path, name, out_name)


if __name__ == "__main__":
    main()
