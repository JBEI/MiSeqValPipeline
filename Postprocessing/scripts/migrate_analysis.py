#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('conffile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('--destination')
parser.add_argument('--name')

args = parser.parse_args()

import os
import shutil
from subprocess import call, Popen

from lxml import etree

jobfiles = {
  'ccs':['GATK/callable.bed','GATK/aligned_reads.bam','GATK/aligned_reads.bam.bai','GATK/snps.gatk.vcf','GATK/snps.gatk.vcf.idx'],
  'roi':['GATK/callable.bed','GATK/aligned_reads.bam','GATK/aligned_reads.bam.bai','GATK/snps.gatk.vcf','GATK/snps.gatk.vcf.idx'],
  'sub':['data/coverage.bed','data/aligned_reads.bam','data/aligned_reads.bam.bai','data/variants.vcf', 'data/variants.gff.gz',],
}

def rsync_if_exists(file,dest):
  if os.path.exists(file):
    if os.path.isdir(file):
      if call(['rsync','-av',file, dest]) == 0:
        return True
    else:
      if call(['rsync','-t',file, dest]) == 0:
        return True
  return False

def mkdir_p(path):
  if not os.path.exists(path): os.mkdir(path)

if __name__=='__main__':
  omask = os.umask(002)
  fh = args.conffile
  parser = etree.XMLParser(remove_blank_text=True)
  tree = etree.parse(fh,parser)
  root = tree.getroot()
  assert root.tag == 'SBAnalysis', "ERROR: Root tag is not SBAnalysis"
  assert len(root.getchildren())==1, "ERROR: Only one analysis can be migrated"
  analysis = root[0]
  assert analysis.tag == 'analysis', "ERROR: incorrect analysis tag"
  
  try:
    #--- Setup analysis root directory  
    old_path = os.path.join(analysis.attrib['location'],analysis.attrib['name'])
    assert os.path.exists(old_path), "Original directory does not exist!"
    assert os.path.exists(args.destination), "Destination directory does not exist!"
    new_path = os.path.join(args.destination, args.name)
    mkdir_p(new_path)

    #--- Copy reference files
    old_refdir = os.path.join(old_path,"references")
    rsync_if_exists(old_refdir,new_path)
    new_refdir = os.path.join(new_path,"references")    

    for pool in analysis:
      old_pool_path = os.path.join(old_path,pool.attrib['name'])
      new_pool_path = os.path.join(new_path,pool.attrib['name'])
      mkdir_p(new_pool_path)
      for job in pool:
        old_job_path = os.path.join(old_pool_path,job.attrib['name'])
        new_job_path = os.path.join(new_pool_path,job.attrib['name'])
        mkdir_p(new_job_path)
        if not rsync_if_exists(os.path.join(old_job_path,'call_summary.txt'),new_job_path):
          print >>sys.stderr, "Not found: %s" % os.path.join(old_job_path,'call_summary.txt')
        for f in jobfiles[job.attrib['protocol']]:
          if not rsync_if_exists(os.path.join(old_job_path,f), new_job_path):
            print >>sys.stderr, "Not found: %s" % os.path.join(old_job_path,f)
    
    # Output configuration file for analysis
    analysis.attrib['reference'] = os.path.join(new_refdir,"seqs/sequence/seqs.fasta")
    analysis.attrib['name']     = args.name
    analysis.attrib['location'] = args.destination
    with open(os.path.join(new_path,'config.xml'),'w') as outh:
      print >>outh, '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'''
      print >>outh, '''<SBAnalysis>'''
      print >>outh, etree.tostring(analysis,pretty_print=True)
      print >>outh, '''</SBAnalysis>'''

  except (IOError,OSError) as e:
    if e.args[1] == 'Read-only file system':
      sys.exit('Destination "%s" is mounted as read-only.' % args.destination)
    else: raise
