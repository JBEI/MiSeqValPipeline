#!/usr/bin/env python

import os

IGV_XML = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="%(urlpath)s/refs/%(ref)s.fasta" locus="All" version="5">
  <Resources>
    <Resource path="%(urlpath)s/%(pool)s/%(bc)s/ccs/%(ref)s.callable.bed" name="ccs_callable"/>  
    <Resource path="%(urlpath)s/%(pool)s/%(bc)s/ccs/%(ref)s.bam" name="ccs_reads"/>
    <Resource path="%(urlpath)s/%(pool)s/%(bc)s/ccs/%(ref)s.vcf" name="ccs_variants"/>    
  </Resources>
</Session>

"""
def get_summary_paths(base,atype):
  from glob import glob
  globstr = '%s/pool?/B???/%s/*.summary.txt' % (base,atype)
  return glob(globstr)

def main(parser):
  args = parser.parse_args()

  aroot = os.path.join(args.analysis_dir, args.analysis_name)
  urlpath = '%s/%s' % (args.url_path,args.analysis_name)
  summaries = get_summary_paths(aroot,args.analysis_type)
  for s in summaries:
    fields = s.replace(aroot,'').split('/')[1:]
    ref  = fields[3].split('.')[0]
    bc   = fields[1]
    pool = fields[0]
    with open('%s.%s.%s.xml' % (ref,bc,pool.strip('pool')),'w') as outh:
      print >>outh, IGV_XML % {'urlpath':urlpath,'pool':pool,'bc':bc,'ref':ref}

if __name__=='__main__':
  import argparse
  import sys

  parser = argparse.ArgumentParser()
  parser.add_argument('--analysis_type',default="ccs")
  parser.add_argument('--analysis_dir',default="/global/dna/shared/synbio/PacBioAnalysis")
  parser.add_argument('--url_path',default="http://synbio.jgi-psf.org:8077/analysis/")
  parser.add_argument('analysis_name')
  main(parser)
