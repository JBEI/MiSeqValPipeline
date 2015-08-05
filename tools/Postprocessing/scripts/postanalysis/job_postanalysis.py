#! /usr/bin/env python

import sys
import os
import gzip

try:
  import cPickle as pickle
except:
  import pickle

from Bio import SeqIO

''' postanalysis imports '''
from covvars import *
from variant import Variant

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Run postanalysis for job.')
  parser.add_argument('--jobdir',default='.')
  parser.add_argument('--reffile',default='../../references.fasta')
  args = parser.parse_args()

  if not os.path.isdir(args.jobdir):   sys.exit('Error: directory "%s" does not exist' % args.pooldir)
  if not os.path.exists(args.reffile): sys.exit('Error: reference file "%s" does not exist' % args.reffile)  

  job_path       = os.path.abspath(args.jobdir)
  reference_file = os.path.abspath(args.reffile)

  # load sequences
  seqs = dict( [(s.id,s) for s in SeqIO.parse(reference_file,'fasta')] )

  summaries = dict(( (name,{'variants':[]}) for name in seqs.keys()))
  ''' GATK variants '''
  print >>sys.stderr, "[ Reading GATK variants ]"
  vlines = [l.strip('\n') for l in open('%s/GATK/snps.gatk.vcf' % job_path,'rU') if not l.startswith('#')]
  for l in vlines:
    v = Variant.from_vcf(l)
    v.caller = 'gatk'
    summaries[v.chrom]['variants'].append(v)

  ''' PacBio variants '''
  print >>sys.stderr, "[ Reading GenCons variants ]"
  glines = [l.strip('\n') for l in gzip.open('%s/data/variants.gff.gz' % job_path,'rb') if not l.startswith('#')]
  for l in glines:
    v = Variant.from_gff(l)
    v.caller = 'gencons'
    summaries[v.chrom]['variants'].append(v)

  ''' coverage variants '''
  print >>sys.stderr, "[ Reading coverage variants ]"
  covdata = parse_covdepth('%s/GATK/covdepth' % job_path)
  covvars = {}
  for ref in covdata.keys():
    assert ref in summaries, "Error: ref %s is not in summaries" % ref
    result = find_variants(covdata[ref],seqs[ref],ref,exclude_edges=True,exclude_overlaps=True)
    summaries[ref]['mean_cov'] = result['mean_cov']
    summaries[ref]['pct_cov'] = result['pct_cov']
    if 'variants' in result:
      covvars[ref] = {'variants':[]}
      for v in result['variants']:
        # v = result['variants'].pop()
        v.caller = 'covvars'
        summaries[ref]['variants'].append(v)
        covvars[ref]['variants'].append(v)

  print >>sys.stderr, "[ Writing results ]"
  with open('%s/job_summary.pickle' % job_path,'wb') as outh:
    pickle.dump(summaries,outh)

  sys.exit(0)
