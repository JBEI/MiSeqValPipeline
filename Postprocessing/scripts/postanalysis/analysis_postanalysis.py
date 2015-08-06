#! /usr/bin/env python

import sys
import os

try:
  import cPickle as pickle
except:
  import pickle

from Bio import SeqIO

''' postanalysis imports '''
import variant,calls

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Run postanalysis for analysis.')
  parser.add_argument('--analysis_dir',default='.')
  parser.add_argument('--reffile',default='references.fasta')
  parser.add_argument('--base_url',default="http://synbio.jgi-psf.org/pbdata")
  parser.add_argument('--poolindex',default='poolindex.txt')
  parser.add_argument('--analysis_root',default='pacbioSB')
  args = parser.parse_args()

  if not os.path.isdir(args.analysis_dir): sys.exit('Error: directory "%s" does not exist' % args.analysis_dir)
  if not os.path.exists(args.reffile):     sys.exit('Error: reference file "%s" does not exist' % args.reffile)  
  if not os.path.exists(args.poolindex):   sys.exit('Error: pool index file "%s" does not exist' % args.poolindex)

  analysis_path  = os.path.abspath(args.analysis_dir)
  reference_file = os.path.abspath(args.reffile)
  poolindex_file = os.path.abspath(args.poolindex)
  
  # load sequences
  seqs = dict( [(s.id,s) for s in SeqIO.parse(reference_file,'fasta')] )
  
  # load job index
  poolIDs = [l.strip() for l in open(poolindex_file,'rU')]

  print >>sys.stderr, '[ Loading pool summaries ]'
  pooldata = {}
  all_ready = True
  for poolname in poolIDs:
    pool_pickle = '%s/%s/pool_summary.pickle' % (analysis_path,poolname)
    if not os.path.exists(pool_pickle):
      all_ready = False
      print >>sys.stderr, '%s is not ready: file "%s" does not exist' % (poolname,pool_pickle)
    else:
      with open(pool_pickle,'rb') as fh:
        pooldata[poolname] = pickle.load(fh)

  if not all_ready:
    sys.exit(1)

  print >>sys.stderr, '\t%s' % '\t'.join([pname for pname in sorted(pooldata.keys())])
  for ref in sorted(seqs.keys()):
    print >>sys.stderr, '%s\t%s' % (ref,    '\t'.join( [pooldata[pname][ref].call for pname in sorted(pooldata.keys())] ) )

  sys.exit(0)