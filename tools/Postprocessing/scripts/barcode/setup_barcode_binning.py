#! /usr/bin/env python
import sys
import argparse
from collections import defaultdict
import re
import subprocess

def get_pool_map(start_pool,npools):
  ''' Returns dict mapping pool numbers to sample names '''
  m = re.match('([a-zA-Z]+)([0-9]+)',start_pool)
  return dict((i+1,'%s%d' % (m.group(1),int(m.group(2))+i)) for i in range(npools))

def decodeList(s,delim=','):
  ''' Converts a string representing a list, possibly containing ranges, to python list'''
  rangeRE = re.compile('^(\d+)-(\d+)$')
  numRE   = re.compile('^(\d+)$')
  numbers = []
  for f in s.split(delim):
    if numRE.match(f):
      numbers.append(int(f))
    else:
      m = rangeRE.match(f)
      assert m is not None
      start = int(m.group(1))
      end   = int(m.group(2))
      if end < start:
        start,end = (end,start)
      numbers.extend(range(start,end+1))
  return sorted(set(numbers))

def main(parser):
  args = parser.parse_args()

  # Iterator for input file
  lines = (l.strip('\n').split('\t') for l in args.barcode_table)
  header = lines.next()
  
  # Get the barcodes contained in each pool
  pool_barcodes = defaultdict(set)
  for l in lines:
    for p in decodeList(l[header.index('Pools')]):
      pool_barcodes[p].add(int(l[header.index('Barcode')]))
  
  # Get the PacBio sample ID for each pool
  # Assume that pool1 has sample ID provided by command
  # Assume that other pools are sequentially numbered
  last_pool = max(int(_) for _ in pool_barcodes.keys())
  pool_sample = get_pool_map(args.sample_id,last_pool)
  
  # Generate the commands and submit
  for p in sorted(pool_barcodes.keys()):
    bclist = ':'.join(str(_) for _ in sorted(pool_barcodes[p]))
    cmd = ["qsub", 
           "-N","%s_binning%s" % (args.analysis_name,p),
           "-v", "samp=%s,bclist=%s" % (pool_sample[p],bclist),
           args.job_sh
          ]
    if args.no_qsub:
      print 'NOT SUBMITTED: %s' % ' '.join(cmd)
    else:
      subprocess.call(cmd)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(prog='setup_barcode_binning', description='Set up binning for barcode')
  parser.add_argument('--job_sh', default="/global/dna/projectdirs/RD/synbio/Repos/JGISB_analysis.git/jobs/process_barcoded_run.sh")
  parser.add_argument('--no_qsub',action='store_true')
  parser.add_argument('sample_id', help='Sample ID for pool 1')
  parser.add_argument('analysis_name', help='Name for analysis')
  parser.add_argument('barcode_table', nargs='?', type=argparse.FileType('rU'), help='table with barcode information', default=sys.stdin)  
  main(parser)
