#! /usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--covfile')
parser.add_argument('--reffile')
# parser.add_argument('--json', default="/dev/null")
parser.add_argument('summary', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

import json
from Bio import SeqIO
from postanalysis.covvars import parse_covdepth, find_variants

''' Load references '''
seqs = [(s.id,s) for s in SeqIO.parse(args.reffile,'fasta')]
sdict = dict(seqs)

''' '''
covout = {}
covdata = parse_covdepth(args.covfile)
for ref,seq in seqs:
  covout[ref] = {}
  result = find_variants(covdata[ref],seq,ref,exclude_edges=True,exclude_overlaps=True)
  covout[ref]['pct_cov']  = result['pct_cov']
  covout[ref]['mean_cov'] = result['mean_cov']
  #if 'variants' in result:
  #  covout[ref]['variants'] = []
  #  for v in result['variants']:
  #    v.caller = 'covvars'
  #    covout[ref]['variants'].append(v.serializable())

'''
with open(args.json,'w') as outh:
  print >>outh, json.dumps(covout)
'''

''' Output summary information '''
print >>args.summary, 'reference\tpct_cov\tmean_cov'
for ref,seq in seqs:
  print >>args.summary, '%s\t%.1f\t%.1f' % (ref,covout[ref]['pct_cov']*100,covout[ref]['mean_cov'])
