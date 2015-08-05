#! /usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--vcffile')
parser.add_argument('--covfile')
parser.add_argument('--reffile')
# parser.add_argument('--json', default="/dev/null")
parser.add_argument('summary', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

import json
from Bio import SeqIO
from postanalysis.covvars import parse_covdepth_gatk, find_variants
from postanalysis.variant import Variant
from postanalysis.genedesign import inputGD, runGD

class Reference:
  def __init__(self,name):
    self.name      = name
    self.pct_cov   = None
    self.mean_cov  = None
    self.variants  = []
    self.dips      = []

  def make_call(self,seq):
    if self.pct_cov < 1.0: return "Incomplete"
    if self.mean_cov < 5.0: return "Low coverage"
    if not self.variants and not self.dips: return "Pass"
    
    result = runGD(inputGD(self.name, self.variants+self.dips, seq))
    if not result['score'] == '-1': return "Fixable:%s,%s" % (result['primer1'],result['primer2'])
    if self.variants: return 'Errors'
    return 'Dips'
      
  def summary(self):
    return '%s\t%.1f\t%.1f\t%d\t%d' % (self.name, 
                                             self.pct_cov*100, self.mean_cov, 
                                             len(self.variants), len(self.dips), 
                                       )

''' Load references '''
seqs = [(s.id,s) for s in SeqIO.parse(args.reffile,'fasta')]
sdict = dict((ref,Reference(name=ref)) for ref,seq in seqs)

''' Analyze coverage '''
covdata = parse_covdepth_gatk(args.covfile)
for ref,seq in seqs:
  result = find_variants(covdata[ref],seq,ref,exclude_edges=True,exclude_overlaps=True)
  sdict[ref].pct_cov = result['pct_cov']
  sdict[ref].mean_cov = result['mean_cov']  
  if 'variants' in result:
    for v in result['variants']:
      v.caller = 'covvars'
      sdict[ref].dips.append(v)

''' Analyze variants '''
vlines = [l.strip('\n') for l in open(args.vcffile,'rU') if not l.startswith('#')]
for l in vlines:
  v = Variant.from_vcf(l)
  v.caller = 'gatk'
  sdict[v.chrom].variants.append(v)

''' Output summary information '''
print >>args.summary, 'ref\tpct_cov\tmean_cov\tnvars\tndips\tcall'
for ref,seq in seqs:
  call = sdict[ref].make_call(seq)
  print >>args.summary, '%s\t%s' % (sdict[ref].summary(),call)
