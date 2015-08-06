#! /usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument('--gfffile')
parser.add_argument('--covfile')
parser.add_argument('--reffile')

parser.add_argument('summary', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

import json
import gzip
from Bio import SeqIO
from postanalysis.covvars import parse_covdepth_samtools, summarize_coverage, find_nocov_variants
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

''' Summarize coverage data '''
covdata = parse_covdepth_samtools(args.covfile,reflens=dict((s[0],len(s[1])) for s in seqs))
for ref,seq in seqs:
  p,m = summarize_coverage(covdata[ref])
  sdict[ref].pct_cov = p
  sdict[ref].mean_cov = m
  ncvars = find_nocov_variants(covdata[ref],chrom=ref,caller='samdepth')
  if ncvars is not None: sdict[ref].dips.extend(ncvars)

''' Analyze variants '''  
glines = [l.strip('\n') for l in gzip.open(args.gfffile,'rb') if not l.startswith('#')]
for l in glines:
    v = Variant.from_gff(l)
    v.caller = 'gencons'
    sdict[v.chrom].variants.append(v)

''' Output summary information '''
print >>args.summary, 'ref\tpct_cov\tmean_cov\tnvars\tndips\tcall'
for ref,seq in seqs:
  call = sdict[ref].make_call(seq)
  print >>args.summary, '%s\t%s' % (sdict[ref].summary(),call)
