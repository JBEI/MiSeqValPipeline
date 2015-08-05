#! /usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--vcffile')
parser.add_argument('--reffile')
# parser.add_argument('--json', default="/dev/null")
parser.add_argument('summary', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

import json
from Bio import SeqIO
from postanalysis.variant import Variant

''' Load references '''
seqs = [(s.id,s) for s in SeqIO.parse(args.reffile,'fasta')]
sdict = dict(seqs)

''' '''
variants = {}
vlines = [l.strip('\n') for l in open(args.vcffile,'rU') if not l.startswith('#')]
for l in vlines:
  v = Variant.from_vcf(l)
  v.caller = 'gatk'
  if v.chrom not in variants:
    variants[v.chrom] = []
  variants[v.chrom].append(v)

''' Output summary information '''
print >>args.summary, 'reference\tvariants'
for ref,seq in seqs:
  if ref in variants:
    print >>args.summary, '%s\t%d' % (ref,len(variants[ref]))
  else:
    print >>args.summary, '%s\t0' % ref
