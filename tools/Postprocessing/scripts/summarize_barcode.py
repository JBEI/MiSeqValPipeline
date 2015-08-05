#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--fspath')
parser.add_argument('--analysis_type')
parser.add_argument('--urlpath')

args = parser.parse_args()

import os
from glob import glob
from collections import defaultdict

from summarize.indexhtml import make_index
from summarize.merge import parse_callsum, call_display
# from summarize.merge import merge_calls, best_calls


# fspath = args.fspath
fspath = "/global/dna/projectdirs/RD/synbio/Analysis/BC082014"
urlpath = "http://synbio.jgi-psf.org:8180/analysis/BC082014"

# analysis_type = args.analysis_type
analysis_type = 'ccs'

prefix_len = len(fspath.split('/'))
summaries = glob('%(fspath)s/pool?/B???/%(analysis_type)s/*.summary.txt' % {'fspath':fspath,'analysis_type':analysis_type})

data = defaultdict(lambda: defaultdict(dict))
allpools = set()
# allrefs = set()
for s in summaries:
  fields = s.split('/')[prefix_len:]
  pool = fields[0]
  bc = fields[1]
  ref = fields[3].split('.')[0]
  allpools.add(pool)
  # allrefs.add(ref)
  data[ref][bc][pool] = parse_callsum(s)

allpools = sorted(allpools)

for ref in sorted(data.keys()):
  for bc in sorted(data[ref].keys()):
    call_per_pool = []
    for p in allpools:
      if p in data[ref][bc]:
        call_per_pool.append(data[ref][bc][p][ref]['call'])
      else:
        call_per_pool.append('X')
    print '%s\t%s\t%s' % (ref,bc,'\t'.join(call_per_pool))
        
        
with open('%s.summary.txt' % basename,'rU') as fh:
  burn = fh.next()
  info = fh.next().split('\t')


IGV_XML = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="%(urlpath)s/refs/%(ref)s.fasta" locus="All" version="5">
  <Resources>
    <Resource path="%(urlpath)s/%(pool)s/%(bc)s/ccs/%(ref)s.callable.bed" name="ccs_callable"/>  
    <Resource path="%(urlpath)s/%(pool)s/%(bc)s/ccs/%(ref)s.bam" name="ccs_reads"/>
    <Resource path="%(urlpath)s/%(pool)s/%(bc)s/ccs/%(ref)s.vcf" name="ccs_variants"/>    
  </Resources>
</Session>

"""


for ref,bcdict in data.iteritems():
  for bc,pdict in bcdict.iteritems():
    for p,call in pdict.iteritems():
      print IGV_XML % {'urlpath':urlpath,'pool':p,'bc':bc,'ref':ref}
      break

    
    for p in allpools:




