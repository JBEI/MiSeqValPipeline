#! /usr/bin/env python

from itertools import chain
import json

from barcode_html import *

import os
from glob import glob
from collections import defaultdict

def get_summary_paths(base,atype):
  globstr = '%s/pool?/B???/%s/*.summary.txt' % (base,atype)
  return glob(globstr)

def parse_callsum(sumfile):
  summaries = {}
  lines = [l.strip('\n').split('\t') for l in open(sumfile,'rU')]
  header = lines[0]
  for l in lines[1:]:
    tmp = dict((h,l[header.index(h)]) for h in header[1:])
    if tmp['call'].startswith('Fixable'):
      p1,p2 = tmp['call'].split(':')[1].split(',')
      tmp['p1'] = p1
      tmp['p2'] = p2
      tmp['call'] = tmp['call'].split(':')[0]
    summaries[l[0]] = tmp
  return summaries

def ff(v):
  return str(int(round(float(v))))

def fp(v):
  return str(round(float(v),1))

def call_display(d):
  call = d['call']
  if call == 'Pass': return {'call':'perfect','display':ff(d['mean_cov'])}
  elif call == 'Fixable': return {'call':'almost','display':ff(d['mean_cov']),'p1':d['p1'],'p2':d['p2']}
  elif call == 'Errors': return {'call':'errors','display':d['nvars']}
  elif call == 'Dips': return {'call':'dips','display':d['ndips']}
  elif call == 'Incomplete': return {'call':'incomplete','display':d['pct_cov']}
  elif call == 'Low coverage': return {'call':'lowcov','display':ff(d['mean_cov'])}
  else: return {'call':'?','display':ff(d['mean_cov'])}

def merge_calls(analysis, reflens, poollist, use_protocol="roi"):
  calltable = {}
  acalls = load_calls(analysis)
  for ref,rlen in reflens:
    calltable[ref] = {}
    for pname in poollist:
      if use_protocol in acalls[pname]:
        calltable[ref][pname] = call_display(acalls[pname][use_protocol][ref])
      else: # use_protocol was not found
        if "roi" in acalls[pname]:
          calltable[ref][pname] = call_display(acalls[pname]["roi"][ref])

  return calltable

def main(parser):
  args = parser.parse_args()
  aroot = os.path.join(args.analysis_dir, args.analysis_name)
  summaries = get_summary_paths(aroot,args.analysis_type)
  
  allpools = set()
  data = defaultdict(lambda: defaultdict(dict))
  for s in summaries:
    fields = s.replace(aroot,'').split('/')[1:]
    ref  = fields[3].split('.')[0]
    pool = fields[0]
    bc   = fields[1]
    allpools.add(pool)
    data[ref][bc][pool] = parse_callsum(s)

  allpools = sorted(allpools)
  
  allrefs = set()
  analysis_outcomes = defaultdict(dict)
  for ref in sorted(data.keys()):
    for bc in sorted(data[ref].keys()):
      newref = '%s|%s' % (ref,bc)
      allrefs.add(newref)
      for p in allpools:
        if p in data[ref][bc]:
          call_info = call_display(data[ref][bc][p][ref])
          call_info['igv'] = 'results/%s.%s.%s.xml' % (ref,bc,p.strip('pool'))
          analysis_outcomes[newref][p] = call_info
  
  allrefs = sorted(allrefs)
  
  htmlstr = HTML_HEAD
  htmlstr += HTML_BODY % {'aname':args.analysis_name}
  htmlstr += JS_INCLUDE % {'aname':args.analysis_name,
                           'baseurl':args.url_path,
                           'pools_json':json.dumps([{'name':p} for p in allpools]),
                           'outcomes_json':json.dumps(analysis_outcomes),
                           'references_json':json.dumps([{'name':r,'length':0} for r in allrefs]),                           
                          }
  htmlstr += JS_RESULT
  
  print >>args.outfile, htmlstr

if __name__=='__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser()
  parser.add_argument('--analysis_type',default="ccs")
  parser.add_argument('--analysis_dir',default="/global/dna/shared/synbio/PacBioAnalysis")
  parser.add_argument('--url_path',default="http://synbio.jgi-psf.org:8077/analysis/")
  parser.add_argument('analysis_name')
  parser.add_argument('outfile',nargs='?',type=argparse.FileType('w'),default=sys.stdout)  
  main(parser)
