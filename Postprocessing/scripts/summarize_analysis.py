#!/usr/bin/env python

import os, sys
import argparse
from urlparse import urljoin
parser = argparse.ArgumentParser()
parser.add_argument('conffile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('--fspath')

args = parser.parse_args()

FSPATH  = '%s/' % args.fspath.rstrip('/')

from lxml import etree
from Bio import SeqIO
from summarize.indexhtml import make_index
from summarize.merge import merge_calls, best_calls
from summarize.excel import create_result_workbook

''' Constants '''
IGVXML = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="%(genome_url)s" locus="All" version="5">
  <Resources>
%(resources)s
  </Resources>
</Session>'''


jobfiles = {
  'ccs':[('callable.bed','callable'),('aligned_reads.bam','reads'),('snps.gatk.vcf','variants'),],
  'bwa_dir':[('callable.bed','callable'),('aligned_reads.bam','reads'),('snps.gatk.vcf','variants'),],  
  'sub':[('coverage.bed','coverage'),('aligned_reads.bam','reads'),('variants.gff.gz','variants'),],
}

''' Functions '''
def igv_xml(genome,resources):
  resstr = ''
  for path,name in resources:
    resstr += '    <Resource path="%s" name="%s"/>\n' % (path,name)
  return IGVXML % {'resources':resstr,'genome_url':genome}

def mkdir_p(path):
  if not os.path.exists(path): os.mkdir(path)

if __name__=='__main__':
  omask = os.umask(002)
  fh = args.conffile
  parser = etree.XMLParser(remove_blank_text=True)
  tree = etree.parse(fh,parser)
  root = tree.getroot()
  assert root.tag == 'SBAnalysis', "ERROR: Root tag is not SBAnalysis"
  assert len(root.getchildren())==1, "ERROR: Only one analysis can be migrated"
  analysis = root[0]
  assert analysis.tag == 'analysis', "ERROR: incorrect analysis tag"

  adir = os.path.join(analysis.attrib['location'],analysis.attrib['name'])
  assert os.path.exists(adir), "Analysis directory does not exist!"

  reffile = analysis.attrib['reference']
  assert os.path.exists(reffile), "Reference file does not exist!"
  refurl = reffile

  try:
    #--- Setup analysis results directory  
    resultdir = os.path.join(adir,'results')
    mkdir_p(resultdir)

    reflens = [(s.id,len(s)) for s in SeqIO.parse(analysis.attrib['reference'],'fasta')]
    refnames = [r[0] for r in reflens]
    poollist = [pool.attrib['name'] for pool in analysis]
  
    #--- Make the IGV files
    igvurls = {}
    for pool in analysis:
      pool_resources = []
      pdir = os.path.join(adir,pool.attrib['name'])
      for job in pool:
        jdir = os.path.join(pdir,job.attrib['name'])
        for f,rname in jobfiles[job.attrib['protocol']]:
          resource_file = os.path.join(jdir,f)
          if os.path.exists(resource_file):
            pool_resources.append((resource_file,'%s_%s' % (job.attrib['protocol'],rname)))
      igvfile = os.path.join(resultdir,'%s.igv.xml' % pool.attrib['name'])
      with open(igvfile,'w') as outh:
        print >>outh, igv_xml(refurl,pool_resources)
      igvurls[pool.attrib['name']] = igvfile


    #--- Summarize calls
    calltable = merge_calls(analysis,reflens,poollist)
    bestbets = best_calls(calltable)

    #--- Create the excel file
    excelfile = os.path.join(resultdir,'%s.xlsx' % analysis.attrib['name'])
    create_result_workbook(refnames, poollist, calltable, bestbets, excelfile)
    
    #--- Make the index.html file
    htmlfile = os.path.join(adir,'index.html')
    with open(htmlfile,'w') as outh:
      print >>outh, make_index(analysis,'localhost',reflens,calltable,bestbets)

  except (IOError,OSError) as e:
    if e.args[1] == 'Read-only file system':
      sys.exit('Destination "%s" is mounted as read-only.' % args.destination)
    else: raise
  finally:
    os.umask(omask)

