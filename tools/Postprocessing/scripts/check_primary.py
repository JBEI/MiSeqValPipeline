#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('conffile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
#parser.add_argument('--archive_path', default="/global/dna/projectdirs/RD/synbio/Analysis", help="Final destination directory for analysis files.")
parser.add_argument('--archive_path', default="/global/dna/shared/synbio/PacBioAnalysis", help="Final destination directory for analysis files.")
args = parser.parse_args()

import os
from lxml import etree
from Bio import SeqIO

ARCHIVE_PATH    = args.archive_path
DEFAULT_TMPPATH = os.path.join(os.path.expandvars("$SCRATCH"),"SYNBIO/production")

from setup import settings, references, inputs, jobs

def check_analysis(name,reffile,dest,copy_ref=True):
  _isready = True
  apath = os.path.join(os.path.abspath(dest),name)
  if os.path.exists(apath):
    _isready = False
    print >>sys.stderr, '[WARNING] Analysis directory "%s" already exists!' % apath
  if not os.path.isfile(reffile):
    _isready = False
    print >>sys.stderr, '[WARNING] Reference file "%s" does not exist!' % reffile
  return apath,_isready

def check_pool(name,samples,apath,useJAMO=True):
  _isready = True
  h5files = []
  for samp in samples:
    f = inputs.find_h5s(samp,useJAMO=useJAMO)
    if f is None:
      _isready = False
  return _isready

if __name__=='__main__':
  fh = args.conffile
  parser = etree.XMLParser(remove_blank_text=True)
  tree = etree.parse(fh,parser)
  root = tree.getroot()
  assert root.tag == 'SBAnalysis', "ERROR: Root tag is not SBAnalysis"

  for analysis in root:
    isready = True  
    assert analysis.tag == 'analysis', "ERROR: incorrect analysis tag"
    aname = analysis.attrib['name']
    # Analysis will be run in tmpdest and will be migrated to final location
    if 'location' in analysis.attrib: # Location is specified in XML file
      if not os.path.exists(analysis.attrib['location']):
        print >>sys.stderr, '[WARNING] Destination directory "%s" does not exist!' % analysis.attrib['location']
        isready = False
      tmpdest = os.path.abspath(analysis.attrib['location'])
    else:
      tmpdest = os.path.abspath(DEFAULT_TMPPATH)
    
    apath,ir = check_analysis(aname, analysis.attrib['reference'], tmpdest)
    isready &= ir
    for pool in analysis:
      pname = pool.attrib['name']
      samples = pool.attrib['samples'].split(',')
      isready &= check_pool(pname,samples,apath,useJAMO=False)

    if isready:
      print >>sys.stderr, "Analysis %s is READY" % aname
      sys.exit(0)
    else:
      print >>sys.stderr, "Analysis %s is NOT READY" % aname
      sys.exit(1)
