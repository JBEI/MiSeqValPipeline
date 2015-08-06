#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('conffile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
#parser.add_argument('--archive_path', default="/global/dna/projectdirs/RD/synbio/Analysis", help="Final destination directory for analysis files.")
parser.add_argument('--archive_path', default="/global/dna/shared/synbio/PacBioAnalysis", help="Final destination directory for analysis files.")
parser.add_argument('--smrtmodule', default="smrtanalysis/2.1.1")
parser.add_argument('--notify', default="mlbendall@lbl.gov", help="Email to be notifed when analysis is complete. Passed to SGE.")

args = parser.parse_args()

import os
from lxml import etree
from Bio import SeqIO


SMRTMODULE      = args.smrtmodule
NOTIFY          = args.notify
ARCHIVE_PATH    = args.archive_path
DEFAULT_TMPPATH = os.path.join(os.path.expandvars("$SCRATCH"),"SYNBIO/production")

USE_JAMO = False

from setup import settings, references, inputs, jobs

def setup_analysis(name,reffile,dest,copy_ref=True):
  ''' Dest is path '''
  apath = os.path.join(os.path.abspath(dest),name)
  assert not os.path.exists(apath), "Analysis directory already exists!"
  os.mkdir(apath)
  assert os.path.isfile(reffile), "Reference file does not exist!"
  if copy_ref:
    newref = os.path.join(apath,'references.fasta')
    with open(newref,'w') as outh:
      for s in SeqIO.parse(reffile,'fasta'):
        z = SeqIO.write(s,outh,'fasta')
    reffile = newref
  #--- Build reference dictionary
  refdict = references.create_reference(seqfile=reffile,destpath=os.path.join(apath,"references"),smrtmodule=SMRTMODULE)
  assert refdict is not None, "Error creating reference"
  #--- 
  with open(os.path.join(apath,'summarize.sh'),'w') as outh:  
    print >>outh, jobs.generate_summarysh(archive_path=ARCHIVE_PATH,analysis_name=name,notify=NOTIFY)
  
  return apath,refdict

def setup_pool(name,samples,apath,useJAMO=True):
  ppath = os.path.join(apath,name)
  assert not os.path.exists(ppath), "Pool directory already exists!"
  os.mkdir(ppath)
  h5files = []
  for samp in samples:
    f = inputs.find_h5s(samp,useJAMO=useJAMO)
    if f is None: sys.exit(1)
    h5files.extend(f)
    ###smrtcells.extend(inputs.find_smrtcells(samp))
  with open(os.path.join(ppath,'input.xml'),'w') as outh:
    print >>outh, inputs.generate_input_xml(h5files)
  return ppath

def setup_job(name, protocol, refdict, ppath, jparams):
  jpath = os.path.join(ppath,name)
  assert not os.path.exists(jpath), "Job directory already exists!"
  os.mkdir(jpath)
  settings_fn = 'settings.%s.xml' % protocol
  with open(os.path.join(jpath,settings_fn),'w') as outh:
    print >>outh, settings.generate_settings(refdict,protocol,**jparams)
  with open(os.path.join(jpath,'job.sh'),'w') as outh:
    print >>outh, jobs.generate_jobsh(protocol,jobname=name,settings_xml=settings_fn,smrtmodule=SMRTMODULE)
  return jpath

JOBCMD = '''cd %(job_path)s && qsub -N SBA_%(aname)s_%(pname)s_%(jname)s job.sh '''
SUMCMD = '''cd %(analysis_path)s && qsub -hold_jid SBA_%(aname)s_* -N SBA_summary_%(aname)s summarize.sh '''
# -m e -M %(email)s

if __name__=='__main__':
  omask = os.umask(002)
  fh = args.conffile
  parser = etree.XMLParser(remove_blank_text=True)
  tree = etree.parse(fh,parser)
  root = tree.getroot()
  assert root.tag == 'SBAnalysis', "ERROR: Root tag is not SBAnalysis"

  for analysis in root:
    assert analysis.tag == 'analysis', "ERROR: incorrect analysis tag"
    aname = analysis.attrib['name']
    # Analysis will be run in tmpdest and will be migrated to final location
    if 'location' in analysis.attrib: # Location is specified in XML file
      assert os.path.exists(analysis.attrib['location']), "Destination directory does not exist!"
      tmpdest = os.path.abspath(analysis.attrib['location'])
    else:
      tmpdest = os.path.abspath(DEFAULT_TMPPATH)
    
    apath,refdict = setup_analysis(aname, analysis.attrib['reference'], tmpdest)
    submit_commands = []
    for pool in analysis:
      pname = pool.attrib['name']
      samples = pool.attrib['samples'].split(',')
      ppath = setup_pool(pname,samples,apath,useJAMO=USE_JAMO)
      for job in pool:
        jname = job.attrib['name']
        jparams = dict((k,v) for k,v in job.attrib.iteritems() if k not in ['name','protocol'])
        jpath = setup_job(jname, job.attrib['protocol'], refdict, ppath, jparams)
        submit_commands.append(JOBCMD % {'job_path':jpath, 'aname':aname, 'pname':pname, 'jname':jname})
    # Output configuration file for analysis
    
    analysis.attrib['location'] = tmpdest
    analysis.attrib['reference'] = os.path.join(tmpdest,'references/seqs/sequence/seqs.fasta')
    with open(os.path.join(apath,'config.xml'),'w') as outh:
      print >>outh, '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'''
      print >>outh, '''<SBAnalysis>'''
      print >>outh, etree.tostring(analysis,pretty_print=True)
      print >>outh, '''</SBAnalysis>'''
    # Output job script
    with open(os.path.join(apath,'launch.sh'),'w') as outh:
      for cmd in submit_commands:
        print >>outh, cmd
      print >>outh, SUMCMD % {'analysis_path':apath,'aname':aname}
  
  os.umask(omask)
