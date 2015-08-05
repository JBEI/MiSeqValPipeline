#! /usr/bin/env python
import sys
import argparse
from collections import defaultdict
import re
import os
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

GPHEADER = """#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=4:00:00
#$ -P gentech-rnd.p
#$ -pe pe_slots 8
#$ -w e
#$ -j y
#$ -R y
#$ -N %(aname)s_align%(jobnum)s"
#$ -hold_jid "%(aname)s_binning*"

umask 0002
export _JAVA_OPTIONS="-Xmx12G"

JOBSH=%(job_sh)s

"""

def main(parser):
  args = parser.parse_args()

  lines = [l.strip('\n').split('\t') for l in args.barcode_table]
  header = lines[0]
  lines = lines[1:]

  ### Setup reference dictionary
  refs = {}
  custom = []
  for l in lines:
    fields = l[header.index('Reference')].strip('>').split('#')
    if len(fields)==1:
      custom.append(fields[0])
      continue
    sname,seq = fields
    if sname in refs:
      assert refs[sname] == seq, "ERROR: %s exists with different sequence" % sname
    else:
      refs[sname] = seq

  refpaths = {}  
  refdir = "refs"
  if not os.path.exists(refdir): os.mkdir(refdir)
  for sname,seq in refs.iteritems():
    fn = os.path.join(refdir,'%s.fasta' % sname.strip('>'))
    refpaths[sname] = os.path.abspath(fn)    
    if os.path.exists(fn): continue
    with open(fn,'w') as outh:
      print >>outh, '>%s\n%s' % (sname,'\n'.join(seq[i:i+60] for i in range(0,len(seq),60)))
    
    if os.path.exists('%s.fai' % fn): os.remove('%s.fai' % fn)
    output = subprocess.call(['samtools','faidx',fn])
    dictfile = '%s.dict' % fn.strip('.fasta')
    if os.path.exists(dictfile): os.remove(dictfile)  
    output = subprocess.call(['picard','CreateSequenceDictionary','R=%s' % fn,'O=%s' % dictfile,])
  
  # Dictionary mapping each pool -> { barcode number -> [references] }
  pool_barcodes = defaultdict(lambda: defaultdict(list))
  for l in lines:
    for p in decodeList(l[header.index('Pools')]):
      pool_barcodes[p][int(l[header.index('Barcode')])].append(l[header.index('Reference')].strip('>').split('#')[0])

  # Get the PacBio sample ID for each pool
  # Assume that pool1 has sample ID provided by command
  # Assume that other pools are sequentially numbered
  last_pool = max(int(_) for _ in pool_barcodes.keys())
  pool_sample = get_pool_map(args.sample_id,last_pool)
  
  # Create variable exports
  varstrings = []
  for pn in pool_barcodes.keys():
    for bn in pool_barcodes[pn].keys():
      for ref in pool_barcodes[pn][bn]:
        if ref in refs:
          params = {'samp':  pool_sample[pn],
                    'pool': 'pool%d' % pn,
                    'bcnum': str(bn),
                    'ref':   refpaths[ref],
                    'aname': args.analysis_name,
                   }
          varstrings.append('export %s' % ' '.join('%s=%s' % t for t in params.iteritems()) )
  
  jobnum = 1
  for i in range(0,len(varstrings),args.alignments_per_job):
    with open('job%d.sh' % jobnum,'w') as outh:
      print >>outh, GPHEADER % {'job_sh':args.job_sh,'aname':args.analysis_name,'jobnum':jobnum}
      for vs in varstrings[i:(i+args.alignments_per_job)]:
        print >>outh, vs
        print >>outh, ". $JOBSH"
    cmd = ['qsub','job%d.sh' % jobnum]
    if args.no_qsub:
      print 'NOT SUBMITTED: %s' %  ' '.join(cmd)
    else:
      subprocess.call(cmd)
    jobnum += 1

if __name__ == '__main__':
  parser = argparse.ArgumentParser(prog='setup_barcode_alignments', description='Set up alignments for barcode')
  parser.add_argument('--job_sh', default="/global/dna/projectdirs/RD/synbio/Repos/JGISB_analysis.git/jobs/barcode_ccs_align.sh")  
  parser.add_argument('--no_qsub',action='store_true')
  parser.add_argument('--alignments_per_job', type=int, default=50, help='Number of alignments per job')
  parser.add_argument('sample_id', help='Sample ID for pool 1')
  parser.add_argument('analysis_name', help='Name for analysis')  
  parser.add_argument('barcode_table', nargs='?', type=argparse.FileType('rU'), help='table with barcode information', default=sys.stdin)  
  main(parser)
