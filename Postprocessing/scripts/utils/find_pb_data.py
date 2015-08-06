#! /usr/bin/env python

import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--jamo', action='store_true')
parser.add_argument('--filetype', default='h5')
parser.add_argument('sample_name')
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
args = parser.parse_args()

if not args.sample_name:
  parser.print_help()
  sys.exit(1)

'''
Constants
'''
JAMO_URL = "https://sdm2.jgi-psf.org"
SMRT_URL = "http://pacbio.jgi-psf.org/smrtportal/api"
SMRT_HEADERS = { "Accept":"application/json","Content-Type":"applications/json"}

# Contains filename patterns for various filetypes
FILE_PATTERNS = {'h5':"Analysis_Results/*.bax.h5",
                 'metadata':"*.metadata.xml"}

def sample_paths_jamo(samp):
  '''
    Find paths to h5 files using JAMO
  '''
  # Load JAMO module and import sdm_curl
  import EnvironmentModules as EnvMod
  EnvMod.module(['load','jamo'])
  # EnvMod.module(['load','hdf5'])
  for x in os.getenv("PYTHONPATH").split(":"):
    if x not in sys.path:
      sys.path.append(x)

  from sdm_curl import Curl  
  curl = Curl(JAMO_URL)
  files = curl.post('api/metadata/query',data={'file_type':'h5',
                                               'metadata.sample_name':{'$regex':samp}})
  if len(files)==0:
    print >>sys.stderr, '[ERROR] No samples matching "%s"' % samp
    sys.exit(1)
  
  h5_paths = (os.path.join(f['file_path'],f['file_name']) for f in files if f['file_name'].endswith('bax.h5'))
  return sorted(h5_paths)

def sample_paths_smrtportal(samp,pattern):
  '''
    Find paths to h5 files using SMRTportal
  '''
  import urllib2
  import json
  from glob import glob
  
  # API query options
  options = {
    "filters": {
      "rules": [{"field":"sampleName","op":"bw","data":samp}]
    },
    "columnNames":["primaryTimestamp","primaryConfigFileName","sampleName","collectionPathUri","runName","inputId"],
    "rows":"0",
  }

  #--- Access API
  args = 'options=%s' % json.dumps(options)
  request  = urllib2.Request('%s/inputs' % SMRT_URL, args, SMRT_HEADERS)
  response = urllib2.urlopen(request)
  the_page = response.read()
  table = json.loads(the_page)
  if table['total']==0:
    print >>sys.stderr, '[ERROR] No samples matching "%s"' % samp
    sys.exit(1)

  h5_paths = []
  run_paths = sorted(set(r['collectionPathUri'] for r in table['rows']))
  for path in run_paths:
    gf = glob(os.path.join(path,pattern))
    if len(gf)==0:
      print >>sys.stderr, '[WARNING] No H5 files found in "%s"' % path
    else:
      h5_paths.extend(gf)
  
  return h5_paths

if args.jamo:
  if args.filetype == 'h5':
    print >>args.outfile, '\n'.join(sample_paths_jamo(args.sample_name))
  else:
    print >>sys.stderr, '[ERROR] File type "%s" not available for JAMO' % args.filetype
    sys.exit(1)
else:
  if args.filetype in FILE_PATTERNS:
    pattern = FILE_PATTERNS[args.filetype]
    print >>args.outfile, '\n'.join(sample_paths_smrtportal(args.sample_name,pattern))
  else:
    print >>sys.stderr, '[ERROR] Unknown file type "%s"' % args.filetype
    sys.exit(1)

"""
from glob import glob
import re
import os

h5_files = [l.strip() for l in open('input.fofn')]
rundirs = sorted(set('/'.join(f.split('/')[:-2]) for f in h5_files))
for rundir in rundirs:
metadata_xml = glob(os.path.join(rundir,"*.metadata.xml"))[0]


movies = sorted(set('.'.join(f.split('.')[:-3]) for f in h5_files))
for movie in movies:
rundir = '/'.join(movie.split('/')[:-2])
"""
