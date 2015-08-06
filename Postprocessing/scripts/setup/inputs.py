import os
import sys
import glob
import urllib2
import json


import EnvironmentModules as EnvMod
EnvMod.module(['load','jamo'])
EnvMod.module(['load','hdf5'])
for x in os.getenv("PYTHONPATH").split(":"):
  if x not in sys.path:
    sys.path.append(x)

"""
Following is the resulting XML file for SMRT Pipe v2.1:

<?xml version="1.0"?>
<pacbioAnalysisInputs>
  <dataReferences>
    <url ref="run:0000000-0000"><location>/mnt/data/2770276/0006/Analysis_Results/m130512_050747_42209_c100524
962550000001823079609281357_s1_p0.2.bax.h5</location></url>
    <url ref="run:0000000-0001"><location>/mnt/data/2770276/0006/Analysis_Results/m130512_050747_42209_c100524
962550000001823079609281357_s1_p0.3.bax.h5</location></url>
    <url ref="run:0000000-0002"><location>/mnt/data/2770276/0006/Analysis_Results/m130512_050747_42209_c100524
962550000001823079609281357_s1_p0.1.bax.h5</location></url>
  </dataReferences>
</pacbioAnalysisInputs>

For SMRT Pipe versions before v2.1:

<?xml version="1.0"?>
<pacbioAnalysisInputs>
 <dataReferences>
    <url ref="run:0000000-0000"><location>/share/data/
    /share/data/run_1 m100923_005722_00122_c15301919401091173_s0_p0.bas.h5
    <url ref="run:0000000-0001"><location>/share/data/
    /share/data/run_2/m100820_063008_00118_c04442556811011070_s0_p0.bas.h5
 </dataReferences>
</pacbioAnalysisInputs>

"""

''' Constants '''
XML_START = '''<?xml version="1.0"?>
<pacbioAnalysisInputs>
  <dataReferences>
'''
XML_ENTRY = '''    <url ref="run:0000000-%(run)04d"><location>%(fn)s</location></url>
'''
XML_END  = '''  </dataReferences>
</pacbioAnalysisInputs>
'''

''' Utility functions '''
def find_smrtcells_jamo_dna(sample_name):
  ''' Returns absolute path for all H5 files matching sample_name '''
  from sdm_curl import Curl
  filenames = []
  curl = Curl('https://sdm2.jgi-psf.org')
  files = curl.post('api/metadata/query',data={'metadata.sample_name':sample_name,'file_type':'h5'})
  if not files:
    return None
  return [os.path.join(f['file_path'],f['file_name']) for f in files]

def find_smrtcells_jamo(sample_name):
  ''' Returns absolute path for all H5 files matching sample_name '''
  from sdm_curl import Curl
  filenames = []
  curl = Curl('https://sdm2.jgi-psf.org')
  files = curl.post('api/metadata/query',data={'metadata.sample_name':sample_name,'file_type':'h5'})
  if not files:
    return None
  locations = list(set(f['metadata']['sdm_smrt_cell']['fs_location'] for f in files))
  return locations

def find_smrtcells(sample_name, pattern = "Analysis_Results/*.ba?.h5"):
  h5files = []
  runs = find_smrtcells_jamo(sample_name)
  assert runs is not None and len(runs)>0, "ERROR: No runs found for %s" % sample_name
  for path in runs:
    gf = glob.glob(os.path.join(path,pattern))
    assert len(gf)>0, "ERROR: No files matching pattern: %s" % os.path.join(path,pattern)
    h5files.extend(gf)
  return h5files

def generate_xml_jamo(filenames):
  ''' Returns XML string to make a smrtpipe input.xml '''
  from subprocess import check_output
  runidx = 0
  xmlstr = ''
  for fn in filenames:
    if 'PulseData' in [l.split()[0] for l in check_output(['h5ls',fn]).split('\n') if l]:
      xmlstr += XML_ENTRY % {'run':runidx,'fn':fn}
      runidx += 1
  return XML_START + xmlstr + XML_END


###
SMRT_URL = "http://pacbio.jgi-psf.org/smrtportal/api"
SMRT_HEADERS = { "Accept":"application/json","Content-Type":"applications/json"}

def find_smrtcells_smrtportal(sample_prefix):
  options = {
    "filters": {
      "rules": [{"field":"sampleName","op":"bw","data":sample_prefix}]
    },
    "columnNames":["primaryTimestamp","primaryConfigFileName","sampleName","collectionPathUri","runName","inputId"],
    "rows":"0",
  }
  
  url = '%s/inputs' % SMRT_URL
  args = 'options=%s' % json.dumps(options)
  headers = SMRT_HEADERS
  
  #--- Access API
  request  = urllib2.Request(url,args,headers)
  response = urllib2.urlopen(request)
  the_page = response.read()
  table = json.loads(the_page)
  if table['total']:
    return list(set(r['collectionPathUri'] for r in table['rows']))
  else:
    return None

def find_h5s(sample_prefix, pattern = "Analysis_Results/*.ba?.h5",useJAMO=True):
  h5files = []
  if useJAMO:
    runs = find_smrtcells_jamo(sample_prefix)
  else:
    runs = find_smrtcells_smrtportal(sample_prefix)
  #assert runs is not None and len(runs)>0, "ERROR: No runs found for %s" % sample_prefix
  if runs is None or len(runs)==0:
    print >>sys.stderr, '[ERROR] No SMRTcells matching "%s"' % sample_prefix
    return None
  for path in runs:
    gf = glob.glob(os.path.join(path,pattern))
    if len(gf)==0:
      print >>sys.stderr, '[ERROR] H5 files not found matching "%s"' % os.path.join(path,pattern)
      return None
    #assert len(gf)>0, "ERROR: No files matching pattern: %s" % os.path.join(path,pattern)
    h5files.extend(gf)
  return h5files

def generate_input_xml(filenames):
  ''' Returns XML string to make a smrtpipe input.xml '''
  from subprocess import check_output
  runidx = 0
  xmlstr = ''
  for fn in filenames:
    if 'PulseData' in [l.split()[0] for l in check_output(['h5ls',fn]).split('\n') if l]:
      xmlstr += XML_ENTRY % {'run':runidx,'fn':fn}
      runidx += 1
  return XML_START + xmlstr + XML_END

