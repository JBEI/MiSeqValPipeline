import json
import fcntl
import shutil
import os

def create_new_track(key,url,type):
  ''' Creates a python dictionary that can be serialized as a jbrowse track description '''
  newtrack = {'label':key, 'key':key}
  if type in ['coverage','alignment']:
    newtrack['storeClass']  = 'JBrowse/Store/SeqFeature/BAM'
    newtrack['urlTemplate'] = url
    newtrack['baiUrlTemplate'] = '%s.bai' % url
    if type == 'coverage':
      newtrack['type'] = 'JBrowse/View/Track/SNPCoverage'
    else:
      newtrack['type'] = 'JBrowse/View/Track/Alignments2'
      newtrack['maxHeight'] = 'JBrowse/View/Track/Alignments2'
      newtrack['style'] = {
            'color_fwd_strand': '#009900',
            'color_rev_strand': '#999900',
          },
  elif type=='vcf':
    newtrack['storeClass'] = 'JBrowse/Store/SeqFeature/VCFTabix'
    newtrack['type'] = 'JBrowse/View/Track/HTMLVariants'
    newtrack['urlTemplate'] = url
    newtrack['tbiUrlTemplate'] = '%s.tbi' % url
  else:
    print "Unknown type: %s" % type
  return newtrack

def insert_track(jobj,newtrack):
  ''' Inserts python track description into python trackList object
      Replaces existing tracks based on "key" field
  '''
  matching = [track for track in jobj['tracks'] if track['key'] == newtrack['key']]
  assert len(matching) <= 1, "too many matching tracks!"
  if len(matching)==1:
    jobj['tracks'].remove(matching[0])
  jobj['tracks'].append(newtrack)
  jobj['tracks'].sort(key=lambda x:x['key'])

def update_tracklist(tracklist,newtracks):
  ''' Inserts list of python tracks into jbrowse tracklist
      1) Opens and locks tracklist file
      2) Loads json as python object
      3) Inserts new tracks
      4) Serializes python object and writes to jbrowse tracklist
  '''
  with open(tracklist,'r+') as fh:
    fcntl.flock(fh,fcntl.LOCK_EX)
    jobj = json.load(fh)
    for nt in newtracks:
      insert_track(jobj,nt)
    fh.seek(0)
    print >>fh, json.dumps(jobj)

def backup_file(filename,overwrite=False):
  ''' Creates backup of file '''
  if not os.path.exists(filename): return False
  if overwrite:
    shutil.copy(filename,'%s.bak' % filename)
    return '%s.bak' % filename
  else:
    index = 1
    while os.path.exists('%s.%02d.bak' % (filename,index)):
      index += 1
    shutil.copy(filename,'%s.%02d.bak' % (filename,index))
    return '%s.%02d.bak' % (filename,index)

def delete_jbrowse(jbrowse_dir,overwrite=True):
  if os.path.exists(jbrowse_dir):
    if not overwrite: return False
    shutil.rmtree(jbrowse_dir)
  return True

'''
import json

def create_new_track(key,url,type):
  newtrack = {'label':key, 'key':key}
  if type in ['coverage','alignment']:
    newtrack['storeClass']  = 'JBrowse/Store/SeqFeature/BAM'
    newtrack['urlTemplate'] = url
    newtrack['baiUrlTemplate'] = '%s.bai' % url
    if type == 'coverage':
      newtrack['type'] = 'JBrowse/View/Track/SNPCoverage'
    else:
      newtrack['type'] = 'JBrowse/View/Track/Alignments2'
      newtrack['maxHeight'] = 'JBrowse/View/Track/Alignments2'
      newtrack['style'] = {
            'color_fwd_strand': '#009900',
            'color_rev_strand': '#999900',
          },
  elif type=='vcf':
    newtrack['storeClass'] = 'JBrowse/Store/SeqFeature/VCFTabix'
    newtrack['type'] = 'JBrowse/View/Track/HTMLVariants'
    newtrack['urlTemplate'] = url
    newtrack['tbiUrlTemplate'] = '%s.tbi' % url
  else:
    print "Unknown type: %s" % type
  return newtrack

def insert_track(jobj,newtrack):
  matching = [track for track in jobj['tracks'] if track['key'] == newtrack['key']]
  assert len(matching) <= 1, "too many matching tracks!"
  if len(matching)==1:
    jobj['tracks'].remove(matching[0])
  jobj['tracks'].append(newtrack)
  jobj['tracks'].sort(key=lambda x:x['key'])


import subprocess
import os

def index_vcf(vcf_file,gz_file=None):
  prev_mask = os.umask(002)
  SAMTOOLS_BIN = '/projectb/projectdirs/RD/synbio/software/samtools/DEFAULT/bin'
  if gz_file is None: gz_file = '%s.gz' % vcf_file
  data = {'sambin':SAMTOOLS_BIN,'vcf_file':vcf_file,'gz_file':gz_file}
  zipcmd = '%(sambin)s/bgzip -c %(vcf_file)s > %(gz_file)s' % data
  subprocess.call(zipcmd,shell=True)
  idxcmd = '%(sambin)s/tabix -p vcf %(gz_file)s' % data
  subprocess.call(idxcmd,shell=True)
  os.umask(prev_mask)
  return gz_file

def delete_vcf_index(vcf_file,gz_file=None):
  if gz_file is None: gz_file = '%s.gz' % vcf_file
  if os.path.exists(gz_file): os.remove(gz_file)
  # also delete index
  if os.path.exists('%s.tbi' % gz_file): os.remove('%s.tbi' % gz_file)







'''