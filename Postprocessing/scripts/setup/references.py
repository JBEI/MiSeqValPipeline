
import os
import sys
from subprocess import check_call, CalledProcessError
REFCMD = '''module load %(smrtmodule)s && export _JAVA_OPTIONS="-Xmx2G" && referenceUploader -c -p %(destpath)s -n %(refname)s -f %(seqfile)s -k 'createSequenceDictionary' -x 'samtools faidx' -s 'sawriter' '''

std_ref_opts = {'smrtmodule':'smrtanalysis/2.1.1',
                'destpath':'./references',
                'refname':'seqs',
                'seqfile':'./references.fasta',
               }


def create_reference(**kwargs):
  omask = os.umask(002)
  opts = dict(std_ref_opts)
  for k,v in kwargs.iteritems():
    opts[k] = v
  
  # Check file and directory
  if not os.path.exists(opts['seqfile']): return None  # no seqfile!
  if not os.path.isdir(opts['destpath']):
    if os.path.isfile(opts['destpath']):  return None  # dest is a file!
    os.makedirs(opts['destpath'])  # create directory if needed
  
  # use absolute paths
  opts['destpath'] = os.path.abspath(opts['destpath'])
  opts['seqfile'] = os.path.abspath(opts['seqfile'])  
  
  print >>sys.stderr, REFCMD % opts
  try:
    check_call(REFCMD % opts, shell=True)
  except CalledProcessError:
    print >>sys.stderr, "Error creating reference. Call was: %s" % (REFCMD % opts)
    return None
  finally:
    os.umask(omask)
  return os.path.join(opts['destpath'],opts['refname'])    


"""
For JBrowse refs
# create results directory
[ ! -d results ] && mkdir results

# initialize jbrowse with references
[ ! -d jbrowse ] && mkdir jbrowse
module purge
module use /projectb/projectdirs/RD/synbio/.modules
module load jbrowse

prepare-refseqs.pl --fasta /projectb/projectdirs/RD/synbio/pacbioSB/testing/A000176/references.fasta --out /projectb/projectdirs/RD/synbio/pacbioSB/testing/A000176/jbrowse || exit $?;
"""