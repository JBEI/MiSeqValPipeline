from cStringIO import StringIO
from subprocess import Popen, PIPE
import re

pathToPipeline = 'path/to/pipeline.py'

def inputGD(refname,variants,seq):
  ''' Creates StringIO (string buffer) to use as input for JGISB_Design_Fixing_Primers.pl
      refname - reference name
      variants - list of variant objects
      seq - biopython SeqRecord
  '''
  vcfin = StringIO()
  for v in sorted(variants,key=lambda x:x.pos):
    print >>vcfin, v.simpleVCF()

  print >>vcfin, '#'
  for i in range(0,len(seq),60):
    print >>vcfin, str(seq[i:i+60].seq)

  return vcfin

outRE = re.compile('^(?P<score>\-?\d+)(\s+(?P<primer1>[ATCGatcg]+)\s+(?P<primer2>[ATCGatcg]+)\s*)?$')

def runGD(vcfin):
  p = Popen(['perl', pathToPipeline+'/tools/DNAssemble/bin/JGISB_Design_Fixing_Primers.pl'],stdout=PIPE,stdin=PIPE,stderr=PIPE)
  stdout,stderr = p.communicate(input=vcfin.getvalue())
  vcfin.close()
  m = outRE.match(stdout)
  if m:
    return m.groupdict()
  else:
    return {'score':'-1'}
    
    
