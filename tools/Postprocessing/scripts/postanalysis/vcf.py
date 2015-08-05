import os
import sys
import subprocess

'''
Old version of make_vcf designed for covvars

def make_vcf(covvars,outfile=None,jbrowseVCF=False):
  prev_mask = os.umask(002)
  if outfile is not None:
    outh = open(outfile,'w')
  else:
    outh = sys.stdout
  typeAsID = jbrowseVCF
  print >>outh, "##fileformat=VCFv4.1"
  print >>outh, "##source=coverageVariantFinder"
  #print >>outh, "reference=file"
  print >>outh, '##INFO=<ID=COVVAR,Number=1,Type=String,Description="Type of coverage variant">'
  print >>outh, '##INFO=<ID=COVQ,Number=.,Type=Integer,Description="Coverage dip quality scores">'  
  print >>outh, '##INFO=<ID=DP,Number=1,Type=Integer,Description="Mean depth of coverage">'
  print >>outh, '##ALT=<ID=DP,Number=1,Type=Integer,Description="Mean depth of coverage">'  
  print >>outh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  for refname in covvars.keys():
    if 'variants' in covvars[refname]:
      for v in sorted(covvars[refname]['variants'],key=lambda x: x.pos):
        print >>outh, v.simpleVCF(typeAsID=typeAsID)
  
  if not outh == sys.stdout: outh.close()
  os.umask(prev_mask)

'''

def make_vcf(variants,outfile=None,jbrowseVCF=False):
  prev_mask = os.umask(002)
  if outfile is not None:
    outh = open(outfile,'w')
  else:
    outh = sys.stdout
  typeAsID = jbrowseVCF
  print >>outh, "##fileformat=VCFv4.1"
  print >>outh, "##source=coverageVariantFinder"
  #print >>outh, "reference=file"
  print >>outh, '##INFO=<ID=COVVAR,Number=1,Type=String,Description="Type of coverage variant">'
  print >>outh, '##INFO=<ID=COVQ,Number=.,Type=Integer,Description="Coverage dip quality scores">'  
  print >>outh, '##INFO=<ID=DP,Number=1,Type=Integer,Description="Mean depth of coverage">'
  print >>outh, '##ALT=<ID=DP,Number=1,Type=Integer,Description="Mean depth of coverage">'  
  print >>outh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  for v in sorted(sorted(variants,key=lambda x:x.pos),key=lambda y:y.chrom):
    print >>outh, v.simpleVCF(typeAsID=typeAsID)
  
  if not outh == sys.stdout: outh.close()
  os.umask(prev_mask)
  if outfile is None:
    return oufile
  else:
    return os.path.abspath(outfile)


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

