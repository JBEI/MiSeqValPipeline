#! /usr/bin/env python

import sys
import os
import re
import json

try:
  import cPickle as pickle
except:
  import pickle

from Bio import SeqIO

''' postanalysis imports '''
import variant,vcf,jbrowse,igv
from calls import CloneCall

def standard_calls(ccs_summary,xlr_summary,ccs_cov=10,xlr_cov=30):
  ''' Standard procedure for calling clones
      Mean coverage and percent covered data from both CCS and XLR analysis are used
      GATK variant calls are used for the CCS analysis, and GenCons variant calls are used
      for the XLR analysis
      Coverage dips are used for the CCS analysis
  '''
  cc = CloneCall(mean_cov=ccs_summary['mean_cov'],pct_cov=ccs_summary['pct_cov'])

  # check whether there is coverage
  if ccs_summary['mean_cov'] < ccs_cov and xlr_summary['mean_cov'] < xlr_cov:
    return cc.set_call('lowcov',['CCS','XLR'])
  elif ccs_summary['mean_cov'] < ccs_cov:
    return cc.set_call('lowcov',['CCS'])
  elif xlr_summary['mean_cov'] < xlr_cov:
    cc.mean_cov = xlr_summary['mean_cov']
    return cc.set_call('lowcov',['XLR'])

  # check whether sequence is complete
  if ccs_summary['pct_cov'] < 1 and xlr_summary['pct_cov'] < 1:
    return cc.set_call('incomplete',['CCS','XLR'])
  elif ccs_summary['pct_cov'] < 1:
    return cc.set_call('incomplete',['CCS'])
  elif xlr_summary['pct_cov'] < 1:
    cc.pct_cov = xlr_summary['pct_cov']
    return cc.set_call('incomplete',['XLR'])

  # check whether variants were called
  ccs_gatk = [v for v in ccs_summary['variants'] if v.caller=='gatk']
  xlr_gencons = [v for v in xlr_summary['variants'] if v.caller=='gencons']
  if ccs_gatk and xlr_gencons:
    cc.variants = ccs_gatk + xlr_gencons
    return cc.set_call('errors',['gatk','gencons'])
  elif ccs_gatk:
    cc.variants = ccs_gatk
    return cc.set_call('errors',['gatk'])
  elif xlr_gencons:
    cc.variants = xlr_gencons
    return cc.set_call('errors',['gencons'])

  dips = [v for v in ccs_summary['variants'] if v.type=='cov_dip']
  if dips:
    cc.variants = dips
    return cc.set_call('dips',['covvars'])
  return cc.set_call('flawless')

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Run postanalysis for pool.')
  parser.add_argument('--pooldir',default='.')
  parser.add_argument('--reffile',default='../references.fasta')
  parser.add_argument('--ccs_cov',default=10,type=int)
  parser.add_argument('--xlr_cov',default=30,type=int)
  parser.add_argument('--base_url',default="http://synbio.jgi-psf.org/pbdata")
  parser.add_argument('--jobindex',default='jobindex.txt')
  parser.add_argument('--vcfout',default='variants.vcf')
  parser.add_argument('--analysis_root',default='pacbioSB')
  parser.add_argument('--tracklist',default='../jbrowse/trackList.json')  
  args = parser.parse_args()

  if not os.path.isdir(args.pooldir):   sys.exit('Error: directory "%s" does not exist' % args.pooldir)
  if not os.path.exists(args.reffile):  sys.exit('Error: reference file "%s" does not exist' % args.reffile)  
  if not os.path.exists(args.jobindex): sys.exit('Error: job index file "%s" does not exist' % args.jobindex)
  
  pool_path              = os.path.abspath(args.pooldir)
  analysis_path,poolname = os.path.split(pool_path)
  reference_file         = os.path.abspath(args.reffile)
  jobindex_file          = os.path.abspath(args.jobindex)
  
  # load sequences
  seqs = dict( [(s.id,s) for s in SeqIO.parse(reference_file,'fasta')] )
  
  # load job index
  jobIDs = dict([l.strip().split('\t') for l in open(jobindex_file,'rU')])


  vcf_variants = []
  calls = {}
  ''' Select the calling method: depends on the protocols used in the pool '''
  if 'ccs' in jobIDs and 'xlr' in jobIDs:   #--- CCS and XLR -> standard_calls
    print >>sys.stderr, "[ Loading job summaries ]"
    all_ready = True    
    ccs_pickle = '%s/%06d/job_summary.pickle' % (pool_path,int(jobIDs['ccs']))
    xlr_pickle = '%s/%06d/job_summary.pickle' % (pool_path,int(jobIDs['xlr']))
    if not os.path.exists(ccs_pickle):
      all_ready = False
      print >>sys.stderr, 'CCS job %06d is not ready: file "%s" does not exist' % (int(jobIDs['ccs']),ccs_pickle)
    if not os.path.exists(xlr_pickle):
      all_ready = False
      print >>sys.stderr, 'XLR job %06d is not ready: file "%s" does not exist' % (int(jobIDs['xlr']),xlr_pickle)

    if not all_ready: sys.exit(1)

    with open(ccs_pickle,'rb') as fh:
      ccs_summaries = pickle.load(fh)

    with open(xlr_pickle,'rb') as fh:
      xlr_summaries = pickle.load(fh)

    for ref in seqs.keys():
      call = standard_calls(ccs_summaries[ref],xlr_summaries[ref])
      calls[ref] = call
      if call.variants is not None:
        vcf_variants.extend(call.variants)

  #--- Here add other calling methods      
  else:
    sys.exit("Unknown protocols for making calls")

  #print >>sys.stderr, "[ Printing calls ]"
  #for ref in seqs.keys():
  #  print '%s\t%s\t%s' % (ref,calls[ref].call,calls[ref].get_reason())
  
  print >>sys.stderr, "[ Running genedesign ]"
  
  import genedesign
  for ref in seqs.keys():
    if calls[ref].call == 'errors':
      vars = calls[ref].variants
      vcfin = genedesign.inputGD(ref,vars,seqs[ref])
      #with open('%s.gdin.txt' % ref,'w') as gdout:
      #  print >>gdout, vcfin.getvalue()
      result = genedesign.runGD(vcfin)
      if result['score'] == '-1':
        pass # not fixable
      else:
        calls[ref].call = 'almost' # use set_call?
        calls[ref].fixscore = int(result['score'])
        calls[ref].primer1 = result['primer1']
        calls[ref].primer2 = result['primer2']        
    print ref
    print calls[ref]
  
  ''' Create the VCF file for the pool '''
  print >>sys.stderr, "[ Creating VCF file ]"
  outvcf = vcf.make_vcf(vcf_variants,'%s/%s' % (pool_path,args.vcfout),jbrowseVCF=True)
  if outvcf is not None: vcf.index_vcf(outvcf)
  
  
  ''' Calculate URLs for files '''
  ap = analysis_path.split('/')
  analysis_suffix = '/'.join(ap[ap.index(args.analysis_root)+1:])
  analysis_url    = '%s/%s' % (args.base_url,analysis_suffix)
  pool_url        = '%s/%s' % (analysis_url,poolname)
  variants_url    = '%s/variants.vcf' % pool_url
  ccs_bam_url     = '%s/%06d/aligned_reads.bam' % (pool_url,int(jobIDs['ccs']))
  xlr_bam_url     = '%s/%06d/aligned_reads.bam' % (pool_url,int(jobIDs['xlr']))
  callable_url    = '%s/%06d/GATK/callable.bed' % (pool_url,int(jobIDs['ccs']))

  ''' Create IGV files '''
  print >>sys.stderr, "[ Creating IGV files ]"
  igv_resources = [{'name':'%s_variants' % poolname,'path':variants_url},
                   {'name':'%s_ccs_reads' % poolname,'path':ccs_bam_url},
                   {'name':'%s_xlr_reads' % poolname,'path':xlr_bam_url},
                   {'name':'%s_callable' % poolname,'path':callable_url},
                  ]

  igvxml = igv.create_igv_session(analysis_url,igv_resources)
  if not os.path.exists('%s/results' % analysis_path): os.mkdir('%s/results' % analysis_path)
  with open('%s/results/%s.igv.xml' % (analysis_path,poolname),'w') as igvout:
    print >>igvout, igvxml.getvalue()

  ''' Add tracks to jbrowse '''
  print >>sys.stderr, "[ Adding tracks to jbrowse ]"
  tracklist_file = os.path.abspath(args.tracklist)
  if not os.path.exists(args.tracklist): sys.exit('Error: trackList file "%s" does not exist' % args.tracklist)    

  jbrowse_resources = [{'key':'%s_variants' % poolname,'url':'%s.gz' % variants_url,'type':'vcf'},
                       {'key':'%s_ccs_reads' % poolname,'url':'%s.gz' % ccs_bam_url,'type':'alignment'},
                       {'key':'%s_ccs_cov' % poolname,'url':'%s.gz' % ccs_bam_url,'type':'coverage'},                       
                       {'key':'%s_xlr_reads' % poolname,'url':'%s.gz' % xlr_bam_url,'type':'alignment'},
                       {'key':'%s_xlr_cov' % poolname,'url':'%s.gz' % xlr_bam_url,'type':'coverage'},                       
                      ]
  newtracks = [jbrowse.create_new_track(**r) for r in jbrowse_resources]

  jbrowse.backup_file(tracklist_file)
  jbrowse.update_tracklist(tracklist_file,newtracks)

  ''' Pickle pool_summary '''
  print >>sys.stderr, "[ Writing results ]"
  with open('%s/pool_summary.pickle' % pool_path,'wb') as outh:
    pickle.dump(calls,outh)

  ''' Write to JSON format '''
  print >>sys.stderr, "[ Writing DB format ]"
  jobj = []
  for ref,call in calls.iteritems():
    jobj.append({'ref':ref, 'clone_call':call.serializable()})
  
  with open('%s/pool_summary.json' % pool_path,'wb') as outh:
    json.dump(jobj,outh)

  sys.exit(0)