#! /usr/bin/env python

#--- import modules ---#
has_scipy = True
try:
  import scipy
  from scipy.stats import poisson
except ImportError as e:
  has_scipy = False
  import sys
  sys.exit("%s" % e)

import itertools

from variant import Variant

#--- utility functions ---#

def intervals(i):
  ''' Combines lists down to continuous intervals '''
  for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
    b = list(b)
    yield b[0][1], b[-1][1]

def remove_overlap(intervals1,intervals2):
  ''' Removes intervals in 1 that overlap with intervals in 2
  '''
  new_intervals = []
  for iv in intervals1:
    found = False
    for niv in intervals2:
      found |= niv[0] <=  iv[0] <= niv[1]
      found |= niv[0] <=  iv[1] <= niv[1]
      found |=  iv[0] <= niv[0] <=  iv[1]
      found |=  iv[0] <= niv[1] <=  iv[1]      
    if not found:
      new_intervals.append(iv)
  return new_intervals

#--- covdepth functions ---#

def parse_covdepth_gatk(infile,as_list=True):
  ''' Parses covdepth output file from GATK
      Returns covdata[reference_name] -> {pos1:cov1,pos2:cov2,...}
      If as_list = True, dictionary is converted to list of coverages with the list index
      as the position, i.e., covdata[reference_name] -> [-1,cov1,cov2,cov3,...covN]
  '''
  refdict = {}
  with open(infile,'rU') as fh:
    header = fh.next()
    for l in fh:
      locus,d1,d2,d3 = l.split('\t')[:4]
      ref,pos = locus.split(':')
      if ref not in refdict:
        refdict[ref] = {}
      refdict[ref][int(pos)] = int(d1)
  if not as_list: return refdict
  # change dict to list
  covdata = {}
  for ref,seq in refdict.iteritems():
    assert max(seq.keys()) == len(seq.keys()), "Different max position (%d) and number of values (%d)" % (max(seq.keys()),len(seq.keys()))
    covdata[ref] = [-1] + [0] * max(seq.keys())
    for pos,depth in seq.iteritems():
      covdata[ref][pos] = depth
  return covdata

def parse_covdepth_samtools(infile,reflens,as_list=True):
  ''' Parses covdepth output file from GATK
      Returns covdata[reference_name] -> {pos1:cov1,pos2:cov2,...}
      If as_list = True, dictionary is converted to list of coverages with the list index
      as the position, i.e., covdata[reference_name] -> [-1,cov1,cov2,cov3,...covN]
  '''
  refdict = {}
  with open(infile,'rU') as fh:
    for l in fh:
      ref,pos,d1 = l.split('\t')
      #ref,pos = locus.split(':')
      if ref not in refdict:
        refdict[ref] = {}
      refdict[ref][int(pos)] = int(d1)
  if not as_list: return refdict
  # change dict to list
  covdata = {}
  for ref,reflen in reflens.iteritems():
    covdata[ref] = [-1] + [0] * reflen
    if ref in refdict:
      for pos,depth in refdict[ref].iteritems():
        covdata[ref][pos] = depth
  return covdata

def local_coverage_score(covlist,window_size=11):
  ''' find local coverage dips
      Returns list of scores for each position in covlist
  '''
  localmeans = [scipy.mean(covlist[max(1,i-window_size):min(i+window_size+1,len(covlist))]) for i,cov in enumerate(covlist)]
  pvals      = [poisson.cdf(cov,localmeans[i]) for i,cov in enumerate(covlist)]
  pvals[0]   = 1
  scores     = [(-10 * scipy.log10(pv)) if pv != 0 else 0 for pv in pvals]
  return scores,localmeans

def adjusted_coverage_score(covlist,window_size=11):
  scores     = []
  localmeans = []
  for i,cov in enumerate(covlist):
    if i==0:
      localmeans.append(0)
      scores.append(0)
      continue
    localmean = scipy.mean(covlist[max(1,i-window_size):min(i+window_size+1,len(covlist))])
    if localmean < 1000:
      pval = poisson.cdf(cov,localmean)
      score = (-10 * scipy.log10(pval)) if pval > 0 else 0
    else:
      if cov / localmean < 0.8: score = 255
      else: score = 0
    localmeans.append(localmean)
    scores.append(score)    
    #print '%d\t%f\t%f\t%f' % (i,cov,localmean,score)
  return scores,localmeans


def summarize_coverage(covlist,min_cov=5): #,seq,chrom,min_cov=5,min_score=30,exclude_edges=False,exclude_overlaps=False):
  ''' identify coverage variants in covlist
      Returns dict with keys 'mean_cov','pct_cov', and 'variants', where dict['variants']
      is a list of Variant objects
  '''
  assert min(covlist[1:]) >= 0
  nocov = [i for i,v in enumerate(covlist) if v < min_cov]
  nocov.remove(0) # take off the -1 at index 0
  mean_cov = scipy.mean(covlist[1:])  
  pct_cov  = 1 - (float(len(nocov)) / (len(covlist) - 1))
  return (pct_cov,mean_cov)

def find_nocov_variants(covlist,chrom='',caller='',min_cov=5):
  variants = []
  assert min(covlist[1:]) >= 0
  nocov = [i for i,v in enumerate(covlist) if v < min_cov]
  nocov.remove(0) # take off the -1 at index 0
  if len(covlist)-1 == len(nocov): return None # entire sequence has no coverage
  nocov_intervals = list(intervals(nocov))
  for iv in nocov_intervals:
    data = {'chrom':chrom,'caller':caller,'pos':iv[0], 'type': 'no_cov'}
    data['length'] = iv[1] - iv[0] + 1
    data['mean_cov'] = scipy.mean(covlist[iv[0]:(iv[1]+1)])
    variants.append(Variant.from_dict(data))
  return variants

def find_variants(covlist,seq,chrom,min_cov=5,min_score=30,exclude_edges=False,exclude_overlaps=False):
  ''' identify coverage variants in covlist
      Returns dict with keys 'mean_cov','pct_cov', and 'variants', where dict['variants']
      is a list of Variant objects
  '''
  assert min(covlist[1:]) >= 0
  assert len(covlist) - 1 == len(seq), "Number of coverage values (%d) is not equal to sequence length (%d)" % (len(covlist)-1,len(seq))
  retval = {}
  nocov = [i for i,v in enumerate(covlist) if v < min_cov]
  nocov.remove(0)
  retval['mean_cov'] = scipy.mean(covlist[1:])  
  retval['pct_cov'] = 1 - (float(len(nocov)) / (len(covlist) - 1))
  if len(nocov) == len(seq):
    return retval
  nocov_intervals = list(intervals(nocov))
  #covscores,localmeans = local_coverage_score(covlist)
  covscores,localmeans = adjusted_coverage_score(covlist)
  covdip = [i for i,v in enumerate(covscores) if v >= min_score]
  covdip_intervals = list(intervals(covdip))
  
  # refine intervals
  if exclude_edges:
    # ignore intervals that overlap the beginning and end of reference
    covdip_intervals = [iv for iv in covdip_intervals if not iv[0]==1 and not iv[1]==(len(covlist)-1)]
  if exclude_overlaps:
    # ignore covdip intervals that overlap with nocov intervals
    covdip_intervals = remove_overlap(covdip_intervals,nocov_intervals)
    # covdip = list(itertools.chain(*[range(v1,v2+1) for v1,v2 in covdip_intervals])
  
  # positions with no coverage are not considered to be coverage dips 
  covdip = [p for p in covdip if p not in nocov]
  
  variants = []
  for iv in nocov_intervals:
    data = {'chrom':chrom, 'pos':iv[0], 'type': 'no_cov'}
    data['length'] = iv[1] - iv[0] + 1
    data['mean_cov'] = scipy.mean(covlist[iv[0]:(iv[1]+1)])
    variants.append(Variant.from_dict(data))
  
  for iv in covdip_intervals:
    data = {'chrom':chrom, 'pos':iv[0], 'type': 'cov_dip'}
    data['length'] = iv[1] - iv[0] + 1
    data['mean_cov'] = scipy.mean(covlist[iv[0]:(iv[1]+1)])
    intscores = covscores[iv[0]:(iv[1]+1)]
    intmeans  = localmeans[iv[0]:(iv[1]+1)]
    data['quality'] = max(intscores)
    data['info'] = {'CovScores':'%s' % ','.join(['%d' % int(round(v)) for v in intscores]),
                    'LocalMeans':'%s' % ','.join(['%d' % int(round(v)) for v in intmeans]),
                   }
    data['ref'] = str(seq[iv[0]:(iv[1]+1)].seq).upper()
    # data['alt'] = data['ref'].lower()    
    variants.append(Variant.from_dict(data))
  
  if variants:
    retval['variants'] = variants
  return retval
