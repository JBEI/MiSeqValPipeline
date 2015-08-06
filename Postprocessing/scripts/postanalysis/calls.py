#--- CloneCall class ---#

from variant import Variant
import json

class CloneCall:
  def __init__(self,call='nocall',mean_cov=None,pct_cov=None,protocols=None):
    self.set_call(call,protocols)
    self.variants = None
    self.mean_cov = mean_cov
    self.pct_cov  = pct_cov

  def set_call(self,call,protocols=None):
    self.call = call
    self.protocols = protocols if protocols is not None else list()
    return self

  def get_reason(self):
    if self.call == 'nocall': return 'not called'
    if self.call == 'lowcov': return '%s have low coverage (%.2f)' % (' and '.join(self.protocols),self.mean_cov)
    if self.call == 'incomplete': return '%s alignments are incomplete (%.3f covered)' % (' and '.join(self.protocols),self.pct_cov)
    if self.call == 'errors': return '%d variants were found by %s' % (len(self.variants),' and '.join(self.protocols))
    if self.call == 'dips': return '%d coverage dips were found by %s' % (len(self.variants),' and '.join(self.protocols))
    if self.call == 'flawless': return 'perfect clone'
    return ''

  def serializable(self):
    jobj = {'call':self.call,'mean_cov':self.mean_cov,'pct_cov':self.pct_cov,'reason':self.get_reason()}
    if self.call == 'almost':
      jobj['fix'] = {'score':self.fixscore,'primer1':self.primer1,'primer2':self.primer2}
    return jobj

  def format_db(self):
    return '%s\t%.3f\t%.3f\t%s' % (self.call,self.mean_cov,self.pct_cov,self.get_reason())

  def format_fix(self):
    return '%d\t%s\t%s' % (self.fixscore,self.primer1,self.primer2)
    
  def __str__(self):
    return '%s\t%.2f\t%.3f\t%s' % (self.call,self.mean_cov,self.pct_cov,self.protocols)
