#--- Variant class ---#
class Variant:
  SFIELDS = ['chrom','type','caller','ref','alt']
  IFIELDS = ['pos','length'] # 'end_pos'
  FFIELDS = ['quality','mean_cov']
  DFIELDS = ['info']
  def __init__(self,data):
    for f in Variant.SFIELDS:
      if f in data: setattr(self,f,data[f])
      else:         setattr(self,f,None)
    for f in Variant.IFIELDS:
      if f in data: setattr(self,f,data[f])
      else:         setattr(self,f,None)
    for f in Variant.FFIELDS:
      if f in data: setattr(self,f,data[f])
      else:         setattr(self,f,None)
    for f in Variant.DFIELDS:
      if f in data: setattr(self,f,data[f])
      else:         setattr(self,f,{})
    # check required fields
    assert self.chrom and self.pos and self.type, "ERROR: missing required fields"
  
  @classmethod
  def from_vcf(cls,line):
    chrom,pos,id,ref,alt,qual,filter,info,format,rgid = line.split('\t')[:10]
    infodict = dict([(kv.split('=')) for kv in info.split(';') if len(kv.split('='))==2])
    data = {'chrom':chrom, 'pos':int(pos)}
    altlen = max([len(a) for a in alt.split(',')])
    if len(ref) == altlen:  data['type'] = 'sub'
    elif len(ref) > altlen: data['type'] = 'del'    
    elif len(ref) < altlen: data['type'] = 'ins'    
    data['length'] = max(len(ref),altlen)
    data['quality'] = float(qual)
    data['ref'] = ref
    data['alt'] = alt
    if id.strip().strip('.'): infodict['id'] = id
    data['info'] = infodict
    if 'DP' in infodict:
      data['mean_cov'] = float(infodict['DP'])
    return cls(data)
  
  @classmethod
  def from_gff(cls,line):
    chrom,f0,type,spos,epos,f1,f2,f3,info = line.split('\t')
    infodict = dict([(kv.split('=')) for kv in info.split(';')])    
    data = {'chrom':chrom,'pos':int(spos)}
    data['type'] = type[:3]
    if data['type']=='ins': # length is calculated differently for insertions
      data['length'] = len(infodict['variantSeq']) #int(infodict['length'])
    else:
      data['length'] = int(epos) - int(spos) + 1
    ### assert int(infodict['length']) == int(epos) - int(spos) + 1, "%d != %d" % (int(infodict['length']),int(epos) - int(spos) + 1)
    data['quality'] = float(infodict['confidence'])
    if 'reference' in infodict:
      data['ref'] = infodict['reference']
    if 'variantSeq' in infodict:
      data['alt']  = infodict['variantSeq']
    if 'coverage' in infodict:
      data['mean_cov']  = float(infodict['coverage'])
    data['info'] = infodict
    return cls(data)
  
  @classmethod
  def from_dict(cls,data):
    return cls(data)
  
  def out(self):
    outs = ''
    for f in Variant.SFIELDS:
      if getattr(self,f) is not None: outs += '%s\t' % getattr(self,f)
      else: outs += '.\t'
      #else:         outs += '.\t' #setattr(self,f,None)
    for f in Variant.IFIELDS:
      if getattr(self,f) is not None: outs += '%d\t' % getattr(self,f)
      else: outs += '.\t'
    for f in Variant.FFIELDS:
      if getattr(self,f) is not None: outs += '%.1f\t' % getattr(self,f)
      else: outs += '.\t'
    for f in Variant.DFIELDS:
      if getattr(self,f) is not None: 
        outs += ';'.join(['%s=%s' % (k,v) for k,v in getattr(self,f).iteritems()])
      else: outs += '.'
      outs += '\t'
    if outs.endswith('\t'): outs = outs[:-1]
    return outs

  def simpleVCF(self,typeAsID=False,altAsN=False):
    outs = []
    outs.append(self.chrom)        # CHROM
    outs.append('%d' % self.pos)   # POS
    #--- ID
    if typeAsID:
      outs.append('%s.%d' % (self.type,self.pos))
    else:
      outs.append('.')
    #---REF
    if self.ref: outs.append('%s' % self.ref)
    else: outs.append('.')
    #---ALT
    if self.alt:
      outs.append(self.alt)
    else:
      outs.append('.')      
    #---QUALITY
    if self.quality: outs.append('%.1f' % self.quality)
    else: outs.append('.')
    #--- FILTER    
    outs.append('.')
    #--- INFO
    infokv = self.info
    if self.mean_cov is not None:
      infokv['DP'] = '%.1f' % self.mean_cov
    if self.type is not None:
      infokv['TYPE'] = '%s' % self.type
    if self.caller is not None:
      infokv['CALLER'] = '%s' % self.caller
    if infokv:
      outs.append(';'.join(['%s=%s' % (k,v) for k,v in infokv.iteritems()]))
    else:
      outs.append('.')
    return '\t'.join(outs)
  
  def pacbioGFF(self):
    outs = []
    outs.append(self.chrom)               # seqid
    outs.append('.')                      # source
    # type
    if self.type is not None:
      if self.type == 'del':   outs.append('deletion')
      elif self.type == 'ins': outs.append('insertion')
      else:                    outs.append('substitution')
    outs.append('%d' % self.pos)          # start
    outs.append('%d' % (self.pos + self.length - 1))  # end
    outs.append('.')                      # score
    outs.append('.')                      # strand
    outs.append('.')                      # phase
    # attributes
    attr_dict = {}
    attr_dict['reference'] =  self.ref
    attr_dict['variantSeq'] =  self.alt
    attr_dict['coverage'] =  '%d' % int(round(self.mean_cov))
    if self.quality is not None:
      attr_dict['confidence'] =  '%d' % int(round(self.quality))
    else:
      attr_dict['confidence'] =  '%d' % 1
    attr_dict['frequency'] =  '%d' % int(round(self.mean_cov))
    attr_dict['length'] =  '%d' % self.length
    if attr_dict:
      attributes = ';'.join(['%s=%s' % kv for kv in attr_dict.iteritems()])
    else:
      attributes = '.'    
    outs.append(attributes)
    return '\t'.join(outs)

  def simpleGFF(self):
    outs = []
    outs.append(self.chrom)               # seqid
    outs.append('cov')                    # source
    outs.append(self.type)                # type
    outs.append('%d' % self.pos)          # start
    outs.append('%d' % (self.pos + self.length))  # end
    # score
    if self.quality is not None:
      outs.append('%.2f' % self.quality)
    else:
      outs.append('.')
    outs.append('.')                      # strand
    outs.append('.') # phase
    # attributes
    attr_dict = {}
    if self.type=='cov_dip': attr_dict['Description'] = 'Coverage dip'
    elif self.type=='no_cov': attr_dict['Description'] = 'No coverage'
    if self.info is not None:    
      attr_dict.update(self.info)
    if attr_dict:
      attributes = ';'.join(['%s=%s' % kv for kv in attr_dict.iteritems()])
    else:
      attributes = '.'
    outs.append(attributes)
    return '\t'.join(outs)
  
  def dbstring(self):
    dbfields = []
    if self.type: dbfields.append('type:%s' % self.type)
    if self.caller: dbfields.append('caller:%s' % self.caller)
    if self.ref: dbfields.append('ref:%s' % self.ref)
    if self.alt: dbfields.append('alt:%s' % self.alt)
    if self.pos: dbfields.append('pos:%d' % self.pos)
    if self.length: dbfields.append('length:%d' % self.length)
    if self.quality: dbfields.append('quality:%d' % self.quality)
    if self.mean_cov: dbfields.append('mean_cov:%d' % self.mean_cov)
    return ';'.join(dbfields)
  
  def serializable(self):
    outd = {}
    for f in Variant.SFIELDS:
      if getattr(self,f) is not None: outd[f] = getattr(self,f)
    for f in Variant.IFIELDS:
      if getattr(self,f) is not None: outd[f] = getattr(self,f)
    for f in Variant.FFIELDS:
      if getattr(self,f) is not None: outd[f] = getattr(self,f)
    for f in Variant.DFIELDS:
      if getattr(self,f): outd[f] = getattr(self,f)
    return outd
  
  
  def __str__(self):
    return self.out()
