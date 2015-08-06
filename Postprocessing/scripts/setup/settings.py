XML_START = '''<?xml version="1.0" encoding="utf-8"?>
<smrtpipeSettings>
'''

XML_PROTOCOL = ''' <protocol>
  <param name="reference">
   <value>%(refpath)s</value>
  </param>
 </protocol>
'''


XML_END = '''</smrtpipeSettings>'''

def dunion(d1,d2):
  return dict(d1.items() + d2.items()) 

def module_xml(module_name,params=None):
  if params is None:
    xmlstr = '  <module name="%s" />\n' % module_name
  else:
    xmlstr = '  <module name="%s">\n' % module_name
    for name,val in params.iteritems():
      xmlstr += '    <param name="%s"><value>%s</value></param>\n' % (name,val)
    xmlstr += '  </module>\n'
  return xmlstr

''' Standard options '''

std_filter_opts = {'minLength':50,'minSubReadLength':50,'readScore':0.75}
roi_filter_opts = {'minFullPasses':2,'minPredictedAccuracy':90}

# Mapping options used by all mapping protocols
std_mapping_opts = { 'maxHits':1,
                     'maxDivergence':30,
                     'minAnchorSize':12,
                     'samBam':'True',
                     'gff2Bed':'True',
                     'placeRepeatsRandomly':'True',
                     'loadPulsesOpts':'bymetric',                     
                   }

''' Parameters
The "align_opts" parameter in the P_Mapping module is passed to pbalign.py
The minLength parameter sets the minimum threshold for aligned read length.
'''
ccs_mapping_opts = dict(std_mapping_opts)
ccs_mapping_opts['align_opts'] = '--seed=1 --minAccuracy=0.75 --minLength=500 --useCcs=denovo --readType=CCS'
ccs_mapping_opts['pulseMetrics'] = 'QualityValue'

sub_mapping_opts = dict(std_mapping_opts)
sub_mapping_opts['align_opts'] = '--seed=1 --minAccuracy=0.75 --minLength=500 --useQuality'
sub_mapping_opts['pulseMetrics'] = 'QualityValue' # 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'

# Consensus options
std_consensus_opts = { 'algorithm':'plurality',
                       'outputConsensus':'True',
                       'makeVcf':'True',
                       'makeBed':'True',
                       'enableMapQVFilter':'True',
                     }


### Production CCS
def _production_ccs(**kwargs):
  xmlstr =  module_xml("P_Fetch")                               # P_Fetch
  xmlstr += module_xml("P_Filter",std_filter_opts)              # P_Filter
  xmlstr += module_xml("P_FilterReports")                       # P_FilterReports
  xmlstr += module_xml("P_Mapping",dunion(ccs_mapping_opts,kwargs))            # P_FilterReports
  xmlstr += module_xml("P_MappingReports")                      # P_MappingReports
  xmlstr += module_xml("P_GenomicConsensus",std_consensus_opts) # P_MappingReports
  xmlstr += module_xml("P_ConsensusReports")                    # P_MappingReports  
  return xmlstr

### Production subreads
def _production_sub(**kwargs):
  xmlstr =  module_xml("P_Fetch")                               # P_Fetch
  xmlstr += module_xml("P_Filter",std_filter_opts)              # P_Filter
  xmlstr += module_xml("P_FilterReports")                       # P_FilterReports
  xmlstr += module_xml("P_Mapping",dunion(sub_mapping_opts,kwargs))            # P_FilterReports
  xmlstr += module_xml("P_MappingReports")                      # P_MappingReports
  xmlstr += module_xml("P_GenomicConsensus",std_consensus_opts) # P_MappingReports
  xmlstr += module_xml("P_ConsensusReports")                    # P_MappingReports  
  return xmlstr

def _reads_of_insert(**kwargs):
  xmlstr =  module_xml("P_Fetch")                               # P_Fetch
  xmlstr += module_xml("P_CCS",roi_filter_opts)                 # P_Filter
  xmlstr += module_xml("P_Mapping",dunion(ccs_mapping_opts,kwargs))
  xmlstr += module_xml("P_MappingReports")                      # P_MappingReports
  xmlstr += module_xml("P_GenomicConsensus",std_consensus_opts)
  xmlstr += module_xml("P_ConsensusReports")
  return xmlstr

def generate_settings(refpath,protocol,**kwargs):
  xmlstr =  XML_START
  xmlstr += XML_PROTOCOL % {'refpath':refpath}
  if protocol == 'ccs':
    xmlstr += _production_ccs(**kwargs)
  elif protocol == 'sub':
    xmlstr += _production_sub(**kwargs)
  elif protocol == 'roi':
    xmlstr += _reads_of_insert(**kwargs)  
  else:
    return None
  xmlstr += XML_END
  return xmlstr


