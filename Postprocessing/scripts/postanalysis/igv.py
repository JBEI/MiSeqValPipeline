from cStringIO import StringIO

''' Create IGV sessions '''
igv_session_header = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="%(analysis_url)s/references/seqs/sequence/seqs.fasta" locus="All" version="5">
    <Resources>'''

igv_session_footer = '''    </Resources>
</Session>
'''

igv_resource_xml = '''        <Resource path="%(path)s" name="%(name)s"/>'''

def create_igv_session(analysis_url,resources):
  igvxml = StringIO()
  print >>igvxml, igv_session_header % {'analysis_url':analysis_url}
  for r in resources:
    print >>igvxml, igv_resource_xml % r
  print >>igvxml, igv_session_footer
  return igvxml
