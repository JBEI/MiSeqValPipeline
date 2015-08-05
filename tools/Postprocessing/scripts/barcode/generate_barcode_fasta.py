#! /usr/bin/env python

PBBARCODES = {'B001':'TCAGACGATGCGTCAT',
              'B002':'CTATACATGACTCTGC',
              'B003':'TACTAGAGTAGCACTC',
              'B004':'TGTGTATCAGTACATG',
              'B005':'ACACGCATGACACACT',
              'B006':'GATCTCTACTATATGC',
              'B007':'ACAGTCTATACTGCTG',
              'B008':'ATGATGTGCTACATCT',
              'B009':'CTGCGTGCTCTACGAC',
              'B010':'GCGCGATACGATGACT',
              'B011':'CGCGCTCAGCTGATCG',
              'B012':'GCGCACGCACTACAGA',
              'B013':'ACACTGACGTCGCGAC',
              'B014':'CGTCTATATACGTATA',
              'B015':'ATAGAGACTCAGAGCT',
              'B016':'TAGATGCGAGAGTAGA',
              'B017':'CATAGCGACTATCGTG',
              'B018':'CATCACTACGCTAGAT',
              'B019':'CGCATCTGTGCATGCA',
              'B020':'TATGTGATCGTCTCTC',
              'B021':'GTACACGCTGTGACTA',
              'B022':'CGTGTCGCGCATATCT',
              'B023':'ATATCAGTCATGCATA',
              'B024':'GAGATCGACAGTCTCG',
              'B025':'CACGCACACACGCGCG',
              'B026':'CGAGCACGCGCGTGTG',
              'B027':'GTAGTCTCGCACAGAT',
              'B028':'GAGACTCTGTGCGCGT',
              'B029':'GCTCGACTGTGAGAGA',
              'B030':'AGAGATGTGTGATGAC',
              'B031':'TACGACTACATATCAG',
              'B032':'TATCTCTGTAGAGTCT',
              'B033':'AGAGAGAGACATGCGC',
              'B034':'ACTCTCGCTCTGTAGA',
              'B035':'TCTATGTCTCAGTAGT',
              'B036':'GCGTATATCTCATGCG',
              'B037':'GTGCGTATGTCGCTAC',
              'B038':'TGCTCGCAGTATCACA',
              'B039':'CTGTGTGTGATAGAGT',
              'B040':'CAGTGAGAGCGCGATA',
              'B041':'GTACATATGCGTCTGT',
              'B042':'GAGACTAGAGATAGTG',
              'B043':'TACGCGTGTACGCAGA',
              'B044':'TGTCACTCATCTGAGT',
              'B045':'GCACATACACGCTCAC',
              'B046':'GCTCGTCGCGCGCACA',
              'B047':'ACAGTGCGCTGTCTAT',
              'B048':'TCACACTCTAGAGCGA',
              'B049':'TCACATATGTATACAT',
              'B050':'CGCTGCGAGAGACAGT',
              'B051':'ACACACAGACTGTGAG',
              'B052':'GCAGACTCTCACACGC',
              'B053':'TGCTCTCGTGTACTGT',
              'B054':'GTGTGAGATATATATC',
              'B055':'CTCAGTGTGACACATG',
              'B056':'TGCGAGCGACTCTATC',
              'B057':'GTCAGCTAGTGTCAGC',
              'B058':'AGATATCATCAGCGAG',
              'B059':'GTGCAGTGATCGATGA',
              'B060':'TGACTCGCTCATAGTC',
              'B061':'ATGCTGATGACGCGCT',
              'B062':'GACAGCATCTGCGCTC',
              'B063':'AGCGTCTGACGTGAGT',
              'B064':'TCGATATACGACGTGC',
              'B065':'TCGTCATACGCTCTAG',
              'B066':'CGACTACGTACAGTAG',
              'B067':'GCGTAGACAGACTACA',
              'B068':'ACAGTATGATGTACTC',
              'B069':'GTCTGATAGATACAGA',
              'B070':'CTGCGCAGTACGTGCA',
              'B071':'TAGATCTCTGACTCAC',
              'B072':'CTGATGCGCGCTGTAC',
              'B073':'CACTCGTGCACGATGC',
              'B074':'TGACAGTATCACAGTG',
              'B075':'GAGATACGCTGCAGTC',
              'B076':'ACGTGAGCTCACTCGC',
              'B077':'ATAGAGAGTGTCTCAG',
              'B078':'CATAGAGAGATAGTAT',
              'B079':'ATCTCGAGATGTAGCG',
              'B080':'ACGATCACTCGTGTCA',
              'B081':'GATCGACTCGAGCATC',
              'B082':'ATGCTCACTACTACAT',
              'B083':'CGTGCACATCTATAGC',
              'B084':'GACTGCACATGCACGA',
              'B085':'TATGACTAGTGTACTA',
              'B086':'GACGTGTCGTAGATAT',
              'B087':'ATAGCGACGCGATATA',
              'B088':'ATCGCTGTGTCTATAG',
              'B089':'TCTCACTGATAGCGTG',
              'B090':'TGTCGTCTATCATGTA',
              'B091':'CACACGAGATCTCATC',
              'B092':'AGATACACATGATACT',
              'B093':'CGTGAGTAGTCAGACG',
              'B094':'TCTCGACTGCACATAT',
              'B095':'TGAGTGACGTGTAGCG',
              'B096':'GTGTGCACTCACACTC',
              'B097':'TACGATCGTAGCTGCT',
              'B098':'TATACACACTCGCTCG',
              'B099':'AGCGCTGCGACACGCG',
              'B100':'GTCGTAGCTGCTGTAT',
              'B101':'CTGTACTAGAGCGTCT',
              'B102':'TCGAGTGTATAGCTCA',
              'B103':'ACTGTGACAGTATGAT',
              'B104':'TGTCTGAGACGCATAC',
              'B105':'CACTCACGTGTGATAT',
              'B106':'ATCGCATCGCAGAGAC',
              'B107':'TACTCATATATGCTAC',
              'B108':'GTCTACGCTCGTCGCG',
              'B109':'TGCGAGACTATCGCGA',
              'B110':'CAGATCTCTCTGATGT',
              'B111':'GTAGAGTGATCGCGTC',
              'B112':'ACGACAGTCAGAGTAT',
              'B113':'ATATATAGCTGATGCG',
              'B114':'TGCTATCTGAGATACT',
              'B115':'CAGCAGATCATGTCGA',
              'B116':'TGCTGCGAGCGCTCTG',
              'B117':'ACTATCGCAGCTCAGT',
              'B118':'CGTCTCTCGTCTGTGC',
              'B119':'GAGTCTCGATATACTA',
              'B120':'TGTCATGTGTACACAC',
              'B121':'TCTGTCGATATACACT',
              'B122':'ACGTGCTCTATAGAGA',
              'B123':'TATCAGCACGACATGC',
              'B124':'GCTCTCACGATATCAG',
              'B125':'TATATGCTCTGTGTGA',
              'B126':'GATAGCTGCTAGCTGA',
              'B127':'TCTCATGTGTGAGCTA',
              'B128':'TCAGATGTGTCGCGAG',
              'B129':'CGTAGCTCAGACACTC',
              'B130':'TCAGAGACACTACGAG',
              'B131':'ATCGAGCAGCAGTCGT',
              'B132':'CGTAGCTCGAGATGAG',
              'B133':'GCTAGTCGATGACAGC',
              'B134':'CATGATGCGAGACGCT',
              'B135':'GTGTAGCGTAGACAGT',
              'B136':'AGCACGTGTGTCGACA',
              'B137':'CTAGACACGCAGTCAC',
              'B138':'TAGCGTGAGAGTGTCG',
              'B139':'GTCTCTCTCTCACGCA',
              'B140':'TGCATAGTAGTGCTCT',
              'B141':'CATATCAGTGCTACAG',
              'B142':'CGACGTCATAGTGCGT',
              'B143':'ACACACTCTATCAGAT',
              'B144':'GCTGTGTGTGCTCGTC',
              'B145':'AGCGTAGCATCTGAGC',
              'B146':'GAGTCTGCACGCGCTA',
              'B147':'AGACGCGAGCGCGTAG',
              'B148':'CTACGATGCTATGTAT',
              'B149':'CGACTAGATCTATCAT',
              'B150':'ATCTCTGTGCGCGCAG',
              'B151':'GCTAGCATGCTCTCAG',
              'B152':'GTCACGATATAGTGAC',
              'B153':'TCTACTGCATGATGTC',
              'B154':'AGTCGTGACTATGCTC',
              'B155':'GTATAGACAGATGTGC',
              'B156':'TAGTGTGCGACTCTGA',
              'B157':'GCACTCAGAGACGCGA',
              'B158':'TCTATCAGCGCTGATG',
              'B159':'ATGTCGCATATATCGC',
              'B160':'CACGACTATATGCTCT',
              'B161':'AGTCACACGCACGCTG',
              'B162':'CATACATCGCGCAGTA',
              'B163':'TGCGAGCGTGCACAGA',
              'B164':'CTCTGACTCGCGTCGA',
              'B165':'CTATCTAGCACTCACA',
              'B166':'ACACGTGATAGCTACG',
              'B167':'GCGATCACTGTACACT',
              'B168':'CGCTAGAGATCTGCTA',
              'B169':'GATACTGACACACTAT',
              'B170':'GAGCTGATGTACATGT',
              'B171':'AGTCGCGTAGCTCATC',
              'B172':'TGTAGAGATACTCACT',
              'B173':'TCGCTGACTCGACACA',
              'B174':'TACATCTCGCTGCGCA',
              'B175':'GTATATATATACGTCT',
              'B176':'TCGCGAGCAGCGACAT',
              'B177':'AGCTCAGTATCATCTG',
              'B178':'ACACAGTAGAGCGAGC',
              'B179':'ACGACGCGCACTGACA',
              'B180':'CTCATAGCGTGTACTC',
              'B181':'GACGACAGACTGCATA',
              'B182':'GTCTGTATAGCTATCT',
              'B183':'TGTCTCGTGCTGAGAC',
              'B184':'CATATGCTCGTGCACT',
              'B185':'ACTACATACTAGATCA',
              'B186':'TGTGCACGACAGCAGT',
              'B187':'ATGATACACGCGCGAC',
              'B188':'TGTCTGATCTGTATCA',
              'B189':'CTCTCGCATACGCGAG',
              'B190':'GAGCGTGTATACAGCG',
              'B191':'GAGCTCATGTAGACAC',
              'B192':'TACATATGTCACGCGC',
              'B193':'ATCGCTCTCATGTCTA',
              'B194':'ACGATGTATCTACGCA',
              'B195':'TCGATACGCACTCGAT',
              'B196':'CACGACACGACGATGT',
              'B197':'CTGCAGCTCACTACTA',
              'B198':'CTATATGAGACGAGTG',
              'B199':'CTCTCGTAGACAGATA',
              'B200':'CGCATGACACGTGTGT',
              'B201':'CACATACTACTACTGA',
              'B202':'AGTCAGATGCGCACTC',
              'B203':'AGCGACGCGAGAGTGC',
              'B204':'ATACACTCATGTGCAC',
              'B205':'GCTACGCTATAGACAT',
              'B206':'TATCTATCGCATATCG',
              'B207':'TCACGTGCAGATATAG',
              'B208':'GCACAGCGTAGCGCAT',
              'B209':'CATGCTACGTCTCTGT',
              'B210':'CTCACGTACGTCACAC',
              'B211':'TCTGAGACACAGACTC',
              'B212':'CTAGTCTCTATCGCAT',
              'B213':'ACGCTCGCTGAGCATA',
              'B214':'ACTCATGTATATGAGT',
              'B215':'AGCGTAGCGCGCGTCA',
              'B216':'TCTCGTCGCAGTCTCT',
              'B217':'GACGAGCGTCTGAGAG',
              'B218':'GTATGATCACTAGTAG',
              'B219':'CTCACACATACACGTC',
              'B220':'GTATCGAGCGTATAGC',
              'B221':'GCTGCGCTGATATGCG',
              'B222':'GTCAGAGCTCTCGTGC',
              'B223':'ATATGACATACACGCA',
              'B224':'CTCGCTCGACGAGCGC',
              'B225':'CGTCATCTATATACAG',
              'B226':'TGTACGCTCTCTATAT',
              'B227':'AGATCGCGCATGTGTA',
              'B228':'GACACAGTGTGTAGTC',
              'B229':'GTGCGCTACAGTCTCT',
              'B230':'CATCGTCTAGCACTCG',
              'B231':'CAGCGCATCTCACGTC',
              'B232':'GTCTCATCATGCTGCG',
              'B233':'ATCGTATAGTCATACA',
              'B234':'AGTGCGCACATGTCAG',
              'B235':'ATCTACGACTAGCAGA',
              'B236':'TCGCGACATATAGATG',
              'B237':'AGATATACTGTCTGAT',
              'B238':'AGTCACTGTCTACTCG',
              'B239':'TATACGAGATACGTGA',
              'B240':'ACATGCGTGACAGTCA',
              'B241':'GTGAGAGTCTGATACT',
              'B242':'GCACGATGTCAGCGCG',
              'B243':'CACGTGCTCGAGAGTC',
              'B244':'GACACTCAGTCTCTCA',
              'B245':'ACAGTAGACTCTCAGA',
              'B246':'ACACTAGATCGCGTGT',
              'B247':'ACGTCAGCACTGCTCT',
              'B248':'CACAGTCGCAGTACGC',
              'B249':'GTGACTCTATGCTATA',
              'B250':'CTCTACATCAGTGCTA',
              'B251':'GATGAGTATAGACACA',
              'B252':'ATCTGAGTCTGACACG',
              'B253':'GCGAGACTCAGCTCTG',
              'B254':'CGTACGACTGCAGCGT',
              'B255':'CGTGTCACTCTGCGTG',
              'B256':'AGCTCTGTCACTAGAC',
              'B257':'GCGAGAGTGAGACGCA',
              'B258':'TCTACTACACTGTACT',
              'B259':'CATCGTCACAGACATA',
              'B260':'GTGCACTCGCGCTCTC',
              'B261':'TGACATCTACACATAC',
              'B262':'GTCGTCTAGATCGACG',
              'B263':'GACATAGCTAGATCGC',
              'B264':'TATATATGTCTATAGA',
              'B265':'CTGTGTATCTGTGTAC',
              'B266':'CGACGCACGATACTAT',
              'B267':'TGATATATACGCGCGT',
              'B268':'CGCGTATGTATGTCGC',
              'B269':'CTCGAGCAGTAGATAC',
              'B270':'CTGTGCTATGTACGCG',
              'B271':'ACTCAGCGCGTACATA',
              'B272':'TGAGATATGCATGATG',
              'B273':'ACTCTATGTCGATGTA',
              'B274':'GCGCGTGCTGCGTCTA',
              'B275':'GATCATGTGAGCATAG',
              'B276':'CATGTAGAGCAGAGAG',
              'B277':'GTGTGTCTCGATGCGC',
              'B278':'CTCGCACGTCGCATAG',
              'B279':'CGAGCTACTCTGACAG',
              'B280':'CGTGAGTATATGTCAT',
              'B281':'ACAGTACTAGTGCGAG',
              'B282':'CTCACTACGCGCGCGT',
              'B283':'GACTCTCTATCGTACT',
              'B284':'TATATACAGAGTCGAG',
              'B285':'TGAGTGAGACATATCA',
              'B286':'GTGACACACAGAGCAC',
              'B287':'CTGCGTATAGATATGA',
              'B288':'GAGAGTGTGAGAGTGT',
              'B289':'CGTCTCTATCTCTCTA',
              'B290':'TACATGTGTCTATGTC',
              'B291':'TCTCGCGCGTGCACGC',
              'B292':'TATGTGTCTGCGCATA',
              'B293':'AGTCTGAGAGAGCTAT',
              'B294':'ACAGTCGAGCGCTGCG',
              'B295':'GAGAGTAGCGTGTACA',
              'B296':'GATATATCGAGTATAT',
              'B297':'GCACACATATCTGATG',
              'B298':'CATCGCGAGTGCGCTC',
              'B299':'ACATATCGTACTCTCT',
              'B300':'AGCACAGTCACATGTC',
              'B301':'GCGCACAGACATCTGT',
              'B302':'ACGCGCTATCTCAGAG',
              'B303':'CTGTAGACATCACACG',
              'B304':'TATCTGAGCGCGAGCA',
              'B305':'CTCTGCTCTGACTCTC',
              'B306':'ACGTAGTGCACACAGA',
              'B307':'TGTATGAGTGTCTGAC',
              'B308':'CTCTGCAGCGATCACT',
              'B309':'ACTGCGAGATACACAC',
              'B310':'TATAGTGCGCAGCGAC',
              'B311':'GATGTGTGCGCAGTGC',
              'B312':'AGACACACACGCACAT',
              'B313':'CACATGTGACTCGACG',
              'B314':'GATCTGTCGTGAGCGT',
              'B315':'ATATAGCGCATAGCTC',
              'B316':'ACTCATCACGTCTCGA',
              'B317':'CTCTCTAGAGTGACAT',
              'B318':'TCACACTGTGCGAGAC',
              'B319':'CGCGCGAGTATCTCGT',
              'B320':'TATCTCTCGAGTCGCG',
              'B321':'TAGATGAGTACACGTA',
              'B322':'CATGTGCGCTCATCAC',
              'B323':'GTATAGCACTCGAGCG',
              'B324':'ACTCTGCTGTCATCGC',
              'B325':'CGCATATCTCACTAGT',
              'B326':'CACTATACACTGCGCT',
              'B327':'CGCACAGATACGCTCT',
              'B328':'CAGATCTCGCGTGACA',
              'B329':'GCGCTCTCTCACATAC',
              'B330':'ACACATCTCGTGAGAG',
              'B331':'AGTAGTGTGATACTAG',
              'B332':'CGAGCATATATATCTC',
              'B333':'CTATACGTATATCTAT',
              'B334':'GTGTATCAGCGAGTAT',
              'B335':'GCTGAGACGACGCGCG',
              'B336':'GCGCAGTGTCACATCA',
              'B337':'TCATACACACAGATAG',
              'B338':'CACTCGACTCTCGCGT',
              'B339':'CACATATCAGAGTGCG',
              'B340':'CGTATACAGTCACGCT',
              'B341':'TGTAGACTAGCGCTGC',
              'B342':'AGCACACATATAGCGC',
              'B343':'GATATCTCGATCTCTG',
              'B344':'TCTCACGAGAGCGCAC',
              'B345':'TGTGCTCTCTACACAG',
              'B346':'TGTCATATGAGAGTGT',
              'B347':'CTGTGTGCTCGCTATG',
              'B348':'TATAGAGCTCTACATA',
              'B349':'CTATACATAGTGATGT',
              'B350':'TCTCTCTATCGCGCTC',
              'B351':'ATAGCGACATCTCTCT',
              'B352':'GCGCGCGCACTCTCTG',
              'B353':'TCTCTCGATATGATAG',
              'B354':'GATCACAGAGATGCTC',
              'B355':'GCTCGCACAGCGCGTC',
              'B356':'CACAGAGACACGCACA',
              'B357':'GCGTGTGTCGAGTGTA',
              'B358':'GTCATCTGTACGCTAT',
              'B359':'CACACGCACTGAGATA',
              'B360':'ACACATATCGCACTAC',
              'B361':'GAGAGCGCTGACTCTG',
              'B362':'ACACGTGTGCTCTCTC',
              'B363':'CGAGTGTGTCTATACT',
              'B364':'GTGATGCATACGTACA',
              'B365':'CTCGTGACGCTGACTG',
              'B366':'TCTGTATCTCTATGTG',
              'B367':'TGTGTCTCTGAGAGTA',
              'B368':'TAGATCTATCATCGTC',
              'B369':'ACATATACAGCGTATC',
              'B370':'CGCTCATATGAGCTCA',
              'B371':'GTCGCGCATAGAGCGC',
              'B372':'TACACACTATGTGCGT',
              'B373':'ATACGCGCGCGCATGC',
              'B374':'GTGCGCGAGAGTATAC',
              'B375':'GCGCTAGTGTGTACGA',
              'B376':'GAGACACGTCGCACAC',
              'B377':'ACAGAGTGTGCAGATA',
              'B378':'TAGAGCGTCTCTCGTA',
              'B379':'TCTATGAGCACTCTCG',
              'B380':'ATGTGTATATAGATAT',
              'B381':'CTCACACTCTCTCACA',
              'B382':'TCAGCGCACTGTGCTG',
              'B383':'GTGCATACATACATAT',
              'B384':'CAGAGAGATATCTCTG',
             }


def decodeList(s,delim=','):
  import re
  rangeRE = re.compile('^(\d+)-(\d+)$')
  numRE   = re.compile('^(\d+)$')
  numbers = []
  for f in s.split(delim):
    if numRE.match(f):
      numbers.append(int(f))
    else:
      m = rangeRE.match(f)
      assert m is not None
      start = int(m.group(1))
      end   = int(m.group(2))
      if end < start:
        start,end = (end,start)
      numbers.extend(range(start,end+1))
  return sorted(set(numbers))

def main(parser):
  args = parser.parse_args()
  bcnums = decodeList(args.barcodes,args.delim)
  for bn in bcnums:
    seqname = 'B%03d' % bn
    if seqname in PBBARCODES:
      print >>args.outfile, '>%s\n%s' % (seqname,PBBARCODES[seqname])

if __name__ == '__main__':
  import sys
  import argparse
  parser = argparse.ArgumentParser(prog='generate_barcode_fasta', description='Create FASTA file with barcode sequences')
  parser.add_argument('--delim', type=str, help='delimiter used between numbers')
  parser.add_argument('barcodes', type=str, help='delimited list of barcode numbers to use')
  parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
  main(parser)