#! /usr/bin/perl

#FAILS 255 IF perl -T

use Test::More;
use Bio::DNAssemble;
use Bio::Seq;

use strict;
use warnings;

my $DNA = Bio::DNAssemble->new();

unless ($DNA->vmatch)
{
    plan skip_all => 'Vmatch support not installed';
}
else
{
    plan tests => 1;
}

my $seq = "TGACTTGTACCGATGAGCTGGCTCTTCTGGGCGAGCTGGCTGATCTTGACGAGCAGACTTCTCCCGA";
$seq .= "CGAGCTGACTTGTGTCGATGAGCTGGCTCTTCTGGGCGAGTTGGCTGATCTTGACGAGCAGACTTCTCC";
$seq .= "CGACGAGCTGACTTGTGTCGATGAGCTGGCTCTTCTGGGCGAACTGGCTGATCTTGACGAGCAGACTTC";
$seq .= "TCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTA";
$seq .= "TTCCTTATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACG";
$seq .= "CAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGG";
$seq .= "TCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGATGTATCAACATG";
$seq .= "TGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACTTC";
$seq .= "ATT";
my $seqobj = Bio::Seq->new( -seq => $seq, -id => "Dicot_A5.05");

my $rhsh = {
  '276' => 'TCCTTATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTG',
  '443' => 'GCCGGTTGATCACTGAGGCCCAATTGATGTATCAACATGTG',
  '206' => 'TCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCG',
  '358' => 'TCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAG',
  '331' => 'CCCCTTTTCTACGCAGCCGAACTCACATCCTGACCAGGCCG',
  '512' => 'TTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACTTCAT',
  '233' => 'TTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCT',
  '437' => 'AAATGGGCCGGTTGATCACTGAGGCCCAATTGATGTATCAA',
  '463' => 'CAATTGATGTATCAACATGTGGTTTTTATAAAAAGAGTCGT',
  '259' => 'TTCCCGAGACTTTCTATTCCTTATTTATACCCGTCCGTATA',
  '458' => 'AGGCCCAATTGATGTATCAACATGTGGTTTTTATAAAAAGA',
  '424' => 'TCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCA',
  '368' => 'GCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTG',
  '451' => 'ATCACTGAGGCCCAATTGATGTATCAACATGTGGTTTTTAT',
  '220' => 'CTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTT',
  '454' => 'ACTGAGGCCCAATTGATGTATCAACATGTGGTTTTTATAAA',
  '316' => 'GAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCAC',
  '395' => 'CGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTT',
  '492' => 'AAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGT',
  '347' => 'CCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGG',
  '208' => 'CCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAA',
  '511' => 'GTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACTTCA',
  '434' => 'TTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGATGTAT',
  '292' => 'TCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCC',
  '378' => 'GGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGG',
  '325' => 'AGAGTGCCCCTTTTCTACGCAGCCGAACTCACATCCTGACC',
  '350' => 'AACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCC',
  '291' => 'GTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGC',
  '374' => 'CTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTC',
  '442' => 'GGCCGGTTGATCACTGAGGCCCAATTGATGTATCAACATGT',
  '199' => 'AGACTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCA',
  '429' => 'ATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGA',
  '226' => 'CTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAG',
  '211' => 'ACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAA',
  '431' => 'TCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGATG',
  '382' => 'TGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGG',
  '337' => 'TTCTACGCAGCCGAACTCACATCCTGACCAGGCCGGGCTTC',
  '409' => 'ATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTT',
  '340' => 'TACGCAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGC',
  '311' => 'AAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAA',
  '241' => 'GTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTA',
  '198' => 'CAGACTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCC',
  '489' => 'TATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCT',
  '389' => 'CCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTA',
  '327' => 'AGTGCCCCTTTTCTACGCAGCCGAACTCACATCCTGACCAG',
  '495' => 'AAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTT',
  '320' => 'TCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACATCC',
  '418' => 'GCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGA',
  '280' => 'TATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATT',
  '273' => 'TATTCCTTATTTATACCCGTCCGTATAGTAGGGTACGCAAG',
  '471' => 'GTATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAG',
  '236' => 'TCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATT',
  '465' => 'ATTGATGTATCAACATGTGGTTTTTATAAAAAGAGTCGTGA',
  '361' => 'TGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTC',
  '202' => 'CTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGT',
  '218' => 'GACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCT',
  '249' => 'AAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTATACC',
  '348' => 'CGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGG',
  '285' => 'ATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGA',
  '412' => 'GTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGAT',
  '385' => 'TGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTC',
  '234' => 'TCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTA',
  '502' => 'GTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAAT',
  '314' => 'GTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTC',
  '307' => 'ACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGC',
  '388' => 'GCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTT',
  '364' => 'CCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAA',
  '355' => 'ACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTC',
  '479' => 'ATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCT',
  '486' => 'TTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATC',
  '326' => 'GAGTGCCCCTTTTCTACGCAGCCGAACTCACATCCTGACCA',
  '509' => 'GAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACTT',
  '367' => 'GGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGT',
  '335' => 'TTTTCTACGCAGCCGAACTCACATCCTGACCAGGCCGGGCT',
  '270' => 'TTCTATTCCTTATTTATACCCGTCCGTATAGTAGGGTACGC',
  '485' => 'TTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAAT',
  '391' => 'GGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATT',
  '288' => 'CCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAG',
  '460' => 'GCCCAATTGATGTATCAACATGTGGTTTTTATAAAAAGAGT',
  '453' => 'CACTGAGGCCCAATTGATGTATCAACATGTGGTTTTTATAA',
  '324' => 'GAGAGTGCCCCTTTTCTACGCAGCCGAACTCACATCCTGAC',
  '244' => 'TCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTT',
  '433' => 'CTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGATGTA',
  '351' => 'ACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCG',
  '410' => 'TGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTG',
  '240' => 'GGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTT',
  '246' => 'GAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTAT',
  '334' => 'CTTTTCTACGCAGCCGAACTCACATCCTGACCAGGCCGGGC',
  '440' => 'TGGGCCGGTTGATCACTGAGGCCCAATTGATGTATCAACAT',
  '488' => 'TTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCC',
  '230' => 'CCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTT',
  '299' => 'AGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTC',
  '430' => 'TTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGAT',
  '377' => 'CGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGG',
  '447' => 'GTTGATCACTGAGGCCCAATTGATGTATCAACATGTGGTTT',
  '201' => 'ACTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGG',
  '379' => 'GCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGC',
  '423' => 'GTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCC',
  '415' => 'GGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCAC',
  '452' => 'TCACTGAGGCCCAATTGATGTATCAACATGTGGTTTTTATA',
  '342' => 'CGCAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCT',
  '266' => 'GACTTTCTATTCCTTATTTATACCCGTCCGTATAGTAGGGT',
  '295' => 'GTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCT',
  '480' => 'TGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTA',
  '341' => 'ACGCAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCC',
  '438' => 'AATGGGCCGGTTGATCACTGAGGCCCAATTGATGTATCAAC',
  '467' => 'TGATGTATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGA',
  '481' => 'GTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAA',
  '474' => 'TCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTT',
  '444' => 'CCGGTTGATCACTGAGGCCCAATTGATGTATCAACATGTGG',
  '214' => 'AGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTC',
  '508' => 'AGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACT',
  '422' => 'GGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCC',
  '221' => 'TTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTC',
  '417' => 'GGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTG',
  '302' => 'AGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTAC',
  '312' => 'AGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAAC',
  '229' => 'TCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACT',
  '507' => 'AAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAAC',
  '405' => 'AGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCC',
  '260' => 'TCCCGAGACTTTCTATTCCTTATTTATACCCGTCCGTATAG',
  '237' => 'CCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTC',
  '503' => 'TGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATC',
  '370' => 'CGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGAT',
  '309' => 'GCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCG',
  '402' => 'TAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGG',
  '315' => 'TGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCA',
  '506' => 'GAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAA',
  '416' => 'GGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACT',
  '338' => 'TCTACGCAGCCGAACTCACATCCTGACCAGGCCGGGCTTCG',
  '380' => 'CCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCT',
  '222' => 'TGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCC',
  '300' => 'GTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCT',
  '286' => 'TACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAG',
  '484' => 'GTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAA',
  '381' => 'CTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTG',
  '305' => 'GTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCA',
  '308' => 'CGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCC',
  '392' => 'GCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTC',
  '254' => 'CCCCTTTCCCGAGACTTTCTATTCCTTATTTATACCCGTCC',
  '496' => 'AGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTT',
  '373' => 'GCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGT',
  '217' => 'TGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCC',
  '328' => 'GTGCCCCTTTTCTACGCAGCCGAACTCACATCCTGACCAGG',
  '239' => 'AGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCT',
  '464' => 'AATTGATGTATCAACATGTGGTTTTTATAAAAAGAGTCGTG',
  '205' => 'CTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTC',
  '269' => 'TTTCTATTCCTTATTTATACCCGTCCGTATAGTAGGGTACG',
  '281' => 'ATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTC',
  '363' => 'ACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTA',
  '399' => 'TTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAA',
  '235' => 'CTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTAT',
  '301' => 'TAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTA',
  '436' => 'GAAATGGGCCGGTTGATCACTGAGGCCCAATTGATGTATCA',
  '497' => 'GAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTG',
  '213' => 'GAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGT',
  '472' => 'TATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGT',
  '456' => 'TGAGGCCCAATTGATGTATCAACATGTGGTTTTTATAAAAA',
  '439' => 'ATGGGCCGGTTGATCACTGAGGCCCAATTGATGTATCAACA',
  '362' => 'GACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCT',
  '265' => 'AGACTTTCTATTCCTTATTTATACCCGTCCGTATAGTAGGG',
  '296' => 'TATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTT',
  '317' => 'AATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACA',
  '411' => 'GGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGA',
  '493' => 'AAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTG',
  '478' => 'CATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTC',
  '384' => 'GTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGT',
  '398' => 'GTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAA',
  '386' => 'GGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCG',
  '407' => 'TGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGG',
  '445' => 'CGGTTGATCACTGAGGCCCAATTGATGTATCAACATGTGGT',
  '200' => 'GACTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAG',
  '366' => 'AGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAG',
  '376' => 'TCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGG',
  '329' => 'TGCCCCTTTTCTACGCAGCCGAACTCACATCCTGACCAGGC',
  '272' => 'CTATTCCTTATTTATACCCGTCCGTATAGTAGGGTACGCAA',
  '298' => 'TAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTT',
  '400' => 'TCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAAT',
  '313' => 'GGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACT',
  '231' => 'CTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTT',
  '243' => 'CTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATT',
  '504' => 'GAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCA',
  '343' => 'GCAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTG',
  '468' => 'GATGTATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAA',
  '287' => 'ACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGA',
  '475' => 'CAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTT',
  '441' => 'GGGCCGGTTGATCACTGAGGCCCAATTGATGTATCAACATG',
  '294' => 'CGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCC',
  '397' => 'AGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGA',
  '413' => 'TCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATC',
  '349' => 'GAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGC',
  '275' => 'TTCCTTATTTATACCCGTCCGTATAGTAGGGTACGCAAGGT',
  '197' => 'GCAGACTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTC',
  '203' => 'TTCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTC',
  '261' => 'CCCGAGACTTTCTATTCCTTATTTATACCCGTCCGTATAGT',
  '459' => 'GGCCCAATTGATGTATCAACATGTGGTTTTTATAAAAAGAG',
  '321' => 'CTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACATCCT',
  '432' => 'CCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTGATGT',
  '284' => 'TATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCG',
  '247' => 'AAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTATA',
  '371' => 'GGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATG',
  '204' => 'TCTCCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCT',
  '289' => 'CCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGT',
  '346' => 'GCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTG',
  '435' => 'TGAAATGGGCCGGTTGATCACTGAGGCCCAATTGATGTATC',
  '401' => 'CTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATG',
  '427' => 'TTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATT',
  '333' => 'CCTTTTCTACGCAGCCGAACTCACATCCTGACCAGGCCGGG',
  '339' => 'CTACGCAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGG',
  '228' => 'ATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGAC',
  '323' => 'CGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACATCCTGA',
  '268' => 'CTTTCTATTCCTTATTTATACCCGTCCGTATAGTAGGGTAC',
  '345' => 'AGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGGT',
  '319' => 'TTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACATC',
  '224' => 'TGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCG',
  '223' => 'GTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCC',
  '404' => 'AAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGC',
  '446' => 'GGTTGATCACTGAGGCCCAATTGATGTATCAACATGTGGTT',
  '282' => 'TTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCT',
  '262' => 'CCGAGACTTTCTATTCCTTATTTATACCCGTCCGTATAGTA',
  '420' => 'TGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGG',
  '212' => 'CGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAG',
  '344' => 'CAGCCGAACTCACATCCTGACCAGGCCGGGCTTCGGCCTGG',
  '352' => 'CTCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGG',
  '494' => 'AAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGT',
  '487' => 'TTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCC',
  '238' => 'CAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCC',
  '251' => 'AGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTATACCCG',
  '426' => 'GTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAAT',
  '369' => 'CCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGA',
  '253' => 'TCCCCTTTCCCGAGACTTTCTATTCCTTATTTATACCCGTC',
  '279' => 'TTATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAAT',
  '448' => 'TTGATCACTGAGGCCCAATTGATGTATCAACATGTGGTTTT',
  '209' => 'CGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAA',
  '498' => 'AGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGG',
  '483' => 'GGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAA',
  '256' => 'CCTTTCCCGAGACTTTCTATTCCTTATTTATACCCGTCCGT',
  '216' => 'CTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCC',
  '357' => 'ATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGA',
  '372' => 'GGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGG',
  '428' => 'TATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAATTG',
  '455' => 'CTGAGGCCCAATTGATGTATCAACATGTGGTTTTTATAAAA',
  '227' => 'TATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGA',
  '336' => 'TTTCTACGCAGCCGAACTCACATCCTGACCAGGCCGGGCTT',
  '457' => 'GAGGCCCAATTGATGTATCAACATGTGGTTTTTATAAAAAG',
  '383' => 'GGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGG',
  '500' => 'TCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTA',
  '264' => 'GAGACTTTCTATTCCTTATTTATACCCGTCCGTATAGTAGG',
  '255' => 'CCCTTTCCCGAGACTTTCTATTCCTTATTTATACCCGTCCG',
  '297' => 'ATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTT',
  '359' => 'CCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGT',
  '277' => 'CCTTATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGA',
  '462' => 'CCAATTGATGTATCAACATGTGGTTTTTATAAAAAGAGTCG',
  '232' => 'TTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGAGACTTTC',
  '414' => 'CGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCA',
  '477' => 'ACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCT',
  '225' => 'GCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTTTCCCGA',
  '207' => 'CCCGACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGA',
  '330' => 'GCCCCTTTTCTACGCAGCCGAACTCACATCCTGACCAGGCC',
  '263' => 'CGAGACTTTCTATTCCTTATTTATACCCGTCCGTATAGTAG',
  '394' => 'TCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCT',
  '505' => 'AGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAA',
  '360' => 'CTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTT',
  '419' => 'CTGGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAG',
  '513' => 'TTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACTTCATT',
  '290' => 'CGTCCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTG',
  '304' => 'GGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGC',
  '476' => 'AACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTC',
  '210' => 'GACGAGCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAA',
  '406' => 'GTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCG',
  '258' => 'TTTCCCGAGACTTTCTATTCCTTATTTATACCCGTCCGTAT',
  '396' => 'GAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTG',
  '482' => 'TGGTTTTTATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAA',
  '510' => 'AGTTTTCTCTAAAAATCCCTTGTGTTTGGTAATCAAACTTC',
  '499' => 'GTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGT',
  '393' => 'CTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTATTCC',
  '449' => 'TGATCACTGAGGCCCAATTGATGTATCAACATGTGGTTTTT',
  '293' => 'CCGTATAGTAGGGTACGCAAGGTGAATTCTCGAGAGTGCCC',
  '274' => 'ATTCCTTATTTATACCCGTCCGTATAGTAGGGTACGCAAGG',
  '365' => 'CAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCGAGTTCTAAA',
  '306' => 'TACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAG',
  '470' => 'TGTATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGA',
  '322' => 'TCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACATCCTG',
  '469' => 'ATGTATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAG',
  '353' => 'TCACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGC',
  '375' => 'TTCGGCCTGGTGGGCCGGCTCGAGTTCTAAAGTGATGGTCG',
  '403' => 'AAAGTGATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGG',
  '466' => 'TTGATGTATCAACATGTGGTTTTTATAAAAAGAGTCGTGAG',
  '252' => 'GTCCCCTTTCCCGAGACTTTCTATTCCTTATTTATACCCGT',
  '310' => 'CAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGA',
  '283' => 'TTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAATTCTC',
  '303' => 'GGGTACGCAAGGTGAATTCTCGAGAGTGCCCCTTTTCTACG',
  '250' => 'AAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTATACCC',
  '421' => 'GGGTCGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGC',
  '501' => 'CGTGAGAAGAGTTTTCTCTAAAAATCCCTTGTGTTTGGTAA',
  '215' => 'GCTGACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCC',
  '450' => 'GATCACTGAGGCCCAATTGATGTATCAACATGTGGTTTTTA',
  '278' => 'CTTATTTATACCCGTCCGTATAGTAGGGTACGCAAGGTGAA',
  '490' => 'ATAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTT',
  '271' => 'TCTATTCCTTATTTATACCCGTCCGTATAGTAGGGTACGCA',
  '387' => 'GGCCGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGT',
  '245' => 'CGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTA',
  '267' => 'ACTTTCTATTCCTTATTTATACCCGTCCGTATAGTAGGGTA',
  '491' => 'TAAAAAGAGTCGTGAGAAGAGTTTTCTCTAAAAATCCCTTG',
  '354' => 'CACATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCT',
  '461' => 'CCCAATTGATGTATCAACATGTGGTTTTTATAAAAAGAGTC',
  '219' => 'ACTTGTGCTATCCTTTCTCCAGGTCTCGAAAAAGTCCCCTT',
  '318' => 'ATTCTCGAGAGTGCCCCTTTTCTACGCAGCCGAACTCACAT',
  '257' => 'CTTTCCCGAGACTTTCTATTCCTTATTTATACCCGTCCGTA',
  '473' => 'ATCAACATGTGGTTTTTATAAAAAGAGTCGTGAGAAGAGTT',
  '248' => 'AAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTATTTATAC',
  '332' => 'CCCTTTTCTACGCAGCCGAACTCACATCCTGACCAGGCCGG',
  '390' => 'CGGCTCGAGTTCTAAAGTGATGGTCGGGGCTGGGTCGTTAT',
  '425' => 'CGTTATTCCTTGAAATGGGCCGGTTGATCACTGAGGCCCAA',
  '356' => 'CATCCTGACCAGGCCGGGCTTCGGCCTGGTGGGCCGGCTCG',
  '408' => 'GATGGTCGGGGCTGGGTCGTTATTCCTTGAAATGGGCCGGT',
  '196' => 'AGCAGACTTCTCCCGACGAGCTGACTTGTGCTATCCTTTCT',
  '242' => 'TCTCGAAAAAGTCCCCTTTCCCGAGACTTTCTATTCCTTAT'
};

my ($left, $right, $id) = (0, 40, 1);
my @plaps;
my $len = length($seq);
while ($left < ($len - 40))
{
  my $ol = substr($seq, $left, $right - $left + 1);
  my $lap = Bio::Seq->new(
      -seq => $ol,
      -id  => $id
  );
  push @plaps, $lap;
  $id++;
  $left++;
  $right = $left + 40;
}
my $arr = $DNA->filter_vmatch(
    -sequences => \@plaps,
    -parent => $seqobj,
    -mismatches => 10);
my %thsh = map {$_->id => $_->seq} @$arr;

is_deeply (\%thsh, $rhsh, "filter vmatch");
                                         

