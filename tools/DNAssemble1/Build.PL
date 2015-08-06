use 5.006;

use strict;
use warnings;

use Module::Build;
use File::Basename;
use Carp;

print "\n";

my $llt = 0;
my ($cpath,  $spath,  $tpath) =  (       q{},                 q{},        q{});
my ($dcpath, $dspath, $dtpath) = ('/etc/DNAssemble/', '/usr/local/bin', '/tmp');
my ($e,  $bl,  $v) =  ( 0,   0,   0);
my ($de, $dbl, $dv) = ('y', 'y', 'y');

my $check = eval
{
  require Bio::DNAssemble::ConfigData;
};
if ($check)
{
  $dcpath = Bio::DNAssemble::ConfigData->config('conf_path')        || $dcpath;
  $dtpath = Bio::DNAssemble::ConfigData->config('tmp_path')         || $dtpath;
  $dspath = Bio::DNAssemble::ConfigData->config('script_path')      || $dspath;
  $dbl = Bio::DNAssemble::ConfigData->config('BLAST_support') ? 'y' : 'n';
  $de =  Bio::DNAssemble::ConfigData->config('EMBOSS_support') ? 'y' : 'n';
  $dv =  Bio::DNAssemble::ConfigData->config('vmatch_support') ? 'y' : 'n';
}

my $DNAB = Module::Build->new
(
    module_name         => 'Bio::DNAssemble',
    license             => 'bsd',
    dist_author         => q{Sarah Richardson <SMRichardson@lbl.gov>},
    dist_version        => '1.00',
    dist_abstract       => 'Functions for the fragmentation of synthetic DNA',
    add_to_cleanup      => [ 'Bio::DNAssemble-*' ],
    create_makefile_pl  => 'traditional',
    build_requires      =>
    {
        'Test::More'      => 0.8,
        'Test::Deep'      => 0.1,
    },
    requires            =>
    {
        'perl'                => 5.006,
        'Bio::Root::Root'     => 0
    },
    configure_requires  =>
    {
      'Module::Build'     => 0.38
    },
    auto_features       =>
    {
      palindrome =>
      {
        description => 'Enable EMBOSS palindrome for hairpin detection',
        requires    =>
        {
          'Bio::Factory::EMBOSS' => 0,
          'Bio::Tools::EMBOSS::Palindrome' => 0
        }
      },
      blast =>
      {
        description => 'Enable BLAST+ for similarity detection',
        requires    =>
        {
          'Bio::Tools::Run::StandAloneBlastPlus' => 0
        }
      },
      vmatch =>
      {
        description => 'Enable vmatch for similarity detection',
        requires    =>
        {
          'Bio::Grep::Backend::Vmatch' => '0.10.6'
        }
      }
    },
    script_files =>
    [
      'bin/JGISB_96_to_384.pl',
      'bin/JGISB_Analyze_PacBio.pl',
      'bin/JGISB_Design_Fixing_Primers.pl',
      'bin/JGISB_CompareConstructs.pl',
      'bin/JGISB_Primer_forEnzfromGB.pl',
      'bin/DNAssemble_001_Classic_SGE.pl',
      'bin/DNAssemble_002_Classic_Local.pl',
      'bin/DNAssemble_003_GBlock_SGE.pl',
      'bin/DNAssemble_004_ScreenedGBlock_SGE.pl',
      'bin/DNAssemble_100_CreateXML.pl',
      'bin/DNAssemble_110_Cognominate.pl',
      'bin/DNAssemble_115_Vector_Flank.pl',
      'bin/DNAssemble_117_GBlockable.pl',
      'bin/DNAssemble_120_Carve.pl',
      'bin/DNAssemble_121_Carve_Simple.pl',
      'bin/DNAssemble_122_Carve_Simple.pl',
      'bin/DNAssemble_130_Compress.pl',
      'bin/DNAssemble_140_Cluster.pl',
      'bin/DNAssemble_150_InspectGBlocks.pl',
      'bin/DNAssemble_155_CorrectGBlocks.pl',
      'bin/DNAssemble_160_Chop.pl',
      'bin/DNAssemble_170_Compress2.pl',
      'bin/DNAssemble_200_Check.pl',
      'bin/DNAssemble_201_Check_GBlock.pl',
      'bin/DNAssemble_300_Summary.pl',
      'bin/DNAssemble_310_Extract.pl',
      'bin/DNAssemble_320_Stitchers.pl',
      'bin/DNAssemble_Classic.pl',
      'bin/DNAssemble_GBlock.pl',
      'bin/DNAssemble_Local.pl',
      'bin/DNAssemble_ScreenedGBlock.pl',
    ],
    get_options =>
    {
      like_last_time =>
      {
        store => \$llt,
        default => 0,
      },
      conf_path =>
      {
        store => \$cpath,
        type => '=s',
      },
      script_path =>
      {
        store => \$spath,
        type => '=s',
      },
      tmp_path =>
      {
        store => \$tpath,
        type => '=s',
      },
      EMBOSS_support =>
      {
        store => \$e,
        type => '=s',
      },
      BLAST_support =>
      {
        store => \$bl,
        type => '=s',
      },
      vmatch_support =>
      {
        store => \$v,
        type => '=s',
      },
    }
);

if ($llt == 0)
{
  if (! $cpath)
  {
    $cpath = $DNAB->prompt('Where should configuration files be installed?', $dcpath);
  }
  if (! $spath)
  {
    $spath = $DNAB->prompt('Where should scripts be installed?', $dspath);
  }
  if (! $tpath)
  {
    $tpath = $DNAB->prompt('Where should DNAssemble write tmp files?', $dtpath);
  }
  if (! $e && $DNAB->feature('palindrome'))
  {
    $e = $DNAB->y_n('Enable EMBOSS palindrome for hairpin detection?', $de);
  }
  if (! $bl && $DNAB->feature('blast'))
  {
    $bl = $DNAB->y_n('Enable BLAST+ for similarity detection?', $dbl);
  }
  if (! $v && $DNAB->feature('vmatch'))
  {
    $v = $DNAB->y_n('Enable vmatch for similarity detection?', $dv);
  }
}
else
{
  $cpath = $cpath || $dcpath;
  $tpath = $tpath || $dtpath;
  $spath = $spath || $dspath;
  $e =     defined $e  ? $e  : $de;
  $bl =    defined $bl ? $bl : $dbl;
  $v =     defined $v  ? $v  : $dv;
}

$DNAB->config_data(conf_path => $cpath);
$DNAB->config_data(tmp_path => $tpath);
$DNAB->config_data(script_path => $spath);
$DNAB->config_data(EMBOSS_support => $e);
$DNAB->config_data(BLAST_support => $bl);
$DNAB->config_data(vmatch_support => $v);

#Prepare configuration directory
my $tcp = $DNAB->config_data('conf_path');
$DNAB->config_data(conf_path => $tcp);
my $confs =
[
  'vectors/pENTR_SD_SmaI.genbank',
  'vectors/pEC-K18mob2.genbank',
  'vectors/pET29b-plus_SMR.genbank',
  'vectors/pET45b-SmaI.genbank',
  'vectors/peU_HSBC_NcoI_RAH_GD.genbank',
  'vectors/pFIL-B-2_ScaI.genbank',
  'vectors/pK18mobsacB.genbank',
  'vectors/pMCC.genbank',
  'vectors/pMCL200.genbank',
  'vectors/pMCL200-T7_LacI.genbank',
  'vectors/pMCL200-T7_LacI_mod.genbank',
  'vectors/pMCL200Y.genbank',
  'vectors/pMCL200Y_wA_pyrG_alcA.genbank',
  'vectors/pMCL200Y_yA_pyroA_alcA.genbank',
  'vectors/pNJ020_GD.genbank',
  'vectors/pNJ022_GD.genbank',
  'vectors/pRS426HygMX.genbank',
  'vectors/pVA1.genbank',
  'vectors/pYJKSD.genbank',
  'vectors/YEp-Leu-Flag_LtOMT.genbank',
  'vectors/YEpADH2-URA-FLAGLink.genbank',
  'vectors/YEpADH2-Trp-FLAGLink.genbank'
];
$DNAB->add_build_element('DNAssemble');
process_conf_files($DNAB, $confs);
$DNAB->install_path(DNAssemble => $DNAB->config_data('conf_path'));
print 'Configuration files will be in ', $DNAB->config_data('conf_path'), "\n";

$DNAB->install_path(script => $DNAB->config_data('script_path'));
print 'Scripts will be in ', $DNAB->config_data('script_path'), "\n";

print 'Temporary files will be written to ', $DNAB->config_data('tmp_path'), "\n";

if ($DNAB->config_data('BLAST_support'))
{
  my $blq = `which blastn`;
  croak ('Either BLAST+ is not installed or it is not on my PATH') unless ($blq);
  my ($y, $blpath) = fileparse($blq);
  $DNAB->config_data('blast_path' => $blpath);
  print "Will use BLAST+ executables found in $blpath\n";
}
if ($DNAB->config_data('vmatch_support'))
{
  my $vl = `which vmatch`;
  croak ('Either vmatch is not installed or it is not on my PATH') unless ($vl);
  my ($y, $vlpath) = fileparse($vl);
  print "Will use vmatch executables found in $vlpath\n";
}

print "\n";

$DNAB->create_build_script();

print "\n";

sub process_conf_files
{
  my $build = shift;
  my $files = shift;
  return unless $files;

  my $conf_dir = File::Spec->catdir($build->blib, 'DNAssemble');
  File::Path::mkpath( $conf_dir );

  foreach my $file (@{$files})
  {
    my $result = $build->copy_if_modified($file, $conf_dir) or next;
    $build->fix_shebang_line($result) unless $build->is_vmsish;
  }
  return;
}