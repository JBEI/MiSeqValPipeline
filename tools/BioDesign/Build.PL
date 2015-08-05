use 5.006;

use strict;
use warnings;

use Module::Build;
use File::Basename;
use Carp;

print "\n";

my $llt = 0;
my ($cpath,  $spath,  $tpath) =  (       q{},                 q{},        q{});
my ($dcpath, $dspath, $dtpath) = ('/etc/GeneDesign/', '/usr/local/bin', '/tmp');
my ($g, $dg) =  ( 0,   q{Y});

my $check = eval
{
  require Bio::GeneDesign::ConfigData;
};
if ($check)
{
  $dcpath = Bio::GeneDesign::ConfigData->config('conf_path')        || $dcpath;
  $dtpath = Bio::GeneDesign::ConfigData->config('tmp_path')         || $dtpath;
  $dspath = Bio::GeneDesign::ConfigData->config('script_path')      || $dspath;
  $dg =     Bio::GeneDesign::ConfigData->config('graphing_support') || $dg;
}

my $GDB = Module::Build->new
(
    module_name         => 'Bio::GeneDesign',
    license             => 'bsd',
    dist_author         => q{Sarah Richardson <SMRichardson@lbl.gov>},
    dist_version        => '5.52',
    dist_abstract       => 'Functions for the design of synthetic genes',
    add_to_cleanup      => [ 'Bio::GeneDesign-*' ],
    create_makefile_pl  => 'traditional',
    build_requires      =>
    {
        'Test::More'      => 0.8,
        'Test::Deep'      => 0.1,
    },
    requires            =>
    {
        'perl'                => 5.006,
        'List::Util'          => 1.22,
        'Math::Combinatorics' => 0.09,
        'Text::Wrap'          => 2009.0305,
        'Exception::Class'    => 0,
        'Bio::Root::Root'     => 0
    },
    configure_requires  =>
    {
      'Module::Build'     => 0.38
    },
    auto_features       =>
    {
      graphing =>
      {
        description => 'Use GD::Graphics to draw graphs',
        requires    =>
        {
          'GD::Graph::lines'  => 1.15,
          'GD::Graph::colour' => 1.10,
          'GD::Image'         => 0
        }
      },
    },
    script_files =>
    [
      'bin/GD_Filter_Enzymes.pl',
      'bin/GD_Generate_RSCU_Table.pl',
      'bin/GD_Graph_Dotplot.pl',
      'bin/GD_Juggle_Codons.pl',
      'bin/GD_Reverse_Translate.pl',
      'bin/GD_Sequence_Subtraction.pl',
      'bin/GD_Graph_RSCU_Values.pl'
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
      graphing_support =>
      {
        store => \$g,
        type => '=s',
      },
    }
);

if ($llt == 0)
{
  if (! $cpath)
  {
    $cpath = $GDB->prompt('Where should configuration files be installed?', $dcpath);
  }
  if (! $spath)
  {
    $spath = $GDB->prompt('Where should scripts be installed?', $dspath);
  }
  if (! $tpath)
  {
    $tpath = $GDB->prompt('Where should GeneDesign write tmp files?', $dtpath);
  }
  if (! $g && $GDB->feature('graphing'))
  {
    $g = $GDB->y_n('Enable GD::Graphics support?', $dg);
  }
}
else
{
  $cpath = $cpath || $dcpath;
  $tpath = $tpath || $dtpath;
  $spath = $spath || $dspath;
  $g = $g || $dg;
}

$GDB->config_data(conf_path => $cpath);
$GDB->config_data(tmp_path => $tpath);
$GDB->config_data(script_path => $spath);
$GDB->config_data(graphing_support => $g) if ($GDB->feature('graphing'));

#Prepare configuration directory
my $tcp = $GDB->config_data('conf_path');
$GDB->config_data(conf_path => $tcp);
my $confs =
[
  'codon_tables/arabidopsis_thaliana.rscu',
  'codon_tables/Bacillus_subtilis.rscu',
  'codon_tables/Caenorhabditis_elegans.rscu',
  'codon_tables/Cglut.rscu',
  'codon_tables/Deinococcus_radiodurans.rscu',
  'codon_tables/Drosophila_melanogaster.rscu',
  'codon_tables/Escherichia_coli.rscu',
  'codon_tables/Flat.rscu',
  'codon_tables/Homo_sapiens.rscu',
  'codon_tables/Mycoplasma_genitalium.ct',
  'codon_tables/Mycoplasma_genitalium.rscu',
  'codon_tables/oryza_sativa.rscu',
  'codon_tables/Saccharomyces_cerevisiae.rscu',
  'codon_tables/Standard.ct',
  'codon_tables/yeast.rscu',
  'enzymes/all_enzymes',
  'enzymes/blunts',
  'enzymes/IIB',
  'enzymes/nonpal',
  'enzymes/nonpal_short',
  'enzymes/nonpal_and_IIB',
  'enzymes/outside',
  'enzymes/standard',
  'enzymes/standard_and_IIB',
  'enzymes/test',
];
$GDB->add_build_element('GeneDesign');
process_conf_files($GDB, $confs);
$GDB->install_path(GeneDesign => $GDB->config_data('conf_path'));
print 'Configuration files will be in ', $GDB->config_data('conf_path'), "\n";

$GDB->install_path(script => $GDB->config_data('script_path'));
print 'Scripts will be in ', $GDB->config_data('script_path'), "\n";

print 'Temporary files will be written to ', $GDB->config_data('tmp_path'), "\n";

print "\n";

$GDB->create_build_script();

print "\n";

sub process_conf_files
{
  my $build = shift;
  my $files = shift;
  return unless $files;

  my $conf_dir = File::Spec->catdir($build->blib, 'GeneDesign');
  File::Path::mkpath( $conf_dir );

  foreach my $file (@{$files})
  {
    my $result = $build->copy_if_modified($file, $conf_dir) or next;
    $build->fix_shebang_line($result) unless $build->is_vmsish;
  }
  return;
}
