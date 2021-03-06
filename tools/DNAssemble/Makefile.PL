# Note: this file was auto-generated by Module::Build::Compat version 0.4007
require 5.006;
use ExtUtils::MakeMaker;
WriteMakefile
(
  'PL_FILES' => {},
  'INSTALLDIRS' => 'site',
  'NAME' => 'Bio::DNAssemble',
  'EXE_FILES' => [
                   'bin/DNAssemble_001_Classic.pl',
                   'bin/DNAssemble_100_Clean.pl',
                   'bin/DNAssemble_120_Carve.pl',
                   'bin/DNAssemble_130_Compress1.pl',
                   'bin/DNAssemble_140_Cluster.pl',
                   'bin/DNAssemble_160_Chop.pl',
                   'bin/DNAssemble_170_Compress2.pl',
                   'bin/DNAssemble_200_Check.pl'
                 ],
  'VERSION_FROM' => 'lib/Bio/DNAssemble.pm',
  'PREREQ_PM' => {
                   'Test::More' => '0.8',
                   'Bio::Root::Root' => 0,
                   'Test::Deep' => '0.1'
                 }
)
;
