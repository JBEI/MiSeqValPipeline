#!/usr/bin/env perl

# SGE commands
#
#$ -l h_rt=14400
#$ -l ram.c=5G
#$ -pe pe_slots 1
#$ -cwd
#$ -V
#$ -P gentech-rnd.p

use EnvironmentModules;
use Bio::DNAssemble;
use Getopt::Long;
use Pod::Usage;
use Readonly;
use English qw( -no_match_vars );
use autodie qw(open close);
use Carp;

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_004_ScreenedGBlock_SGE_' . $VERSION;
my $DNAS = '_001';

module('load dnassemble');
module('load genedesign');
module('load openmpi');
module('load taskfarmermq/2.4');

local $OUTPUT_AUTOFLUSH = 1;

# Get command line arguments
#
my %p = ();
GetOptions (
  'help'                 => \$p{HELP},
  'input=s'              => \$p{INPUT},
  'assemblyTm:i'         => \$p{ASSEMBLYTM},
  'slots:i'              => \$p{SLOTS},
  'name:s'               => \$p{NAME},
                         
  'anonkey=s'            => \$p{ANONKEY},
  'procdir=s'            => \$p{PROCDIR},
  'logfile=s'            => \$p{LOGFILE},
  
  'SKIP_WARNINGS'        => \$p{SKIP_WARNINGS},
  'DESTINATION_VECTOR'   => \$p{DESTINATION_VECTOR},
  'GC_GLOBAL_MIN:i'      => \$p{GC_GLOBAL_MIN},
  'GC_GLOBAL_MAX:i'      => \$p{GC_GLOBAL_MAX},
  'GC_REGIONAL_WINDOW'   => \$p{GC_REGIONAL_WINDOW},
  'GC_REGIONAL_MIN:i'    => \$p{GC_REGIONAL_MIN},
  'GC_REGIONAL_MAX:i'    => \$p{GC_REGIONAL_MAX},
  'GC_LOCAL_WINDOW'      => \$p{GC_LOCAL_WINDOW},
  'GC_LOCAL_MIN:i'       => \$p{GC_LOCAL_MIN},
  'GC_LOCAL_MAX:i'       => \$p{GC_LOCAL_MAX},
  'GC_TERMINAL_WINDOW'   => \$p{GC_TERMINAL_WINDOW},
  'GC_TERMINAL_MIN:i'    => \$p{GC_TERMINAL_MIN},
  'GC_TERMINAL_MAX:i'    => \$p{GC_TERMINAL_MAX},
  'HP_A_LENGTH:i'        => \$p{HP_A_LENGTH},
  'HP_T_LENGTH:i'        => \$p{HP_T_LENGTH},
  'HP_C_LENGTH:i'        => \$p{HP_C_LENGTH},
  'HP_G_LENGTH:i'        => \$p{HP_G_LENGTH},
  'DIMER_LIMIT:i'        => \$p{DIMER_LIMIT},
  'TRIMER_LIMIT:i'       => \$p{TRIMER_LIMIT},
  'DR_LOCAL_WINDOW:i'    => \$p{DR_LOCAL_WINDOW},
  'DR_LOCAL_RATIO:i'     => \$p{DR_LOCAL_RATIO},
  'DR_REGIONAL_WINDOW:i' => \$p{DR_REGIONAL_WINDOW},
  'DR_REGIONAL_RATIO:i'  => \$p{DR_REGIONAL_RATIO},
  'DR_GLOBAL_RATIO:i'    => \$p{DR_GLOBAL_RATIO},
  'DR_SOLO_RATIO:i'      => \$p{DR_SOLO_RATIO},
  'DR_LENGTH:i'          => \$p{DR_LENGTH},
  'IR_LENGTH_MIN:i'      => \$p{IR_LENGTH_MIN},
  'IR_LENGTH_MAX:i'      => \$p{IR_LENGTH_MAX},
  'IR_WIDTH:i'           => \$p{IR_WIDTH},

  'GBLOCK_LENGTH_MIN:i'  => \$p{GBLOCK_LENGTH_MIN},
  'GBLOCK_LENGTH_MAX:i'  => \$p{GBLOCK_LENGTH_MAX},
  'GBLOCK_LENGTH:i'      => \$p{GBLOCK_LENGTH},
  'GBLOCK_OVERLAP:i'     => \$p{GBLOCK_OVERLAP},
  'GBLOCK_OVERLAP_MIN:i' => \$p{GBLOCK_OVERLAP_MIN},

  'OLIGO_LENGTH_MIN:i'   => \$p{OLIGO_LENGTH_MIN},
  'OLIGO_LENGTH_MAX:i'   => \$p{OLIGO_LENGTH_MAX},
  'OLIGO_LENGTH:i'       => \$p{OLIGO_LENGTH},
  'OLIGO_OVERLAP_MIN:i'  => \$p{OLIGO_OVERLAP_MIN},
);
if ($p{HELP})
{
  pod2usage(-verbose=>99);
}


################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $DNA = Bio::DNAssemble->new();

# Non user configurable variables - these are defaults for building block and
# oligo size, etc. I also define the scripts to be used at the different stages
# here. I don't like hard coding $BIN. I will talk to Kirsten and Doug about
# that.
#

Readonly my $BIN        => '/usr/common/jgi/synthbio/DNAssemble/bin/';
Readonly my $CLEAN      => $BIN . 'DNAssemble_100_CreateXML.pl';
Readonly my $COGNOM     => $BIN . 'DNAssemble_110_Cognominate.pl';
Readonly my $VECFLANK   => $BIN . 'DNAssemble_115_Vector_Flank.pl';
Readonly my $GBLOCKING  => $BIN . 'DNAssemble_117_GBlockable.pl';
Readonly my $CARVE      => $BIN . 'DNAssemble_122_Carve_Simple.pl';
Readonly my $INSPECT    => $BIN . 'DNAssemble_150_InspectGBlocks.pl';
Readonly my $COMPRESSA  => $BIN . 'DNAssemble_130_Compress1.pl';
Readonly my $CHECK      => $BIN . 'DNAssemble_201_Check_GBlock.pl';
Readonly my $SUMMARY    => $BIN . 'DNAssemble_300_Summary.pl';
Readonly my $EXTRACT    => $BIN . 'DNAssemble_310_Extract.pl';

my $SGEFLAGS = "-l ram.c=5.25G -pe pe_fill $p{SLOTS} -cwd -V";
my $subregex = qr{Your [ ] job [ ] (\d+) [ ] [(]["]}msx;

# The input file must exist. We will not do deep parsing but it must pass an
# initial GeneDesign load to prove formatting.
#
croak "\nDNA_ERROR: You must supply an input file.\n\n" if (! $p{INPUT});
croak "\nDNA_ERROR: $p{INPUT} does not exist.\n\n" if (! -e $p{INPUT});
my ($iterator, $filename, $filesuffix) = $DNA->import_seqs($p{INPUT});

# We need a processing directory and a logfile.
#
croak "\nDNA_ERROR: No processing directory named!\n\n" if (! $p{PROCDIR});
croak "\nDNA_ERROR: No processing directory exists!\n\n" if (!-e $p{PROCDIR});
croak "\nDNA_ERROR: No logfile named!\n\n" if (! $p{LOGFILE});
$p{PROCDIR} = $DNA->endslash($p{PROCDIR});
my $FILEPREF = $p{PROCDIR} . $p{ANONKEY};
_writelog("Processing in $p{PROCDIR}");

# Every path in this array will be unlinked at the end of execution.
#
my @CLEANUP = ();

################################################################################
#################################### XMLING ####################################
################################################################################
my $LASTFILE = $p{INPUT};
my $NEXTFILE = $FILEPREF . q{_100.xml};
my %cleanopts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  name     => $p{NAME},
  logfile  => $p{LOGFILE},
);
_runscript($CLEAN, \%cleanopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
################################ COGNOMINATING  ################################
################################################################################
$NEXTFILE = $FILEPREF . q{_110.xml};
my %cogopts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  logfile  => $p{LOGFILE},
  anonkey  => $p{ANONKEY},
);
_runscript($COGNOM, \%cogopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
############################### VECTOR FLANKING  ###############################
################################################################################
$NEXTFILE = $FILEPREF . q{_115.xml};
my %vecopts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  logfile  => $p{LOGFILE},
  vector   => $p{DESTINATION_VECTOR},
  overlap  => $p{GBLOCK_OVERLAP},
);
_runscript($VECFLANK, \%vecopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
################################## GBLOCKING  ##################################
################################################################################
$NEXTFILE = $FILEPREF . q{_117.xml};
my %gblopts = (
  input             => $LASTFILE,
  output            => $NEXTFILE,
  logfile           => $p{LOGFILE},
  gblock_length_min => $p{GBLOCK_LENGTH_MIN},
);
_runscript($GBLOCKING, \%gblopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
################################## PRECARVING ##################################
################################################################################
my %precarveopts = (
  GBLOCK_LENGTH      => $p{GBLOCK_LENGTH},
  GBLOCK_LENGTH_MAX  => $p{GBLOCK_LENGTH_MAX},
  GBLOCK_LENGTH_MIN  => $p{GBLOCK_LENGTH_MIN},
  GBLOCK_OVERLAP     => $p{GBLOCK_OVERLAP},
  GBLOCK_OVERLAP_MIN => $p{GBLOCK_OVERLAP_MIN},
  GC_TERMINAL_MIN    => $p{GC_TERMINAL_MIN},
  GC_TERMINAL_MAX    => $p{GC_TERMINAL_MAX},
  GC_TERMINAL_WINDOW => $p{GC_TERMINAL_WINDOW},
  ASSEMBLYTM         => $p{ASSEMBLYTM},
  DR_LENGTH          => $p{DR_LENGTH},
);
my $CMD_CARVE = _assemblecmd($CARVE, \%carveopts);

# Take each chunk that was split up in the previous step and create a taskfarmer
# work order for it.
#
my $carvedesign = $DNA->load_design_from_file($LASTFILE);
my @chlist;
my @carvelogs;
my $carvetasklist = q{};
my @constructs = $carvedesign->get_constructs();
foreach my $construct ( sort {$b->seqlen <=> $a->seqlen} @constructs )
{
  my $prefix = $p{PROCDIR} . $construct->id();
  my $inp = $prefix . q{.xml};
  my $log = $prefix . q{_LOG.txt};
  my $out = $prefix . q{_122.xml};
  my $cmd = $CMD_CARVE . q{ -input } . $inp . q{ -logfile } . $log;
  $cmd .= q{ -output } . $out;
  $carvetasklist .= $cmd . q{:} . $log . q{,} . $out . q{:0} . "\n";
  _constructdump($construct, $inp);
  push @chlist, $out;
  push @CLEANUP, $inp;
  push @carvelogs, $log;
}
_taskfarm($carvetasklist, 'CARVE');

$NEXTFILE = $FILEPREF . q{_122.xml};
_combinefiles($LASTFILE, \@chlist, $NEXTFILE, 0);
push @CLEANUP, @chlist;
$LASTFILE = $NEXTFILE;

my $templog = $p{PROCDIR} . 'templog';
my $carvecatcmd = "cat $p{LOGFILE} @carvelogs > $templog";
$DNA->safeopen($carvecatcmd);
my $carvecpcmd = "cp $templog $p{LOGFILE}";
$DNA->safeopen($carvecpcmd);
push @CLEANUP, @carvelogs, $templog;


################################################################################
################################## INSPECTING ##################################
################################################################################
$NEXTFILE = $FILEPREF . q{_150.xml};
my %inspectopts = (
  input              => $LASTFILE,
  output             => $NEXTFILE,
  logfile            => $p{LOGFILE},
  SKIP_WARNINGS      => $p{SKIP_WARNINGS},
  GC_GLOBAL_MIN      => $p{GC_GLOBAL_MIN},
  GC_GLOBAL_MAX      => $p{GC_GLOBAL_MAX},
  GC_REGIONAL_WINDOW => $p{GC_REGIONAL_WINDOW},
  GC_REGIONAL_MIN    => $p{GC_REGIONAL_MIN},
  GC_REGIONAL_MAX    => $p{GC_REGIONAL_MAX},
  GC_LOCAL_WINDOW    => $p{GC_LOCAL_WINDOW},
  GC_LOCAL_MIN       => $p{GC_LOCAL_MIN},
  GC_LOCAL_MAX       => $p{GC_LOCAL_MAX},
  GC_TERMINAL_WINDOW => $p{GC_TERMINAL_WINDOW},
  GC_TERMINAL_MIN    => $p{GC_TERMINAL_MIN},
  GC_TERMINAL_MAX    => $p{GC_TERMINAL_MAX},
  HP_A_LENGTH        => $p{HP_A_LENGTH},
  HP_T_LENGTH        => $p{HP_T_LENGTH},
  HP_C_LENGTH        => $p{HP_C_LENGTH},
  HP_G_LENGTH        => $p{HP_G_LENGTH},
  DIMER_LIMIT        => $p{DIMER_LIMIT},
  TRIMER_LIMIT       => $p{TRIMER_LIMIT},
  DR_LOCAL_WINDOW    => $p{DR_LOCAL_WINDOW},
  DR_GLOBAL_RATIO    => $p{DR_GLOBAL_RATIO},
  DR_SOLO_RATIO      => $p{DR_SOLO_RATIO},
  DR_LENGTH          => $p{DR_LENGTH},
  IR_LENGTH_MIN      => $p{IR_LENGTH_MIN},
  IR_LENGTH_MAX      => $p{IR_LENGTH_MAX},
  IR_WIDTH           => $p{IR_WIDTH},
);
_runscript($INSPECT, \%inspectopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


exit;


################################################################################
################################# SUBROUTINES  #################################
################################################################################
# Dump a string to a filepath.
#
sub _filedump
{
  my ($path, $data) = @_;
  open my $FH, '>', $path;
  print {$FH} $data;
  close $FH;
  return;
}

# given path to executable, a set of arguments as hash, and the file that should
# be created on success, run a command
#
sub _runscript
{
  my ($script, $args, $nextfile) = @_;
  my $command = _assemblecmd($script, $args);
  _writelog($command);
  $DNA->safeopen($command);

  if (! -e $nextfile)
  {
    my $errormsg = "\nDNA_ERROR: No $nextfile result! Dying.\n";
    _writelog($errormsg);
    croak $errormsg;
  }
  return 1;
}

# given a path to an executable and a set of arguments as a hash, generate a
# command
#
sub _assemblecmd
{
  my ($script, $args) = @_;
  my $command = $script;
  my @flags = sort keys %{$args};
  foreach my $flag (@flags)
  {
    $command .= q{ -} . $flag . q{ } . $args->{$flag};
  }
  return $command;
}

# given a list of tasks and a job name, start a bunch of taskfarmer workers
# and a client
#
sub _taskfarm
{
  my ($jobstr, $jobname) = @_;

  # Create the queue file that lists the tasks for the queue
  #
  my $queuefile = $FILEPREF . '_tasklist_' . $jobname . '.txt';
  _filedump($queuefile, $jobstr);
  push @CLEANUP, $queuefile;

  # Create worker request files that create workers for the queue.
  #
  my $workerfile = $FILEPREF . '_worker_' . $jobname . '.sh';
  my $queuename = $p{ANONKEY} . '_q_' . $jobname;
  my $workercmd = "\#!/bin/bash -l\nmpirun -n $p{SLOTS} ";
  $workercmd .= "tfmq-worker -q $queuename\n";
  _filedump($workerfile, $workercmd);
  push @CLEANUP, $workerfile;

  # Submit the request for the workers and capture the sge job id.
  #
  my $jname = $p{ANONKEY} . '_' . $jobname;
  my $workerqsub = "qsub $SGEFLAGS -N $jname ";
  $workerqsub .= '-o ' . $p{PROCDIR} . $jname . q{_o\$JOB_ID.txt };
  $workerqsub .= '-e ' . $p{PROCDIR} . $jname . q{_e\$JOB_ID.txt };
  $workerqsub .= $workerfile;
  _writelog($workerqsub);
  my $workerparse = $DNA->safeopen($workerqsub);
  _writelog($workerparse);

  # Start a client. Capture the job id and wait for completion.
  #
  my $clientcmd = "tfmq-client -i $queuefile -q $queuename";
  _writelog($clientcmd);
  $DNA->safeopen($clientcmd);

  return 1;
}

# write to the log
#
sub _writelog
{
  my ($data) = @_;
  open my $LOGH, '>>', $p{LOGFILE};
  print {$LOGH} "\n\n";
  print {$LOGH} $data;
  print {$LOGH} "\n\n";
  close $LOGH;
  return;
}

# dump a construct to an individual file
#
sub _constructdump
{
  my ($construct, $path) = @_;
  my $newdesign = Bio::DNAssemble::Design->new();
  $newdesign->add_construct($construct);
  $newdesign->dump_xml($path);
  return;
}

# merge two xml files to a newfile. Bring along old constructs if asked.
#
sub _combinefiles
{
  my ($oldpath, $pathsarr, $newpath, $oldswitch) = @_;
  my $olddesign = $DNA->load_design_from_file($oldpath);
  my $newdesign = $DNA->create_design(-name => $olddesign->name());
  
  my @addpaths = grep {-e $_} @{$pathsarr};
  my %skips = ();
  foreach my $docloc (@addpaths)
  {
    my $adddesign = $DNA->load_design_from_file($docloc);
    my @addconstructs = $adddesign->get_constructs();
    foreach my $addconstruct (@addconstructs)
    {
      $skips{$addconstruct->id()} = $addconstruct;
    }
    my @adderrors = $adddesign->get_errors();
    $newdesign->add_error($_) foreach (@adderrors);

    my @addwarnings = $adddesign->get_warnings();
    $newdesign->add_warning($_) foreach (@addwarnings);
  }
  
  my @oldconstructs = $olddesign->get_constructs();
  my %seens;
  foreach my $oldconstruct (@oldconstructs)
  {
    my $oldid = $oldconstruct->id();
    if (exists $skips{$oldid})
    {
      $newdesign->add_construct($skips{$oldid});
    }
    elsif ($oldswitch)
    {
      $newdesign->add_construct($oldconstruct);
    }
    $seens{$oldid}++;
  }
  foreach my $id (sort keys %skips)
  {
    next if (exists $seens{$id});
    $newdesign->add_construct($skips{$id});
  }
  $newdesign->dump_xml($newpath);
}

__END__

=head1 NAME

DNAssemble_004_ScreenedGBlock_SGE.pl

=head1 VERSION

1.00

=head1 DESCRIPTION

004_ScreenedGBlock runs the Prescreened GBlock design pipeline.
100-110-121-130-201-300-310

=head1 USAGE

=head1 OPTIONS

  You don't need to know.

=head1 DEPENDENCIES

If you are running this on genepool, you have erred.
Run DNAssemble_ScreenedGBlock.pl instead.

=head1 BUGS AND LIMITATIONS

NO BUGS! HAHAHAHAHAHA!

=head1 INCOMPATIBILITIES

As long as you use it EXACTLY like I showed you, and on genepool, there are no
incompatibilities.

=head1 AUTHOR

Sarah Richardson
SMRichardson@lbl.gov

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2014, DNAssemble developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of the Joint Genome Institute, the Lawrence Berkeley National
Laboratory, the Department of Energy, and the JGI developers may not be used to
endorse or promote products derived from this software without specific prior
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
