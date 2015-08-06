#!/usr/bin/env perl

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
my $DNAV = 'DNAssemble_001_Classic_Local_' . $VERSION;
my $DNAS = '_002';

local $OUTPUT_AUTOFLUSH = 1;

# Get command line arguments
#
my %p = ();
GetOptions (
  'help'            => \$p{HELP},
  'input=s'         => \$p{INPUT},
  'cloningvector=s' => \$p{CLONVEC},
  'destvector=s'    => \$p{DESTVEC},
  'assemblyTm:i'    => \$p{ASSEMBLYTM},
  'screenpi:f'      => \$p{PERCID},
  'name:s'          => \$p{NAME},

  'anonkey=s'       => \$p{ANONKEY},
  'procdir=s'       => \$p{PROCDIR},
  'logfile=s'       => \$p{LOGFILE},

  'bblength:i'      => \$p{TARBBLEN},
  'bbmaxlength:i'   => \$p{MAXBBLEN},
  'bbminlength:i'   => \$p{MINBBLEN},
  'bblap:i'         => \$p{TARBBLAP},

  'ollength:i'      => \$p{TAROLLEN},
  'olmaxlength:i'   => \$p{MAXOLLEN},
  'olminlength:i'   => \$p{MINOLLEN},
  'poolsize:i'      => \$p{POOLSIZE},
  'maxpoolnum:i'    => \$p{MAXPOOLNUM},
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
Readonly my $MINBBLAP   =>   10;

Readonly my $CLEAN      => 'DNAssemble_100_CreateXML.pl';
Readonly my $COGNOM     => 'DNAssemble_110_Cognominate.pl';
Readonly my $CARVE      => 'DNAssemble_120_Carve.pl';
Readonly my $COMPRESSA  => 'DNAssemble_130_Compress1.pl';
Readonly my $CLUSTER    => 'DNAssemble_140_Cluster.pl';
Readonly my $CHOP       => 'DNAssemble_160_Chop.pl';
Readonly my $COMPRESSB  => 'DNAssemble_170_Compress2.pl';
Readonly my $CHECK      => 'DNAssemble_200_Check.pl';
Readonly my $SUMMARY    => 'DNAssemble_300_Summary.pl';
Readonly my $STITCHERS  => 'DNAssemble_320_Stitchers.pl';

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
################################### CLEANING ###################################
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
################################### CARVING  ###################################
################################################################################
$NEXTFILE = $FILEPREF . q{_120.xml};
my %carveopts = (
  input             => $LASTFILE,
  output            => $NEXTFILE,
  logfile           => $p{LOGFILE},
  cloningvector     => $p{CLONVEC},
  destinationvector => $p{DESTVEC},
  length            => $p{TARBBLEN},
  maxlength         => $p{MAXBBLEN},
  minlength         => $p{MINBBLEN},
  lap               => $p{TARBBLAP},
  stitch            => $p{ASSEMBLYTM},
);
_runscript($CARVE, \%carveopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;

################################################################################
############################## FIRST  COMPRESSING ##############################
################################################################################
$NEXTFILE = $FILEPREF . q{_130.xml};
my %fcompressopts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  logfile  => $p{LOGFILE},
);
_runscript($COMPRESSA, \%fcompressopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
################################## CLUSTERING ##################################
################################################################################
$NEXTFILE = $FILEPREF . q{_140.xml};
my %clusteropts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  logfile  => $p{LOGFILE},
  percid   => $p{PERCID},
);
_runscript($CLUSTER, \%clusteropts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
################################### CHOPPING ###################################
################################################################################
$NEXTFILE = $FILEPREF . q{_160.xml};
my %chopopts = (
  input       => $LASTFILE,
  output      => $NEXTFILE,
  logfile     => $p{LOGFILE},
  olilen      => $p{TAROLLEN},
  maxolilen   => $p{MAXOLLEN},
  minolilen   => $p{MINOLLEN},
  laptemp     => $p{ASSEMBLYTM},
  poolsize    => $p{POOLSIZE},
  maxpoolnum  => $p{MAXPOOLNUM},
);
_runscript($CHOP, \%chopopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
############################## SECOND COMPRESSING ##############################
################################################################################
$NEXTFILE = $FILEPREF . q{_170.xml};
my %scompressopts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  logfile  => $p{LOGFILE},
);
_runscript($COMPRESSB, \%scompressopts, $NEXTFILE);
$LASTFILE = $NEXTFILE;


################################################################################
################################### CHECKING ###################################
################################################################################
$NEXTFILE = $FILEPREF . q{_200.xml};
my %checkopts = (
  input    => $LASTFILE,
  output   => $NEXTFILE,
  logfile  => $p{LOGFILE},
);
_runscript($CHECK, \%checkopts, $NEXTFILE);

################################################################################
################################# SUMMARIZING  #################################
################################################################################
my $summary = $p{PROCDIR} . $p{ANONKEY} . q{_Summary.txt};
my %summaryopts = (
  input    => $NEXTFILE,
  output   => $summary,
  logfile  => $p{LOGFILE},
);
_runscript($SUMMARY, \%summaryopts, $summary);


################################################################################
################################## STITCHING  ##################################
################################################################################
my $stitch = $p{PROCDIR} . $p{ANONKEY} . q{_StitchingOligos.xls};
my %stichingopts = (
  input    => $NEXTFILE,
  output   => $stitch,
  logfile  => $p{LOGFILE},
);
_runscript($STITCHERS, \%stichingopts, $stitch);


################################################################################
################################## CONCLUDING ##################################
################################################################################
$DNA->cleanup(\@CLEANUP);

my $closingremarks = q{};
my $finalxml = $p{PROCDIR} . $p{ANONKEY} . q{_Design.xml};
$DNA->safeopen("cp $NEXTFILE $finalxml");
$closingremarks .= "Wrote $finalxml\n";
$closingremarks .= "\n\n" . $DNA->attitude() . " brought to you by $DNAV\n\n";
_writelog($closingremarks);


exit;


################################################################################
################################# SUBROUTINES  #################################
################################################################################
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


__END__

=head1 NAME

DNAssemble_001_Classic_Local.pl

=head1 VERSION

1.00

=head1 DESCRIPTION

001_Local runs the classic design pipeline.
100-110-120-130-140-160-170-200

=head1 USAGE

=head1 OPTIONS

  You don't need to know.

=head1 DEPENDENCIES

If you are running this on genepool, you have erred. Run DNAssemble_Classic.pl
instead.

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
