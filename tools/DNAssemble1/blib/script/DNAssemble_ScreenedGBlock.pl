#!/usr/bin/env perl

use Bio::DNAssemble;
use Getopt::Long;
use Pod::Usage;
use Readonly;
use English qw( -no_match_vars );
use autodie qw(open close mkdir);
use Carp;
use Env qw(SCRATCH);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_ScreenedGBlock_' . $VERSION;
my $DNAS = '_000';

local $OUTPUT_AUTOFLUSH = 1;

# Get command line arguments
#
my %p = ();
GetOptions (
  'help'                 => \$p{HELP},
  'INPUT=s'              => \$p{INPUT},
  'DESTINATION_VECTOR=s' => \$p{DESTINATION_VECTOR},
  'ANONKEY=s'            => \$p{ANONKEY},
  'NAME:s'               => \$p{NAME},
  'ASSEMBLYTM:i'         => \$p{ASSEMBLYTM},
                         
  'EMAIL:s'              => \$p{EMAIL},
  'SLOTS:i'              => \$p{SLOTS},

  'SKIP_WARNINGS'        => \$p{SKIP_WARNINGS},
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

# Non user configurable variables
# I also define the scripts to be used at the different stages
# here. I don't like hard coding $BIN. I will talk to Kirsten and Doug about
# that.
#
Readonly my $SLOTS      =>   16;

Readonly my $BIN        => '/usr/common/jgi/synthbio/DNAssemble/bin/';
Readonly my $GBLOCK     => $BIN . 'DNAssemble_004_ScreenedGBlock_SGE.pl';

$p{SLOTS} = $p{SLOTS} || $SLOTS;

# The input file must exist. We will not do deep parsing but it must pass an
# initial GeneDesign load to prove formatting.
#
croak "\nDNA_ERROR: You must supply an input file.\n\n" if (! $p{INPUT});
croak "\nDNA_ERROR: $p{INPUT} does not exist.\n\n" if (! -e $p{INPUT});
my ($iterator, $filename, $filesuffix) = $DNA->import_seqs($p{INPUT});

# The anon name must be provided and be less than 150 characters. No spaces,
# commmas, or periods are allowed.
#
croak "\nDNA_ERROR: No anonymization key given!\n\n" if (! $p{ANONKEY});
if ($p{ANONKEY} =~ /[\_\,\.\s]/)
{
  croak "\nDNA_ERROR: anonymization key has bad characters\n\n";
}
if (length $p{ANONKEY} > 150)
{
  croak "\nDNA_ERROR: anonymization key is too long\n\n";
}
my $NAME = $p{NAME} || $p{ANONKEY};

# We create a processing directory inside a temporary directory. We will try to
# use $SCRATCH, otherwise it will be created in the cwd. The name will be tied
# to the anonymization key. We will also start a logfile.
#
my $tmpdir = $SCRATCH || q{./};
$tmpdir = $DNA->endslash($tmpdir);
my $PROCDIR = $tmpdir . $p{ANONKEY} . q{/};
mkdir $PROCDIR if (! -e $PROCDIR);
my $LOGPATH = $PROCDIR . $p{ANONKEY} . q{_LOG.txt};

# The assembly Tm is the melting temperature of all oligos in the design.
# I should probably check for sanity here but if I can't trust you, who can I
# trust?
#
croak "\nDNA_ERROR: No assembly Tm defined!\n\n" if (! $p{ASSEMBLYTM});

# A destination vector must be defined.
# It must be accessible to dnassemble. I am sorry about that, I
# will work on that next.
#
croak "\nDNA_ERROR: You must name a destination vector\n\n" if (! $p{DESTINATION_VECTOR});
my $destvec = $DNA->load_vector(-name => $p{DESTINATION_VECTOR});

# GBlock and oligo length parameters and sanity check
#

Readonly my $D_GBLOCK_LENGTH_MAX => 1000;
$p{GBLOCK_LENGTH_MAX} = $p{GBLOCK_LENGTH_MAX} || $D_GBLOCK_LENGTH_MAX;
Readonly my $D_GBLOCK_LENGTH_MIN => 140;
$p{GBLOCK_LENGTH_MIN} = $p{GBLOCK_LENGTH_MIN} || $D_GBLOCK_LENGTH_MIN;
Readonly my $D_GBLOCK_LENGTH => 980;
$p{GBLOCK_LENGTH} = $p{GBLOCK_LENGTH} || $D_GBLOCK_LENGTH;
Readonly my $D_GBLOCK_MINOVERLAP => 25;
Readonly my $D_GBLOCK_OVERLAP => 40;
$p{GBLOCK_OVERLAP} = $p{GBLOCK_OVERLAP} || $D_GBLOCK_OVERLAP;

croak "\n DNA_ERROR: gblock size is outside of allowable range.\n"
  if ($p{GBLOCK_LENGTH} < $p{GBLOCK_LENGTH_MIN} || $p{GBLOCK_LENGTH} > $p{GBLOCK_LENGTH_MAX});
croak "\n DNA_ERROR: chewback overlap is too small.\n"
  if ($p{GBLOCK_OVERLAP} < $D_GBLOCK_MINOVERLAP);
  
Readonly my $D_OLIGO_LENGTH_MAX => 200;
$p{OLIGO_LENGTH_MAX} = $p{OLIGO_LENGTH_MAX} || $D_OLIGO_LENGTH_MAX;
Readonly my $D_OLIGO_LENGTH_MIN => 45;
$p{OLIGO_LENGTH_MIN} = $p{OLIGO_LENGTH_MIN} || $D_OLIGO_LENGTH_MIN;
Readonly my $D_OLIGO_LENGTH => 150;
$p{OLIGO_LENGTH} = $p{OLIGO_LENGTH} || $D_OLIGO_LENGTH;
Readonly my $D_OLIGO_MINOVERLAP => 10;
Readonly my $D_OLIGO_OVERLAP    => 40;
$p{OLIGO_OVERLAP} = $p{OLIGO_OVERLAP} || $D_OLIGO_OVERLAP;

croak "\n DNA_ERROR: oligo size is outside of allowable range.\n"
  if ($p{OLIGO_LENGTH} < $p{OLIGO_LENGTH_MIN} || $p{OLIGO_LENGTH} > $p{OLIGO_LENGTH_MAX});
croak "\n DNA_ERROR: chewback overlap is too small.\n"
  if ($p{OLIGO_OVERLAP} < $D_OLIGO_MINOVERLAP);

# GBlock sequence parameters and sanity check
#
$p{SKIP_WARNINGS} = $p{SKIP_WARNINGS} || 1;

## GLOBAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_GLOBAL_MIN => 25;
Readonly my $D_GC_GLOBAL_MAX => 69;
$p{GC_GLOBAL_MIN} = $p{GC_GLOBAL_MIN} || $D_GC_GLOBAL_MIN;
$p{GC_GLOBAL_MAX} = $p{GC_GLOBAL_MAX} || $D_GC_GLOBAL_MAX;
die 'DNA_ERROR: Global GC percentage arguments do not parse.'
  if ($p{GC_GLOBAL_MIN} >= $p{GC_GLOBAL_MAX});

## REGIONAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_REGIONAL_WINDOW => 100;
$p{GC_REGIONAL_WINDOW} = $p{GC_REGIONAL_WINDOW} || $D_GC_REGIONAL_WINDOW;
Readonly my $D_GC_REGIONAL_MIN => 28;
Readonly my $D_GC_REGIONAL_MAX => 77;
$p{GC_REGIONAL_MIN} = $p{GC_REGIONAL_MIN} || $D_GC_REGIONAL_MIN;
$p{GC_REGIONAL_MAX} = $p{GC_REGIONAL_MAX} || $D_GC_REGIONAL_MAX;
die 'DNA_ERROR: Regional GC percentage arguments do not parse.'
  if ($p{GC_REGIONAL_MIN} >= $p{GC_REGIONAL_MAX});

## LOCAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_LOCAL_WINDOW => 20;
$p{GC_LOCAL_WINDOW} = $p{GC_LOCAL_WINDOW} || $D_GC_LOCAL_WINDOW;
Readonly my $D_GC_LOCAL_MIN => 15;
Readonly my $D_GC_LOCAL_MAX => 90;
$p{GC_LOCAL_MIN} = $p{GC_LOCAL_MIN} || $D_GC_LOCAL_MIN;
$p{GC_LOCAL_MAX} = $p{GC_LOCAL_MAX} || $D_GC_LOCAL_MAX;
die 'DNA_ERROR: Local GC percentage arguments do not parse.'
  if ($p{GC_LOCAL_MIN} >= $p{GC_LOCAL_MAX});
die 'DNA_ERROR: Regional and local window sizes do not parse.'
  if ($p{GC_LOCAL_WINDOW} >= $p{GC_REGIONAL_WINDOW});

## TERMINAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_TERMINAL_WINDOW => 30;
$p{GC_TERMINAL_WINDOW} = $p{GC_TERMINAL_WINDOW} || $D_GC_TERMINAL_WINDOW;
Readonly my $D_GC_TERMINAL_MIN => 24;
Readonly my $D_GC_TERMINAL_MAX => 76;
$p{GC_TERMINAL_MIN} = $p{GC_TERMINAL_MIN} || $D_GC_TERMINAL_MIN;
$p{GC_TERMINAL_MAX} = $p{GC_TERMINAL_MAX} || $D_GC_TERMINAL_MAX;
die 'DNA_ERROR: Terminal GC percentage arguments do not parse.'
  if ($p{GC_TERMINAL_MIN} >= $p{GC_TERMINAL_MAX});

## HOMOPOLYMER LENGTHS ARGUMENTS
Readonly my $HP_MIN_LENGTH => 5;
Readonly my $D_HP_A_LENGTH => 12;
$p{HP_A_LENGTH} = $p{HP_A_LENGTH} || $D_HP_A_LENGTH;
die 'DNA_ERROR: The homopolymer limit for A is too low.'
  if ($p{HP_A_LENGTH} < $HP_MIN_LENGTH);
Readonly my $D_HP_T_LENGTH => 12;
$p{HP_T_LENGTH} = $p{HP_T_LENGTH} || $D_HP_T_LENGTH;
die 'DNA_ERROR: The homopolymer limit for T is too low.'
  if ($p{HP_T_LENGTH} < $HP_MIN_LENGTH);
Readonly my $D_HP_C_LENGTH => 8;
$p{HP_C_LENGTH} = $p{HP_C_LENGTH} || $D_HP_C_LENGTH;
die 'DNA_ERROR: The homopolymer limit for C is too low.'
  if ($p{HP_C_LENGTH} < $HP_MIN_LENGTH);
Readonly my $D_HP_G_LENGTH => 8;
$p{HP_G_LENGTH} = $p{HP_G_LENGTH} || $D_HP_G_LENGTH;
die 'DNA_ERROR: The homopolymer limit for G is too low.'
  if ($p{HP_G_LENGTH} < $HP_MIN_LENGTH);

## REPEAT ARGUMENTS
Readonly my $D_DIMER_LIMIT => 9;
$p{DIMER_LIMIT} = $p{DIMER_LIMIT} || $D_DIMER_LIMIT;
Readonly my $D_TRIMER_LIMIT => 5;
$p{TRIMER_LIMIT} = $p{TRIMER_LIMIT} || $D_TRIMER_LIMIT;

##LOCAL REPEAT ARGUMENTS
Readonly my $D_DR_LOCAL_WINDOW => 70;
$p{DR_LOCAL_WINDOW} = $p{DR_LOCAL_WINDOW} || $D_DR_LOCAL_WINDOW;
Readonly my $D_DR_LOCAL_RATIO => 90;
$p{DR_LOCAL_RATIO} = $p{DR_LOCAL_RATIO} || $D_DR_LOCAL_RATIO;

## REGIONAL REPEAT ARGUMENTS
Readonly my $D_DR_REGIONAL_WINDOW => 500;
$p{DR_REGIONAL_WINDOW} = $p{DR_REGIONAL_WINDOW} || $D_DR_REGIONAL_WINDOW;
Readonly my $D_DR_REGIONAL_RATIO => 60;
$p{DR_REGIONAL_RATIO} = $p{DR_REGIONAL_RATIO} || $D_DR_REGIONAL_RATIO;

Readonly my $D_DR_GLOBAL_RATIO => .69;
$p{DR_GLOBAL_RATIO} = $p{DR_GLOBAL_RATIO} || $D_DR_GLOBAL_RATIO;
Readonly my $D_DR_SOLO_RATIO => .4;
$p{DR_SOLO_RATIO} = $p{DR_SOLO_RATIO} || $D_DR_SOLO_RATIO;
Readonly my $D_DR_LENGTH => 8;
$p{DR_LENGTH} = $p{DR_LENGTH} || $D_DR_LENGTH;
Readonly my $D_IR_LENGTH_MIN => 16;
$p{IR_LENGTH_MIN} = $p{IR_LENGTH_MIN} || $D_IR_LENGTH_MIN;
Readonly my $D_IR_LENGTH_MAX => 19;
$p{IR_LENGTH_MAX} = $p{IR_LENGTH_MAX} || $D_IR_LENGTH_MAX;
Readonly my $D_IR_WIDTH => 100;
$p{IR_WIDTH} = $p{IR_WIDTH} || $D_IR_WIDTH;


################################################################################
############################## SUBMIT THE COMMAND ##############################
################################################################################
my %screengblockopts = (
  input              => $p{INPUT},
  anonkey            => $p{ANONKEY},
  procdir            => $PROCDIR,
  logfile            => $LOGPATH,
  assemblyTm         => $p{ASSEMBLYTM},
  slots              => $p{SLOTS},
  name               => $NAME,
  DESTINATION_VECTOR => $p{DESTINATION_VECTOR},
  
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
  DR_LOCAL_RATIO     => $p{DR_LOCAL_RATIO},
  DR_REGIONAL_WINDOW => $p{DR_REGIONAL_WINDOW},
  DR_REGIONAL_RATIO  => $p{DR_REGIONAL_RATIO},
  DR_GLOBAL_RATIO    => $p{DR_GLOBAL_RATIO},
  DR_SOLO_RATIO      => $p{DR_SOLO_RATIO},
  DR_LENGTH          => $p{DR_LENGTH},
  IR_LENGTH_MIN      => $p{IR_LENGTH_MIN},
  IR_LENGTH_MAX      => $p{IR_LENGTH_MAX},
  IR_WIDTH           => $p{IR_WIDTH},

  GBLOCK_LENGTH_MIN  => $p{GBLOCK_LENGTH_MIN},
  GBLOCK_LENGTH_MAX  => $p{GBLOCK_LENGTH_MAX},
  GBLOCK_LENGTH      => $p{GBLOCK_LENGTH},
  GBLOCK_OVERLAP     => $p{GBLOCK_OVERLAP},
  GBLOCK_OVERLAP_MIN => $D_GBLOCK_MINOVERLAP,
  OLIGO_LENGTH_MIN   => $p{OLIGO_LENGTH_MIN},
  OLIGO_LENGTH_MAX   => $p{OLIGO_LENGTH_MAX},
  OLIGO_LENGTH       => $p{OLIGO_LENGTH},
  OLIGO_OVERLAP_MIN  => $p{OLIGO_OVERLAP_MIN},
);

my $command = $GBLOCK;
my @flags = sort keys %screengblockopts;
foreach my $flag (@flags)
{
  $command .= q{ -} . $flag . q{ } . $screengblockopts{$flag};
}

my $submit = 'qsub -N ' . $p{ANONKEY} . q{ };
$submit .= '-M ' . $p{EMAIL} . ' -m abe ' if ($p{EMAIL});
$submit .= '-o ' . $PROCDIR . $p{ANONKEY} . q{_o\$JOB_ID.txt };
$submit .= '-e ' . $PROCDIR . $p{ANONKEY} . q{_e\$JOB_ID.txt };
$submit .= $command;
print "$submit\n\n";
my $parse = $DNA->safeopen($submit);
print "$parse\n\n";

open my $LOGH, '>', $LOGPATH;
print {$LOGH} "\n\n";
print {$LOGH} $submit;
print {$LOGH} "\n\n";
print {$LOGH} $parse;
print {$LOGH} "\n\n";
close $LOGH;


exit;

__END__

=head1 NAME

DNAssemble_ScreenedGBlock.pl

=head1 VERSION

1.00

=head1 DESCRIPTION

Carves a set of chunks into building blocks, all embarassingly parallel like.

=head1 USAGE

First load the module
module load dnassemble

=head1 REQUIRED ARGUMENTS

  --input : A sequence file to use as an input. This file may contain one
          or more construct sequences that are destined for building block and
          oligo design.
  --destination_vector: A vector known to DNAssemble that is to serve as the
          destination vector
  --anonkey: What anonymization key are we using? 150 characters max,
          no spaces, periods, or underscores
  --assemblyTm: The target melting temperature of assembly primers
          Usually we use 70


=head1 OPTIONS

  --name  : A name for the design. If not supplied, it defaults to the anonkey.
  --email : Your email, if you want genepool to tell you when the job is done.
  --slots : How many nodes on Genepool to use for distributed jobs (16)
  --help : Display this message

Options for GBlocks sequence content:

  --SKIP_WARNINGS : Do not report things that will not cause outright rejection
                    default is 0 - will not skip
  --GC_GLOBAL_MIN : Minimum GC percentage allowable in a gblock (25)
  --GC_GLOBAL_MAX : Maximum GC percentage allowable in a gblock (69)
  --GC_REGIONAL_WINDOW : Regional window size for GC calculations (100)
  --GC_REGIONAL_MIN : Minimum GC percentage allowable in a regional window (28)
  --GC_REGIONAL_MAX : Maximum GC percentage allowable in a regional window (77)
  --GC_LOCAL_WINDOW : Local window size for GC percentage calculations (20)
  --GC_LOCAL_MIN : Minimum GC percentage allowable, over any local window (15)
  --GC_LOCAL_MAX : Maximum GC percentage allowable, over any local window (90)
  --GC_TERMINAL_WINDOW : Number of bases at each end to check for GC (20)
  --GC_TERMINAL_MIN : Minimum GC percentage allowable in terminals (24)
  --GC_TERMINAL_MAX : Maximum GC percentage allowable in terminals (76)
  --HP_A_LENGTH   : Maximum allowable length of an A homopolymer (12)
  --HP_T_LENGTH   : Maximum allowable length of a T homopolymer (12)
  --HP_C_LENGTH   : Maximum allowable length of a C homopolymer (8)
  --HP_G_LENGTH   : Maximum allowable length of a G homopolymer (8)
  --DIMER_LIMIT   : Maximum number of direct dimer repeats allowable (9)
  --TRIMER_LIMIT  : Maximum number of direct dimer repeats allowable (9)
  --DR_LOCAL_WINDOW : Window for direct repeat calculations (70)
  --DR_LOCAL_RATIO: Percent of bases involved in repeats in a local window (90)
  --DR_REGIONAL_WINDOW : Window for direct repeat calculations (500)
  --DR_REGIONAL_RATIO : Percent bases involved in repeats in a regional window (60)
  --DR_GLOBAL_RATIO : Percent bases involved in repeats (69)
  --DR_SOLO_RATIO : Percent of bases involved in any iteration of a repeat (.4)
  --DR_LENGTH     : Length of a direct repeat that is cause for concern (8)
  --IR_LENGTH_MIN : Length of an inverted repeat that could be cause for concern (16)
  --IR_LENGTH_MAX : Maximum allowable length of an inverted repeat (19)
  --IR_WIDTH      : Window size for inverted repeat (100)

Options for GBlock and oligo sizes:

  --GBLOCK_LENGTH_MIN  : Minimum allowable length for a gblock (140)
  --GBLOCK_LENGTH_MAX  : Maximum allowable length for a gblock (1000)
  --GBLOCK_LENGTH      : Target length for a gblock (980)
  --GBLOCK_OVERLAP     : Length of overlap with each other and vector (40)
  --OLIGO_LENGTH_MIN   : Minimum allowable length for an oligo (45)
  --OLIGO_LENGTH_MAX   : Maximum allowable length for an oligo (200)
  --OLIGO_LENGTH       : Target length for an oligo (150)
  --OLIGO_OVERLAP_MIN  : Minimum overlap between oligos (10)

=head1 DEPENDENCIES

If you are running this on genepool, you must load the module dnassemble first.
If you didn't do that, I am surprised that you can read this.

This application requires taskfarmermq v2.4, which requires the python module.

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
