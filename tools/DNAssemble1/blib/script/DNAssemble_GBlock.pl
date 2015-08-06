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
my $DNAV = 'DNAssemble_GBlock_' . $VERSION;
my $DNAS = '_000';

local $OUTPUT_AUTOFLUSH = 1;

# Get command line arguments
#
my %p = ();
GetOptions (
  'help'            => \$p{HELP},
  'input=s'         => \$p{INPUT},
  'cloningvector=s' => \$p{CLONVEC},
  'anonkey=s'       => \$p{ANONKEY},
  'name:s'          => \$p{NAME},
  'assemblyTm:i'    => \$p{ASSEMBLYTM},

  'email:s'         => \$p{EMAIL},
  'slots:i'         => \$p{SLOTS},

  'bblength:i'      => \$p{TARBBLEN},
  'bbmaxlength:i'   => \$p{MAXBBLEN},
  'bbminlength:i'   => \$p{MINBBLEN},
  'bblap:i'         => \$p{TARBBLAP},
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
Readonly my $TARBBLEN   =>  950;
Readonly my $MAXBBLEN   => 1000;
Readonly my $MINBBLEN   =>  140;
Readonly my $TARBBLAP   =>   60;
Readonly my $MINBBLAP   =>   50;
Readonly my $SLOTS      =>   16;

Readonly my $BIN        => '/usr/common/jgi/synthbio/DNAssemble/bin/';
Readonly my $GBLOCK     => $BIN . 'DNAssemble_003_GBlock_SGE.pl';

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

# A cloning vector must be defined.
# It must be accessible to genedesign. I am sorry about that, I
# will work on that next.
#
croak "\nDNA_ERROR: You must use a cloning vector\n\n" if (! $p{CLONVEC});
my $clonvec = $DNA->load_vector(-name => $p{CLONVEC});

# Building block size parameters and sanity check
#
$p{TARBBLEN} = $p{TARBBLEN} || $TARBBLEN;
$p{MAXBBLEN} = $p{MAXBBLEN} || $MAXBBLEN;
$p{MINBBLEN} = $p{MINBBLEN} || $MINBBLEN;
$p{TARBBLAP} = $p{TARBBLAP} || $TARBBLAP;
croak "\n DNA_ERROR: building block size is outside of allowable range.\n"
  if ($p{TARBBLEN} < $p{MINBBLEN} || $p{TARBBLEN} > $p{MAXBBLEN});
croak "\n DNA_ERROR: chewback overlap is too small.\n"
  if ($p{TARBBLAP} < $MINBBLAP);


################################################################################
############################## SUBMIT THE COMMAND ##############################
################################################################################
my %classicopts = (
  input         => $p{INPUT},
  anonkey       => $p{ANONKEY},
  procdir       => $PROCDIR,
  logfile       => $LOGPATH,
  assemblyTm    => $p{ASSEMBLYTM},
  slots         => $p{SLOTS},
  name          => $NAME,

  cloningvector => $p{CLONVEC},
  bblength      => $p{TARBBLEN},
  bbmaxlength   => $p{MAXBBLEN},
  bbminlength   => $p{MINBBLEN},
  bblap         => $p{TARBBLAP},
);

my $command = $GBLOCK;
my @flags = sort keys %classicopts;
foreach my $flag (@flags)
{
  $command .= q{ -} . $flag . q{ } . $classicopts{$flag};
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

DNAssemble_GBlock.pl

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

  --cloningvector: A vector known to GeneDesign that is to serve as the
          cloning vector

  --anonkey: What anonymization key are we using? 150 characters max,
          no spaces, periods, or underscores

  --assemblyTm: The target melting temperature of assembly primers
          Usually we use 70


=head1 OPTIONS

  --name  : A name for the design. If not supplied, it defaults to the anonkey.

  --email : Your email, if you want genepool to tell you when the job is done.

  --slots : How many nodes on Genepool to use for distributed jobs.
          Defaults to 16.

  --help : Display this message

Options for Building Blocks:

  --BBlength: the length in bp of building blocks
          Defaults to 950

  --BBmaxlength: the maximum length a building block is allowed to be
          Defaults to 1000

  --BBminlength: the minimum length a building block is allowed to be
          Defaults to 140

  --BBlap: the target overlap between building blocks
          Defaults to 60, minimum is 50


=head1 DEPENDENCIES

If you are running this on genepool, you must load the modules usearch and
dnassemble first. If you didn't do that, I am surprised that you can read this.

This application requires taskfarmermq v2.1, which requires the python module.

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
