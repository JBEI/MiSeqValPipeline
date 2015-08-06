#!/usr/bin/env perl

use Bio::DNAssemble;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use Readonly;
use POSIX;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_115_Vector_Flank_' . $VERSION;
my $DNAS = '_115';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
  'input=s'                => \$p{INPUT},
  'output:s'               => \$p{OUTPUT},
  'destination_vector=s'   => \$p{DESTINATION_VECTOR},
  'logfile:s'              => \$p{LOGPATH},
  'overlap:i'              => \$p{OVERLAP},
  'help'                   => \$p{HELP}
);


################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>'NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE')
  if ($p{HELP});

my $DNA = Bio::DNAssemble->new();

my $OUTH;
if ($p{LOGPATH})
{
  open $OUTH, '>>', $p{LOGPATH};
}
else
{
  $OUTH = *STDOUT;
}

print {$OUTH} "\n\n******* $DNAV WORKING\n\n";

#The input file must exist and be a format we care to read.
my $design = $DNA->load_design_from_file($p{INPUT});
my $filename = $design->filename();
$filename = $1 if ($filename =~ m{(.+)\_}msix);
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.xml};

#The vector must be provided.
if (! $p{DESTINATION_VECTOR})
{
  print {$OUTH} "\nDNA_ERROR: vector not provided.\n";
  print {$OUTH} "\n\n******* $DNAV BAILED\n\n";
  close $OUTH;
  croak;
}

#An overlap length must be provided
if (! $p{OVERLAP})
{
  print {$OUTH} "\nDNA_ERROR: overlap not provided.\n";
  print {$OUTH} "\n\n******* $DNAV BAILED\n\n";
  close $OUTH;
  croak;
}

my $destvec = $DNA->load_vector(-name => $p{DESTINATION_VECTOR});
my $fiveflank = $destvec->chew5($p{OVERLAP});
my $threeflank = $destvec->chew3($p{OVERLAP});

################################################################################
############################### VECTOR FLANKING  ###############################
################################################################################
my $count = 0;
my @constructs = $design->get_constructs();
foreach my $construct (@constructs)
{
  my $currseq = $construct->sequence();
  my $newseq = $fiveflank . $currseq . $threeflank;
  $construct->sequence($newseq);
  $construct->vector($p{DESTINATION_VECTOR});
  my @features = $construct->get_features();
  foreach my $feature (@features)
  {
    my $newstart = $feature->start + $p{OVERLAP};
    my $newend = $feature->end + $p{OVERLAP};
    $feature->start($newstart);
    $feature->end($newend);
  }
  $count++;
}

print {$OUTH} "$count constructs nicknamed.\n";

################################################################################
################################## REPORTING ###################################
################################################################################
$design->dump_xml($p{OUTPUT});

print {$OUTH} "\n\n";
print {$OUTH} "Wrote $p{OUTPUT}\n\n";
print {$OUTH} $DNA->attitude() . " brought to you by $DNAV\n\n";

print {$OUTH} "\n\n******* $DNAV FINISHED\n\n";
close $OUTH;
exit;


################################################################################
################################# SUBROUTINES ##################################
################################################################################

__END__

=head1 NAME

  DNAssemble_115_Vector_Flank.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

    -d,  --destination_vector : what vector to use

  Optional arguments:

    -o,  --output : a filepath for dumping output

    -l,  --logfile : A logfile to update

    -h,  --help : display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, GeneDesign developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of Johns Hopkins, the Joint Genome Institute, the Lawrence Berkeley
National Laboratory, the Department of Energy, and the GeneDesign developers may
not be used to endorse or promote products derived from this software without
specific prior written permission.

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