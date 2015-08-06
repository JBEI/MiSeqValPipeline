#!/usr/bin/env perl

use Bio::DNAssemble;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_310_Extract_' . $VERSION;
my $DNAS = '_300';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'             => \$p{INPUT},
      'output:s'            => \$p{OUTPUT},
      'type=s'              => \$p{TYPE},
      'logfile:s'           => \$p{LOGPATH},
      'help'                => \$p{HELP}
);


################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>'NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE')
  if ($p{HELP});

my $DNA = Bio::DNAssemble->new();
my $GD = $DNA->GD();

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
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.fasta};

#The user must specify some types.
my %GTYPES = map {$_ => 1} qw(assembly_oligo stitching_oligo building_block
  intermediate chunk);
croak "\n DNA_ERROR: no construct types named for extraction.\n"
  if (! $p{TYPE});
my @types = split q{,}, $p{TYPE};
@types = grep {exists $GTYPES{$_}} @types;
croak "\n DNA_ERROR: no valid construct types named for extraction.\n"
  if (! scalar @types);


################################################################################
################################# EXTRACTING  ##################################
################################################################################
my @report;
foreach my $type (@types)
{
  push @report, $design->constructs_to_fasta(-kind => $type, -sort => 'name');
}


################################################################################
################################## REPORTING ###################################
################################################################################
open my $REPH, '>', $p{OUTPUT};
print $REPH @report;
close $REPH;

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

  DNAssemble_310_Extract.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

    -t,  --type : a comma separated list of types to pull out.
          assembly_oligo,stitching_oligo,building_block,intermediate,chunk

  Optional arguments:

    -o,  --output : a filepath for dumping output

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