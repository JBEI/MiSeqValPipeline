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
my $DNAV = 'DNAssemble_117_GBlockable_' . $VERSION;
my $DNAS = '_117';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
  'input=s'                => \$p{INPUT},
  'output:s'               => \$p{OUTPUT},
  'logfile:s'              => \$p{LOGPATH},
  'gblock_length_min:i'    => \$p{GBLOCK_LENGTH_MIN},
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

################################################################################
################################## MARKING UP ##################################
################################################################################
my $count = 0;
my $ucount = 0;
my @constructs = $design->get_constructs();
foreach my $construct (@constructs)
{
  my $seqlen = $construct->seqlen();
  my @features = $construct->get_features(-primary_tag => 'ultramerable');
  my $mask = '0' x $seqlen;
  my %ublocks = ();
  if (! scalar @features)
  {
    my $gblock = Bio::SeqFeature::Generic->new(
      -start    => 1,
      -end      => $seqlen,
      -primary  => 'gblockable',
    );
    $construct->add_features([$gblock]);
  }
  foreach my $ultrablock (@features)
  {
    if ($seqlen - $ultrablock->end < $p{GBLOCK_LENGTH_MIN})
    {
      $ultrablock->end($seqlen);
    }
    if ($ultrablock->start - 1 < $p{GBLOCK_LENGTH_MIN})
    {
      $ultrablock->start(1);
    }
    my $ulen = $ultrablock->end - $ultrablock->start + 1;
    $ublocks{$ultrablock->start} = $ultrablock;
    substr $mask, $ultrablock->start - 1, $ulen, '1' x $ulen;
  }
  ##THIS IS WHERE WE FIX THE RANGES - SOMETIME IN THE FUTURE I HOPE
  my @ranges = @{find_deserts($mask)};
  if (scalar @ranges == 1)
  {
    my ($start, $end) = @{$ranges[0]};
    my $gblock = Bio::SeqFeature::Generic->new(
      -start    => $start,
      -end      => $end,
      -primary  => 'gblockable',
    );
    $construct->add_features([$gblock]);
  }
  else
  {
    foreach my $range (@ranges)
    {
      my ($start, $end) = @{$range};
      my $dlen = $end - $start + 1;
      if ($dlen < $p{GBLOCK_LENGTH_MIN})
      {
        print "Desert $start..$end is too small to use.. THIS SHOULD BE FIXED\n";
      }
    }
  }
  $count++;
}

print {$OUTH} "$count constructs screened.\n";

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
sub find_deserts
{
  my ($mask) = @_;
  my @ranges;

  my $len   = length $mask;
  my $init  = substr $mask, 0, 1;
  my $start = $init == 0 ? 1 : undef;
  my $end  = $init == 0 ? 0 : 1;

  for my $x (0 .. $len-1)
  {
    my $stat = substr $mask, $x, 1;
    #moving from feature to desert
    if ($stat == 0 && $end != 0)
    {
      $start = $x+1;
    }
    #moving from desert to feature
    elsif ($end == 0 && $stat != 0)
    {
      push @ranges, [$start, $x];
    }
    $end = $stat;
  }
  return \@ranges;
}

__END__

=head1 NAME

  DNAssemble_117_GBlockable.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

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