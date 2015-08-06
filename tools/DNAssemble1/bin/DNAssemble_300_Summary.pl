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
my $DNAV = 'DNAssemble_300_Summary_' . $VERSION;
my $DNAS = '_300';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'             => \$p{INPUT},
      'output:s'            => \$p{OUTPUT},
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
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.txt};


################################################################################
################################# SUMMARIZING ##################################
################################################################################

my @report;

# Pull all of the errors
#
my $errors = $design->errors_as_hash();
foreach my $sid (sort keys %{$errors})
{
  push @report, "ERROR for $sid : " . $errors->{$sid} . "\n";
}

# Pull all of the warnings
#
my $warnings = $design->warnings_as_hash();
foreach my $sid (sort keys %{$warnings})
{
  push @report, "WARNING for $sid : " . $warnings->{$sid} . "\n";
}

push @report, "\n\n";

# Pull all of the assembly oligos
#
my @cons = $design->get_constructs(-kind => 'assembly_oligo');
my $aolcount = 0;
my $aolbasecount = 0;
my ($minaollen, $maxaollen) = (undef, undef);
foreach my $aol (@cons)
{
  $aolcount++;
  my $seqlen = $aol->seqlen();
  $minaollen = $seqlen if (! defined $minaollen || $seqlen < $minaollen);
  $maxaollen = $seqlen if (! defined $maxaollen || $seqlen > $maxaollen);
  $aolbasecount +=  $aol->seqlen();
}
push @report, "Assembly Oligos:\n";
push @report, "\t$aolcount oligos ($aolbasecount bp)\n";
if ($aolcount > 0)
{
  my $avgaollen = int ($aolbasecount / $aolcount);
  push @report, "\t\tAverage length: $avgaollen\n";
  push @report, "\t\tMinimum length: $minaollen\n";
  push @report, "\t\tMaximum length: $maxaollen\n";
}
push @report,  "\n\n";

# Pull all of the stitching oligos
#
@cons = $design->get_constructs(-kind => 'stitching_oligo');
my $solcount = 0;
my $solbasecount = 0;
my ($minsollen, $maxsollen) = (undef, undef);
foreach my $sol (@cons)
{
  $solcount++;
  my $seqlen = $sol->seqlen();
  $minsollen = $seqlen if (! defined $minsollen || $seqlen < $minsollen);
  $maxsollen = $seqlen if (! defined $maxsollen || $seqlen > $maxsollen);
  $solbasecount +=  $seqlen;
}
push @report, "Stitching Oligos:\n";
push @report, "\t$solcount oligos ($solbasecount bp)\n";
if ($solcount > 0)
{
  my $avgsollen = int ($solbasecount / $solcount);
  push @report, "\t\tAverage length: $avgsollen\n";
  push @report, "\t\tMinimum length: $minsollen\n";
  push @report, "\t\tMaximum length: $maxsollen\n";
}
push @report,  "\n\n";

# Pull all of the building blocks
#
@cons = $design->get_constructs(-kind => 'building_block');
my $bbcount = 0;
my $bbbasecount = 0;
my $bbbcount = 0;
my $bbbbasecount = 0;
my ($minbblen, $maxbblen) = (undef, undef);
my ($minbbblen, $maxbbblen) = (undef, undef);
foreach my $bb (@cons)
{
  my $seqlen = $bb->seqlen();
  my $method = $bb->method();
  if ($method eq 'magic')
  {
    $bbbcount++;
    $bbbbasecount +=  $seqlen;
    $minbbblen = $seqlen if (! defined $minbbblen || $seqlen < $minbbblen);
    $maxbbblen = $seqlen if (! defined $maxbbblen || $seqlen > $maxbbblen);
    next;
  }
  $minbblen = $seqlen if (! defined $minbblen || $seqlen < $minbblen);
  $maxbblen = $seqlen if (! defined $maxbblen || $seqlen > $maxbblen);
  $bbcount++;
  $bbbasecount +=  $seqlen;
}
push @report, "Building Blocks:\n";
push @report, "\t$bbcount building blocks ($bbbasecount bp)\n";
if ($bbcount > 0)
{
  my $avgbblen = int ($bbbasecount / $bbcount);
  push @report, "\t\tAverage length: $avgbblen\n";
  push @report, "\t\tMinimum length: $minbblen\n";
  push @report, "\t\tMaximum length: $maxbblen\n";
}
if ($bbbcount)
{
  push @report, "\t$bbbcount magic building blocks ($bbbbasecount bp)\n";
  my $avgbbblen = int ($bbbbasecount / $bbbcount);
  push @report, "\t\tAverage length: $avgbbblen\n";
  push @report, "\t\tMinimum length: $minbbblen\n";
  push @report, "\t\tMaximum length: $maxbbblen\n";
}
push @report,  "\n\n";

# Pull all of the intermediates
#
@cons = $design->get_constructs(-kind => 'intermediate');
my $intcount = 0;
my $intbasecount = 0;
my $bintcount = 0;
my $bintbasecount = 0;
my ($minintlen, $maxintlen) = (undef, undef);
my ($minbintlen, $maxbintlen) = (undef, undef);
my %iscreengroups = ();
foreach my $int (@cons)
{
  my $seqlen = $int->seqlen();
  my $method = $int->method();
  $iscreengroups{$int->id} = $int->screening_group();
  
  if ($method eq 'magic')
  {
    $bintcount++;
    $bintbasecount += $seqlen;
    $minbintlen = $seqlen if (! defined $minbintlen || $seqlen < $minbintlen);
    $maxbintlen = $seqlen if (! defined $maxbintlen || $seqlen > $maxbintlen);
    next;
  }
  $intcount++;
  $intbasecount += $seqlen;
  $minintlen = $seqlen if (! defined $minintlen || $seqlen < $minintlen);
  $maxintlen = $seqlen if (! defined $maxintlen || $seqlen > $maxintlen);
}
push @report, "Intermediates:\n";
push @report, "\t$intcount intermediates ($intbasecount bp)\n";
if ($intcount > 0)
{
  my $avgintlen = int ($intbasecount / $intcount);
  push @report, "\t\tAverage length: $avgintlen\n";
  push @report, "\t\tMinimum length: $minintlen\n";
  push @report, "\t\tMaximum length: $maxintlen\n";
}
if ($bintcount)
{
  push @report, "\t$bintcount magic intermediates ($bintbasecount bp)\n";
  my $avgbintlen = int ($bintbasecount / $bintcount);
  push @report, "\t\tAverage length: $avgbintlen\n";
  push @report, "\t\tMinimum length: $minbintlen\n";
  push @report, "\t\tMaximum length: $maxbintlen\n";
}
push @report, "\tScreening Groups:\n";
foreach my $intid (sort keys %iscreengroups)
{
  push @report, "\t$intid\t" . $iscreengroups{$intid} . "\n";
}
push @report,  "\n\n";

# Pull all of the chunks
#
@cons = $design->get_constructs(-kind => 'chunk');
my $chcount = 0;
my $chbasecount = 0;
my $bchcount = 0;
my $bchbasecount = 0;
my ($minchlen, $maxchlen) = (undef, undef);
my ($minbchlen, $maxbchlen) = (undef, undef);
foreach my $ch (@cons)
{
  my $seqlen = $ch->seqlen();
  my $method = $ch->method();
  if ($method eq 'magic')
  {
    $bchcount++;
    $bchbasecount += $seqlen;
    $minbchlen = $seqlen if (! defined $minbchlen || $seqlen < $minbchlen);
    $maxbchlen = $seqlen if (! defined $maxbchlen || $seqlen > $maxbchlen);
    next;
  }
  $chcount++;
  $chbasecount += $seqlen;
  $minchlen = $seqlen if (! defined $minchlen || $seqlen < $minchlen);
  $maxchlen = $seqlen if (! defined $maxchlen || $seqlen > $maxchlen);
}
push @report, "Chunks:\n";
push @report, "\t$chcount chunks ($chbasecount bp)\n";
if ($chcount > 0)
{
  my $avgchlen = int ($chbasecount / $chcount);
  push @report, "\t\tAverage length: $avgchlen\n";
  push @report, "\t\tMinimum length: $minchlen\n";
  push @report, "\t\tMaximum length: $maxchlen\n";
}
if ($bchcount)
{
  push @report, "\t$bchcount magic chunks ($bchbasecount bp)\n";
  my $avgbchlen = int ($bchbasecount / $bchcount);
  push @report, "\t\tAverage length: $avgbchlen\n";
  push @report, "\t\tMinimum length: $minbchlen\n";
  push @report, "\t\tMaximum length: $maxbchlen\n";
}
push @report,  "\n\n";

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

  DNAssemble_300_Summary.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

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