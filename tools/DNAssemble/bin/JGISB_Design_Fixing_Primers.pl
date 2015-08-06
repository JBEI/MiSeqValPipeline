#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use English qw(-no_match_vars);
use Carp;
use CGI qw(:standard);
use Bio::GeneDesign;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_Design_Fixing_Primers_$VERSION";

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'    => \$p{HELP},
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $GD = Bio::GeneDesign->new();
my $raw = undef;
while (<>)
{
  $raw .= $_;
}

if ($raw !~ m{\#} || ! defined $raw)
{
  die "\n JGISB_ERROR: Can't understand input!\n";
}

my @arr = split m{\#}x, $raw;
my $reference = $arr[1];
$reference =~ s{[\s\n]}{}xg;
my $reflen = length $reference;

#Pull variants and determine range
my @variants = split m{\n}x, $arr[0];
my ($min, $max) = ($reflen, 0);
foreach my $variant (@variants)
{
  my @atts = split m{\s}x, $variant;
  my $offset = $atts[1];
  if ($offset < $min)
  {
    $min = $offset;
  }
  if ($offset > $max)
  {
    $max = $offset;
  }
}
if (abs($reflen - $min) < $reflen / 10)
{
  print "-1\n";
  exit;
}
elsif (abs($reflen - $max) < $reflen / 10)
{
  print "-1\n";
  exit;
}
elsif (abs($max - $min) > 5)
{
  print "-1\n";
  exit;
}

my ($bit, $upol, $dnol) = fixing_primers($min, $max, $reference);
print "$bit\t$upol\t$dnol\n";


exit;

################################################################################
################################# SUBROUTINES ##################################
################################################################################

sub fixing_primers
{
  my ($fstart, $fstop, $refseq) = @_;
  my $bit = 1;
  my $uplen = 10;
  my $upstart = ($fstart - 1) - $uplen;
  my $upol = substr $refseq, $upstart, $uplen;
  my $upmelt = $GD->melt($upol);
  while ($upmelt < 58)
  {
    $upstart = $upstart - 1;
    $uplen++;
    $upol = substr $refseq, $upstart, $uplen;
    $upmelt = $GD->melt($upol);
  }
  while ($upmelt > 60.5)
  {
    $upstart++;
    $uplen = $uplen - 1;
    $upol = substr $refseq, $upstart, $uplen;
    $upmelt = $GD->melt($upol);
  }

  my $upcount = $GD->count($upol);
  $bit++ if $upcount->{GCp} > 58;
  $upol .= lc substr $refseq, $upstart + $uplen, 10;
  $upol = $GD->complement(-sequence => $upol, -reverse => 1);
  my $dnstart = $fstop;
  my $dnlen = 10;
  my $dnol = substr $refseq, $dnstart, $dnlen;
  my $dnmelt = $GD->melt($dnol);
  while ($dnmelt < 58)
  {
    $dnlen++;
    $dnol = substr $refseq, $dnstart, $dnlen;
    $dnmelt = $GD->melt($dnol);
  }
  while ($dnmelt > 60.5)
  {
    $dnlen = $dnlen - 1;
    $dnol = substr $refseq, $dnstart, $dnlen;
    $dnmelt = $GD->melt($dnol);
  }

  my $dncount = $GD->count($dnol);
  $bit++ if $dncount->{GCp} > 58;
  my $homu = $GD->contains_homopolymer(-sequence => $upol);
  my $homd = $GD->contains_homopolymer(-sequence => $dnol);
  $bit++ if ($homu || $homd);
  my $morednol = lc substr $refseq, $dnstart - 10, 10;
  $dnol = $morednol . $dnol;

  return ($bit, $upol, $dnol);
}

__END__

=head1 NAME

  JGISB_Design_Fixing_Primers.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

  Given a reference sequence and a list of variants, determine if one set of
  fixing oligos can be made to remove the variants

=head1 USAGE

  You can either pipe input in
    cat flyfix.txt | JGISB_Design_Fixing_Primers.pl

  or provide a file as an argument
    JGISB_Design_Fixing_Primers.pl flyfix.txt


  The input should resemble that from a vcf file for the top half. Then a hash,
  then the reference sequence (the right sequence)

    GH48_129	1525	.	T	IG	1.00	0	NS=1;DP=1
    GH48_129	1528	.	G	IG	1.00	0	NS=1;DP=1
    GH48_129	1530	.	A	IT	1.00	0	NS=1;DP=1
    #
    ATACACATACCCGGGATTTAGGTGACACTATAGAATACACGGAATTCGCGTTTTTATTTT
    TAATTTTCTTTCAAATACTTCTAGCTAGAGTATTTTTACAACAATTACCAACAACAACAA
    CAAACAACAACAACATTACATTTTACATTCTACAACTACAGCCGCGATCGCCATGGCAGA
    ...


  But note that the third through nth columns in the variant list are not
  necessary! Matt, Angie, you can generate fixes by faking the lines. Like so
  for a fix at base 1320:

    UGHIhatethis  1320
    #
    ATACACATACCCGGGATTTAGGTGACACTATAGAATACACGGAATTCGCGTTTTTATTTT
    TAATTTTCTTTCAAATACTTCTAGCTAGAGTATTTTTACAACAATTACCAACAACAACAA
    CAAACAACAACAACATTACATTTTACATTCTACAACTACAGCCGCGATCGCCATGGCAGA
    ...

  As long as the second column (tab delimited) is an offset inside the reference
  sequence, fixing primers should be generated.

  A -1 return means that no primers could be made.
  Otherwise, you get a number, primer1, primer2.
  The higher the number, the crappier the primers (1 is best)

=head1 ARGUMENTS

Required arguments:

  A variant/reference input file.

Optional arguments:

  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

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