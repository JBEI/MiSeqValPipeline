#! /usr/bin/perl -T

use Test::More tests => 1;

use Bio::DNAssemble;
use Bio::Seq;

use strict;
use warnings;

my $DNA = Bio::DNAssemble->new();

my $seq = "ATGGACAGATCTTGGAAGCAGAAGCTGAACCGCGACACCGTGAAGCTGACCGAGGTGATGACCTGGA";
my $busted = "ATGGACAGATCTTGGAAXGCAGAAGCTGAAvCCGCGACACCGTGA AGCTGACCGAGGCCTGGA";
my $fixed = "ATGGACAGATCTTGGAAGCAGAAGCTGAACCGCGACACCGTGAAGCTGACCGAGGCCTGGA";

#TESTING purify
subtest "purify" => sub
{
  plan tests => 4;
  
  my ($tnotambig, $notarr) = $DNA->purify($seq);
  my ($tisambig, $yesarr) = $DNA->purify($busted);
  my $rarr = [' ', 'X', 'v'];
  is_deeply($yesarr, $rarr, "sequence was not pure");
  is_deeply($notarr, [], "sequence was pure");
  is($tnotambig, $seq, "sequence remains the same");
  is($tisambig, $fixed, "sequence repaired");
};
