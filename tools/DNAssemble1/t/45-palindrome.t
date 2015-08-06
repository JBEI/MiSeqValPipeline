#! /usr/bin/perl

#FAILS 25 IF perl -T

use Test::More;
use Bio::DNAssemble;
use Bio::Seq;

use strict;
use warnings;

my $DNA = Bio::DNAssemble->new();

eval {require Bio::Factory::EMBOSS};
if($@)
{
    plan skip_all => 'Bio::Factory::EMBOSS not installed';
}
elsif (! $DNA->EMBOSS)
{
    plan skip_all => 'EMBOSS support not installed';
}
else
{
    plan tests => 1;
}

my $pal = Bio::Seq->new(
  -seq => "ATAACTTCGTATAATGTACATTATACGAAGTTAT",
  -id => "loxPsym"
);
my $nopal = Bio::Seq->new(
  -seq => "ACATTATACGAAGTTAT",
    -id => "halfloxPsym"
);

my $rarr = [$nopal];

my $tarr = $DNA->filter_palindromes(-sequences => [$pal, $nopal]);

is_deeply ($tarr, $rarr, "filter palindromes");
                                         

