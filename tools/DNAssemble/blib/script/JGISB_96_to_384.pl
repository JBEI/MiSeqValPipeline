#!/usr/bin/env perl

use IPC::Open2;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_96_to_384_$VERSION";

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'number=s'          => \$p{NUMBER},
);
pod2usage(-verbose=>99) if ($p{HELP});

my %map = (1 => 'A', 2 => 'B', 3 => 'C', 4 => 'D', 5 => 'E', 6 => 'F', 7 => 'G',
  8 => 'H', 9 => 'I', 10 => 'J', 11 => 'K', 12 => 'L', 13 => 'M', 14 => 'N',
  15 => 'O', 16 => 'P');

################################################################################
################################ SANITY  CHECK #################################
################################################################################
die "\n JGISB_ERROR: You must supply an input number greater than 1.\n"
  if (! $p{NUMBER} || $p{NUMBER} <= 1);


################################################################################
################################### WRITING ####################################
################################################################################
my ($col, $row) = (1, 1);
my $wellcount = 0;
my $firstwell = undef;
my $lastwell = undef;
my ($rcol, $rrow) = (1, 1);
my $rwell = $map{$rrow} . $rcol;
my $rcount = 1;
my $dpcount = 1;
my @flags = (12, 24, 36);
open (my $REPORT, '>', "plate$dpcount.csv") || die "Can't write: $OS_ERROR\n";
foreach my $number (1 .. $p{NUMBER})
{
  #print "loading $number\n";
  my $y = $row;
  while ($y < $row + 16)
  {
    my $well = $map{$y} . $col;
    print $REPORT "$rwell,$well,1800\n";
    $firstwell = $well if (! defined $firstwell);
    $lastwell = $well;
    $wellcount++;
    $y += 2;
    if ($wellcount % 22 == 0)
    {
      $rcol++;
      if ($rcol > 24)
      {
        $rcol = 1;
        $rrow++;
      }
      $rwell = $map{$rrow} . $rcol;
      $rcount++;
    }
  }
  if ($number == $flags[0])
  {
    $col = 2;
    $row = 1;
  }
  elsif ($number == $flags[1])
  {
    $col = 1;
    $row = 2;
  }
  elsif ($number == $flags[2])
  {
    $col = 2;
    $row = 2;
  }
  else
  {
    $col += 2;
  }
  if ($number % 48 == 0)
  {
    $dpcount++;
    close $REPORT;
    open ($REPORT, '>plate' . $dpcount . '.csv')
      || die "Can't write: $OS_ERROR\n";
    $col = 1;
    $row = 1;
    @flags = map {$_ += 48} @flags;
  }
}

close $REPORT;
print "$wellcount from $firstwell to $lastwell\n";
print "$rcount reageant wells from $firstwell to $rwell\n";

print "formula says: ", sprintf "%.1f", $p{NUMBER} * 8 / 22, " reagent wells\n";

print "\n\n";
print "brought to you by $JGISBV\n\n";

exit;


__END__

=head1 NAME

  JGISB_96_to_384.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 USAGE


=head1 ARGUMENTS

Required arguments:

  -i,   --input :

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