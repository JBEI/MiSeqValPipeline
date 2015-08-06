#!/usr/bin/env perl

use Bio::DNAssemble;
use English qw( -no_match_vars );
use Spreadsheet::WriteExcel;
use Getopt::Long;
use Pod::Usage;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_320_Stitchers_' . $VERSION;
my $DNAS = '_320';

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

$p{PLATENAME} = $filename . "_StitchingOligos";
my $xlfile = $p{OUTPUT} || $p{PLATENAME} . ".xls";
my $xlpath = $xlfile;


################################################################################
################################# CONFIGURING ##################################
################################################################################
my %plates = (384 => {cols => 24, rows => "P"}, 96 => {cols => 12, rows => "H"});
my $WC = 96;

my $order = Spreadsheet::WriteExcel->new($xlpath);
$order->compatibility_mode();

my $univF = $order->set_custom_color(15, 51, 51, 153);
my $univR = $order->set_custom_color(24, 153, 51, 0);
my %fonta = (font => 'Arial', size => 10, color => 'white', bold => 1);
my $formata = $order->add_format(%fonta, bg_color => $univF);
my $formatb = $order->add_format(%fonta, bg_color => $univR);

my $platecount = 1;
my $currplate = $p{PLATENAME} . q{_} . $platecount;
my $mapsheet = $order->add_worksheet('platemaps');
$mapsheet->set_column(0, $plates{$WC}->{cols}-1, 20);
$mapsheet->write(0, 0, "$currplate map");
my $maprow = 1;
my $mapcol = 0;

my $fragsheet = $order->add_worksheet('fragmaps');
$fragsheet->set_column(0, $plates{$WC}->{cols}-1, 20);
my $fragpcount = 1;
$fragsheet->write(0, 0, "fragment plate $fragpcount");
my $fragcol = 0;

my $sinsheet = $order->add_worksheet('singlets');
$sinsheet->write(0, 0, 'intermediate name');
$sinsheet->write(0, 1, 'fragment');
my $sinrow = 1;

my $intsheet = $order->add_worksheet('products');
$intsheet->write(0, 0, 'intermediate name');
$intsheet->write(0, 1, 'fragments');
my $introw = 1;

my $currsheet = $order->add_worksheet($currplate);
my $xlrow = 0;
my $xlcol = 0;
my $len = 0;

my %platemaps = ($currplate => {});

my $currrow = "A";
my $currcol = 1;
my $olcount = 0;
my $currwellcount = 0;

################################################################################
################################# EXTRACTING  ##################################
################################################################################
my @bbs = $design->get_constructs(-kind => 'building_block');
my %bbhsh = map {$_->id => $_} @bbs;
my %seenbbs = ();
my @ints = $design->get_constructs(-kind => 'intermediate');
foreach my $int (@ints)
{
  my $intcol = 1;
  my $intname = $int->id();
  my @pools = $int->get_pools(-method => 'fragamp');
  my @subconstructs;
  foreach my $pool (@pools)
  {
    my $list = $pool->subconstruct();
    my @sublist = split q{; }, $list;
    push @subconstructs, @sublist;
  }
  my @bbs = grep {exists $bbhsh{$_}} @subconstructs;
  my $bbcount = scalar @bbs;
  if ($bbcount == 1)
  {
    $sinsheet->write($sinrow, 0, $intname);
    $sinsheet->write($sinrow, 1, $bbs[0]);
    $sinrow++;
    next;
  }
  my @stitchers = grep {! exists $bbhsh{$_}} @subconstructs;
  my $thiscount = scalar @stitchers;
  
  ##check to see if it will fit in this plate
  if ($thiscount + $currwellcount > $WC)
  {
    while ($currrow le $plates{$WC}->{rows})
    {
      while ($currcol <= $plates{$WC}->{cols})
      {
        my $well = wellname($currrow, $currcol);
        $platemaps{$currplate}->{$well} = q{.};
        $currsheet->write($xlrow, $xlcol + 0, $well);
        $currsheet->write($xlrow, $xlcol + 1, 'empty');
        $currsheet->write($xlrow, $xlcol + 2, 'empty');
        $xlrow++;
        $currcol++;
      }
      $currrow++;
      $currcol = 1;
    }
#    print "Need $wellwidth and have $currwellcount; start new plate\n";
    $currrow = "A";
    $currcol = 1;
    $platecount++;
    $currplate = $p{PLATENAME} . q{_} . $platecount;
    $maprow+= 3;
    $mapcol = 0;
    $fragcol = 0;
    $fragpcount++;
    $mapsheet->write($maprow - 1, $mapcol, "$currplate map");
    $fragsheet->write($maprow - 1, $mapcol, "fragment plate $fragpcount");
    $currsheet = $order->add_worksheet($currplate);
    $xlrow = 0;
    $xlcol = 0;
    $currwellcount = 0;
  }
  
  my $passint = 1;
  my $well = undef;
  $intsheet->write($introw, 0, $intname);
  while (scalar @stitchers)
  {
    my $leftol = shift @stitchers;
    $leftol = undef if ($passint == 1);
    my $rightol = shift @stitchers;
    $rightol = undef if ($passint == $bbcount);
    my $bb = shift @bbs;
    if (exists $seenbbs{$bb})
    {
      $intsheet->write($introw, $intcol, $bb, $formata);
      $intcol++;
      next;
    }
    $intsheet->write($introw, $intcol, $bb);
    $intcol++;
    $seenbbs{$bb}++;
    $fragsheet->write($maprow, $fragcol, $bb, $formatb);
    
    if ($leftol)
    {
      my @leftcons = $design->get_constructs(-id => $leftol);
      my $leftcon = $leftcons[0];
      my $leftseq = $leftcon->sequence();
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol + 0, $well);
      $currsheet->write($xlrow, $xlcol + 1, $leftol);
      $currsheet->write($xlrow, $xlcol + 2, $leftseq);
      #print "$bb $leftol in $well\n";
      $mapsheet->write($maprow, $mapcol, $leftol);
      $len += length $leftseq;
      $olcount++;
    }
    else
    {
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol + 0, $well);
      $currsheet->write($xlrow, $xlcol + 1, "empty");
      $currsheet->write($xlrow, $xlcol + 2, "empty");
      $mapsheet->write($maprow, $mapcol, "UL", $formata);
    }
    $xlrow++;
    $currcol++;
    $mapcol++;
    $currwellcount++;
    if ($rightol)
    {
      my @rightcons = $design->get_constructs(-id => $rightol);
      my $rightcon = $rightcons[0];
      my $rightseq = $rightcon->sequence();
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol + 0, $well);
      $currsheet->write($xlrow, $xlcol + 1, $rightol);
      $currsheet->write($xlrow, $xlcol + 2, $rightseq);
      #print "$bb $rightol in $well\n";
      $mapsheet->write($maprow, $mapcol, $rightol);
      $len += length $rightseq;
      $olcount++;
    }
    else
    {
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol + 0, $well);
      $currsheet->write($xlrow, $xlcol + 1, "empty");
      $currsheet->write($xlrow, $xlcol + 2, "empty");
      $mapsheet->write($maprow, $mapcol, "UR", $formatb);
    }
    $xlrow++;
    $currcol++;
    $fragcol++;
    $mapcol++;
    $currwellcount++;
    $passint++;
  
    if ($currcol > $plates{$WC}->{cols})
    {
      $currcol = 1;
      $currrow++;
      $maprow++;
      $fragcol = 0;
      $mapcol = 0;
    }
  }
  $introw++;
}


################################################################################
################################## REPORTING ###################################
################################################################################
while ($currrow le $plates{$WC}->{rows})
{
  while ($currcol <= $plates{$WC}->{cols})
  {
    my $well = wellname($currrow, $currcol);
    $platemaps{$currplate}->{$well} = ".";
    $currsheet->write($xlrow, $xlcol+0, $well);
    $currsheet->write($xlrow, $xlcol+1, "empty");
    $currsheet->write($xlrow, $xlcol+2, "empty");
    $xlrow++;
    $currcol++;
  }
  $currrow++;
  $currcol = 1;
}
$order->close();

print {$OUTH} "\n\n";
print {$OUTH} "Wrote $xlpath\n\n";
print {$OUTH} "Generated $platecount $WC well plates containing $olcount ";
print {$OUTH} "oligos and $len bases.\n\n";
print {$OUTH} $DNA->attitude() . " brought to you by $DNAV\n\n";

print {$OUTH} "\n\n******* $DNAV FINISHED\n\n";
close $OUTH;
exit;


################################################################################
################################# SUBROUTINES ##################################
################################################################################
sub wellname
{
  my ($row, $num) = @_;
  $num = "0" . $num while length($num) < 2;
  return $row . $num;
}

__END__


=head1 NAME

  DNAssemble_320_Stitchers.pl

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