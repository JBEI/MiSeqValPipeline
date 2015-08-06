#!/usr/bin/env perl

use Bio::GeneDesign;
use Getopt::Long;
use Pod::Usage;
use Spreadsheet::WriteExcel;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_Primer_forEnzfromGB_$VERSION";

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT},
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();

#The input file must exist and be a format we care to read.
die "\n JGISB_ERROR: You must supply an input file.\n"
  if (! $p{INPUT});

my ($iter, $filename, $suffix) = $GD->import_seqs($p{INPUT});


################################################################################
################################# CONFIGURING ##################################
################################################################################

$p{PLATENAME} = 'Batch30' . "_IntStitch";
my $xlfile = $p{OUTPUT} || $p{PLATENAME} . ".xls";
my $xlpath = $xlfile;

my %plates = (384 => {cols => 24, rows => "P", rownum => 16}, 96 => {cols => 12, rows => "H", rownum => 8});
my $WC = 96;

my $stemp = 65;
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
my $mapfix = 1;
my $mapcol = 0;

my $fragsheet = $order->add_worksheet('fragmaps');
$fragsheet->set_column(0, $plates{$WC}->{cols}-1, 20);
my $fragpcount = 1;
$fragsheet->write(0, 0, "fragment plate $fragpcount");
my $fragcol = 0;

my $ordersheet = $order->add_worksheet($currplate);
my $orderrow = 0;
my $ordercol = 0;
my $len = 0;

my %platemaps = ($currplate => {});

my $currrow = 'A';
my $currcol = 1;
my $olcount = 0;
my $currwellcount = 0;

################################################################################
################################## ARRAYING  ###################################
################################################################################
my %chhsh;
while ( my $obj = $iter->next_seq() )
{
  my $chid = $obj->id;
  my @bbs = grep {$_->primary_tag eq 'building_block'} $obj->get_SeqFeatures();
  my %sethsh;
  foreach my $bb (@bbs)
  {
    my $bbid = join q{}, $bb->get_tag_values('name');
    my $setn = join q{}, $bb->get_tag_values('subset');
    my $seq  = $obj->subseq($bb->start, $bb->end);
    if (exists $sethsh{$setn})
    {
      push @{$sethsh{$setn}}, $seq;
    }
    else
    {
      $sethsh{$setn} = [$seq];
    }
  }
  $chhsh{$chid} = \%sethsh;
}

my $well = undef;
foreach my $chid (sort keys %chhsh)
{
  my %sethsh = %{$chhsh{$chid}};
  my @setlist = sort keys %sethsh;
  my $setnum = scalar @setlist;
  my $setreg = '0';
  if ((2 * $setnum) + $currwellcount > $WC)
  {
    while ($currrow le $plates{$WC}->{rows})
    {
      $well = wellname($currrow, $currcol);
      $ordersheet->write($orderrow, $ordercol + 0, $well);
      $ordersheet->write($orderrow, $ordercol + 1, 'empty');
      $ordersheet->write($orderrow, $ordercol + 2, 'empty');
      $currcol++;
      $orderrow++;
      $well = wellname($currrow, $currcol);
      $ordersheet->write($orderrow, $ordercol + 0, $well);
      $ordersheet->write($orderrow, $ordercol + 1, 'empty');
      $ordersheet->write($orderrow, $ordercol + 2, 'empty');
      $orderrow++;
      $currrow++;
      $currcol--;
    }
    $currrow = "A";
    $currcol = 1;
    $platecount++;
    $currplate = $p{PLATENAME} . q{_} . $platecount;
    $maprow  = $mapfix + 12;
    $mapfix = $maprow;
    $mapcol = 0;
    $fragcol = 0;
    $fragpcount++;
    $mapsheet->write($mapfix - 1, $mapcol, "$currplate map");
    $fragsheet->write($mapfix - 1, $mapcol, "fragment plate $fragpcount");
    $ordersheet = $order->add_worksheet($currplate);
    $orderrow = 0;
    $ordercol = 0;
    $currwellcount = 0;
  }
  foreach my $set (@setlist)
  {
    $setreg++;
    my @bbarr = @{$sethsh{$set}};
    my ($lprimer, $exrprimer) = $GD->make_amplification_primers(
        -sequence    => $bbarr[0],
        -temperature => $stemp,
    );
    my ($exlprimer, $rprimer) = $GD->make_amplification_primers(
        -sequence    => $bbarr[-1],
        -temperature => $stemp,
    );
    my $intid = $chid . q{_} . $set;
    if ($setreg == 1)
    {
      $lprimer = undef;
    }
    elsif ($setreg == $setnum)
    {
      $rprimer = undef;
    }
    
    my ($lname, $rname) = ($intid . q{.F}, $intid . q{.R});
    if ($lprimer)
    {
      $well = wellname($currrow, $currcol);
      $ordersheet->write($orderrow, $ordercol + 0, $well);
      $ordersheet->write($orderrow, $ordercol + 1, $lname);
      $ordersheet->write($orderrow, $ordercol + 2, $lprimer);
      $mapsheet->write($maprow, $mapcol, $lname);
      $len += length $lprimer;
      $olcount++;
    }
    else
    {
      $well = wellname($currrow, $currcol);
      $ordersheet->write($orderrow, $ordercol + 0, $well);
      $ordersheet->write($orderrow, $ordercol + 1, "empty");
      $ordersheet->write($orderrow, $ordercol + 2, "empty");
      $mapsheet->write($maprow, $mapcol, "UL", $formata);
    }
    $orderrow++;
    $mapcol++;
    $currcol++;
    $currwellcount++;
    if ($rprimer)
    {
      $well = wellname($currrow, $currcol);
      $ordersheet->write($orderrow, $ordercol + 0, $well);
      $ordersheet->write($orderrow, $ordercol + 1, $rname);
      $ordersheet->write($orderrow, $ordercol + 2, $rprimer);
      $mapsheet->write($maprow, $mapcol, $rname);
      $len += length $rprimer;
      $olcount++;
    }
    else
    {
      $well = wellname($currrow, $currcol);
      $ordersheet->write($orderrow, $ordercol + 0, $well);
      $ordersheet->write($orderrow, $ordercol + 1, "empty");
      $ordersheet->write($orderrow, $ordercol + 2, "empty");
      $mapsheet->write($maprow, $mapcol, "UR", $formatb);
    }
    $fragsheet->write($maprow, $fragcol, $intid);
    $orderrow++;
    $currrow++;
    $currcol--;
    $maprow++;
    $mapcol--;
    $currwellcount++;
  
    if ($currrow gt $plates{$WC}->{rows})
    {
      $currcol += 2;
      $currrow = 'A';
      $maprow = $mapfix;
      $mapcol += 2;
      $fragcol ++;
    }
  }
}
while ($currrow le $plates{$WC}->{rows})
{
  $well = wellname($currrow, $currcol);
  $ordersheet->write($orderrow, $ordercol + 0, $well);
  $ordersheet->write($orderrow, $ordercol + 1, 'empty');
  $ordersheet->write($orderrow, $ordercol + 2, 'empty');
  $currcol++;
  $orderrow++;
  $well = wellname($currrow, $currcol);
  $ordersheet->write($orderrow, $ordercol + 0, $well);
  $ordersheet->write($orderrow, $ordercol + 1, 'empty');
  $ordersheet->write($orderrow, $ordercol + 2, 'empty');
  $orderrow++;
  $currrow++;
  $currcol--;
}

$order->close();

print  "\n\n";
print  "Wrote $xlpath\n\n";
print  "Generated $platecount $WC well plates containing $olcount ";
print  "oligos and $len bases.\n\n";

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

JGISB_Primer_forEnzfromGB.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 USAGE

=head1 ARGUMENTS

Required arguments:

  -i,   --input : A genbank file to use as an input

Optional arguments:

  -o,   --output: Where should files be written
  -s,   --skip_adds: Do not add suffixes or prefixes to construct sequences
  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012, JGI Syn Bio developers
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