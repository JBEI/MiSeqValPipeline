#!/usr/bin/env perl

use Bio::GeneDesign;
use Spreadsheet::WriteExcel;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $GDV = "JGISB_Order_Stitch_$VERSION";
my %plates = (384 => {cols => 24, rows => "P"}, 96 => {cols => 12, rows => "H"});
my $mpn = 6;

local $| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT},
      'output=s'          => \$p{OUTPUT},
      'wellcount=i'       => \$p{WC},
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();

die "\n JGISB_ERROR: You must supply an input genbank file.\n"
  if (! $p{INPUT});

my ($iter, $filename, $suffix) = $GD->import_seqs($p{INPUT});

#The output path must exist, and we'll need it to end with a slash
$p{OUTPUT} = $p{OUTPUT} || ".";
$p{OUTPUT} .= "/" if (substr($p{OUTPUT}, -1, 1) !~ /[\/]/);
die "\n JGISB_ERROR: $p{OUTPUT} does not exist.\n"
  if ($p{OUTPUT} && ! -e $p{OUTPUT});

die "\n JGISB_ERROR: $p{WC} doesn't make sense.\n"
  if ($p{WC} && ! exists($plates{$p{WC}}));

################################################################################
################################# CONFIGURING ##################################
################################################################################
$p{OUTPUT} = $p{OUTPUT} || $GD->{output_dir};
$p{WC} = $p{WC} || 96;
my $onreg = qr/.+\..+\.(\d+)\./;

$p{PLATENAME} = $p{NAME} || $filename . "_stitch";
my $xlfile = $p{PLATENAME} . ".xls";
my $xlpath = $p{OUTPUT} . $xlfile;

my $order = Spreadsheet::WriteExcel->new($xlpath);
$order->compatibility_mode();

my $univF = $order->set_custom_color(15, 51, 51, 153);
my $univR = $order->set_custom_color(24, 153, 51, 0);
my %fonta = (font => 'Arial', size => 10, color => 'white', bold => 1);
my $formata = $order->add_format(%fonta, bg_color => $univF);
my $formatb = $order->add_format(%fonta, bg_color => $univR);

################################################################################
################################## ARRAYING  ###################################
################################################################################
my $platecount = 1;
my $currplate = $p{PLATENAME} . "_" . $platecount;

my $mapsheet = $order->add_worksheet("platemaps");
$mapsheet->set_column(0, $plates{$p{WC}}->{cols}-1, 20);
$mapsheet->write(0, 0, "$currplate map");
my $maprow = 1;
my $mapcol = 0;

my $fragsheet = $order->add_worksheet("fragmaps");
$fragsheet->set_column(0, $plates{$p{WC}}->{cols}-1, 20);
my $fragpcount = 1;
$fragsheet->write(0, 0, "fragment plate $fragpcount");
my $fragrow = 1;
my $fragcol = 0;
my $currfrow = 1;

my $currsheet = $order->add_worksheet($currplate);
my $xlrow = 0;
my $xlcol = 0;
my $len = 0;

my %platemaps = ($currplate => {});

my $currrow = "A";
my $currcol = 1;
my $olcount = 0;
my $currwellcount = 0;

my @cfails;
my @bfails;

my $maxbbcount = 0;
#Pull cons to assess building block counts and vectors
while ( my $obj = $iter->next_seq() )
{
  my $collec = $obj->annotation;
  my @keys = $collec->get_all_annotation_keys();
  my @comments = $collec->get_Annotations('comment');
  my $destvec = undef;
  my $origname = undef;
  foreach my $comment (@comments)
  {
    my $ctext = $comment->as_text;
    if ($ctext =~ m{destination_vector \= (.+)}ms)
    {
      $destvec = $1;
    }
    elsif ($ctext =~ m{original_name \= (.+)}ms)
    {
      $origname = $1;
    }
  }
  
  my $conid = $obj->id;
  my @bbs = grep {$_->primary_tag eq "building_block"} $obj->get_SeqFeatures;
  unless (scalar(@bbs))
  {
    warn ("JGISB_WARNING: $conid has no annotated building blocks... skipping!\n");
    push @cfails, $conid;
    next;
  }
  my $bbcount = scalar(@bbs);
  $maxbbcount = $bbcount if ($bbcount > $maxbbcount);
  my @clonvecs = map {join("", $_->get_tag_values("vector"))} @bbs;
  my $clonvec = $clonvecs[0];
  $destvec = $clonvec unless ($destvec);
  if (scalar(@bbs) == 1)
  {
    print ("$conid has only one building block... skipping!\n");
    next;
  }
  
  my $uflag = 0;
  $uflag++ if ($destvec eq $clonvec);
  my $rflag = 0;
  
  #print " working on $conid 's building blocks:\n";
  my $wellwidth = $bbcount * 2;

  my $bbnum = 0;
  my $well;
  ##check to see if it will fit in this plate
  if ($wellwidth + $currwellcount > $p{WC})
  {
    while ($currrow le $plates{$p{WC}}->{rows})
    {
      while ($currcol <= $plates{$p{WC}}->{cols})
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
#    print "Need $wellwidth and have $currwellcount; start new plate\n";
    $currrow = "A";
    $currcol = 1;
    $platecount++;
    $currplate = $p{PLATENAME} . "_" . $platecount;
    $maprow+= 3;
    $mapcol = 0;
    $fragcol = 0;
    $fragpcount++;
    $mapsheet->write($maprow-1, $mapcol, "$currplate map");
    $fragsheet->write($maprow-1, $mapcol, "fragment plate $fragpcount");
    $currsheet = $order->add_worksheet($currplate);
    $xlrow = 0;
    $xlcol = 0;
    $currwellcount = 0;
  }
  foreach my $bb (@bbs)
  {
    my $bbname = join("", $bb->get_tag_values("name"));
    $rflag++ if ($bb->has_tag("subset"));
    $bbnum++;
    my $blank = 0;
  
    ##FORWARD
    if ($bb->has_tag('vector') && $bbnum != 1)
    {
      my $lname = $bbname . "_sF";
      my $lseq = join("", $bb->get_tag_values('stitch_left'));
      $lseq =~ s/\s//g;
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol+0, $well);
      $currsheet->write($xlrow, $xlcol+1, $lname);
      $currsheet->write($xlrow, $xlcol+2, $lseq);
      #print "$bbname $lname in $well\n";
      $mapsheet->write($maprow, $mapcol, $lname);
      $len += length($lseq);
      $olcount++;
    }
    else
    {
      $blank++;
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol+0, $well);
      $currsheet->write($xlrow, $xlcol+1, "empty");
      $currsheet->write($xlrow, $xlcol+2, "empty");
      $fragsheet->write($maprow, $fragcol, $bbname, $formata);
      if ($uflag && $bbnum == 1 )
      {
        #print "$bbname UF in $well\n";
        $mapsheet->write($maprow, $mapcol, "UF ($destvec)", $formata);
      }
      if ($rflag)
      {
        $mapsheet->write($maprow, $mapcol, "UF ($clonvec)", $formata);
      }
    }
    $xlrow++;
    $currcol++;
    $mapcol++;
    $currwellcount++;
  
  
    ##REVERSE
    if ($bb->has_tag('vector') && $bbnum != $bbcount)
    {
      my $rname = $bbname . "_sR";
      my $rseq = join("", $bb->get_tag_values("stitch_right"));
      $rseq =~ s/\s//g;
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol+0, $well);
      $currsheet->write($xlrow, $xlcol+1, $rname);
      $currsheet->write($xlrow, $xlcol+2, $rseq);
      #print "$bbname $rname in $well\n";
      $mapsheet->write($maprow, $mapcol, $rname);
      $len += length($rseq);
      $olcount++;
    }
    else
    {
      $blank++;
      $well = wellname($currrow, $currcol);
      $currsheet->write($xlrow, $xlcol+0, $well);
      $currsheet->write($xlrow, $xlcol+1, "empty");
      $currsheet->write($xlrow, $xlcol+2, "empty");
      $fragsheet->write($maprow, $fragcol, $bbname, $formatb);
      if ($uflag && $bbnum == $bbcount)
      {
        #print "$bbname UR in $well\n";
        $mapsheet->write($maprow, $mapcol, "UR ($destvec)", $formatb);
      }
      if ($rflag)
      {
        $mapsheet->write($maprow, $mapcol, "UR ($clonvec)", $formatb);
      }
    }
    if ($blank == 0)
    {
      $fragsheet->write($maprow, $fragcol, $bbname);
    }
    $xlrow++;
    $currcol++;
    $fragcol++;
    $mapcol++;
    $currwellcount++;
  
    if ($currcol > $plates{$p{WC}}->{cols})
    {
      $currcol = 1;
      $currrow++;
      $maprow++;
      $fragcol = 0;
      $mapcol = 0;
    }
  }
}

while ($currrow le $plates{$p{WC}}->{rows})
{
  while ($currcol <= $plates{$p{WC}}->{cols})
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
print "Wrote $xlpath\n\n";

print "\nGenerated $platecount $p{WC} well plates containing $olcount oligos ";
print "and $len bases.\n\n";

print "JGISB_WARNING: Constructs (@cfails) were not processed.\n\n" if (scalar @cfails);
print "JGISB_WARNING: Building blocks (@bfails) were not ordered.\n\n" if (scalar @bfails);

print "brought to you by $GDV\n\n";

exit;

sub wellname
{
  my ($row, $num) = @_;
  $num = "0" . $num while length($num) < 2;
  return $row . $num;
}


__END__

=head1 NAME

  JGISB_Order_Stitch.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  Converts a genbank file containing stitching oligos to an excel order sheet
  
=head1 USAGE

=head1 ARGUMENTS

Required arguments:
  
  -i,   --input : A genbank file to use as an input

Optional arguments:

  -o,   --output: Where should files be written
  -w,   --wellcount : The number of wells in a plate (default 96)
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