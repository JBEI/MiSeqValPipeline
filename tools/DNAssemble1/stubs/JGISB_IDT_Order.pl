#!/usr/bin/env perl

use Bio::GeneDesign;
use Spreadsheet::WriteExcel;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $JGISBV = "JGISB_IDT_Order_$VERSION";
my %TEMPL = (384 => {cols => 24, rows => "P"}, 96 => {cols => 12, rows => "H"});
my $mpn = 6;

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT},
      'output=s'          => \$p{OUTPUT},
      'batch=s'            => \$p{BATCH},
      'wellcount=i'       => \$p{WC},
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

#The output path must exist, and we'll need it to end with a slash
$p{OUTPUT} = $p{OUTPUT} || ".";
$p{OUTPUT} .= "/" if (substr($p{OUTPUT}, -1, 1) !~ /[\/]/);
die "\n GDERROR: $p{OUTPUT} does not exist.\n"
  if ($p{OUTPUT} && ! -e $p{OUTPUT});

$p{BATCH} = $p{BATCH} || '1';
my $batch = $GD->pad($p{BATCH}, 3);
  
#The wellcount should be 96 or 384.
$p{WC} = $p{WC} || 96;
die "\n JGISB_ERROR: $p{WC} doesn't make sense.\n"
  if ($p{WC} && ! exists($TEMPL{$p{WC}}));

################################################################################
################################# CONFIGURING ##################################
################################################################################
$p{PLATENAME} = $p{NAME} || $filename . "_plate";
my $xlfile = $p{PLATENAME} . ".xls";
my $xlpath = $p{OUTPUT} . $xlfile;

my $order = Spreadsheet::WriteExcel->new($xlpath);
$order->compatibility_mode();

my $header = $order->add_format();
$header->set_bold();
$header->set_bg_color('yellow');

################################################################################
################################## ARRAYING  ###################################
################################################################################
my $xlrow = 1;
my $len = 0;
my $currrow = "A";
my $currcol = 1;
my $olcount = 0;
my $wellcount = 0;
my @cfails;
my @bfails;
my $plates = [];
my $currplate = {};

while ( my $obj = $iter->next_seq() )
{
  my $conid = $obj->id;
  my @seqfeats = $obj->get_SeqFeatures;
  my @bbs = grep {$_->primary_tag eq "building_block"} @seqfeats;
  unless (scalar(@bbs))
  {
    warn ("JGISB_WARNING: $conid has no annotated building blocks... skipping!\n");
    push @cfails, $conid;
    next;
  }
  my @ols = grep {$_->primary_tag eq "assembly_oligo"} @seqfeats;
  unless (scalar(@ols))
  {
    warn ("JGISB_WARNING: $conid has no annotated assembly oligos... skipping!\n");
    push @cfails, $conid;
    next;
  }
  foreach my $bb (@bbs)
  {
    my $bbname = join("", $bb->get_tag_values("name"));
    my $tbbname = quotemeta $bbname;
    unless ($bb->has_tag("pca_pool_count"))
    {
      warn ("JGISB_WARNING: $bbname in $conid has no pca pool count - skipping\n");
      push @bfails, $bbname;
      next;
    }
    my $poolnum = join("", $bb->get_tag_values("pca_pool_count"));
    my @oligos = grep {join("", $_->get_tag_values("name")) =~ /^$tbbname/} @ols;
    foreach my $pool (1.. $poolnum)
    {
      my @pools = grep {join("", $_->get_tag_values("name")) =~ /^$tbbname\.($pool)\./} @oligos;
      foreach my $oligo (@pools)
      {
        my $olname = join("", $oligo->get_tag_values("name"));
        my $olseq = join("", $oligo->get_tag_values("sequence"));
        $olseq =~ s/\s//g;
        my $well = wellname($currrow, $currcol);
        $currplate->{$xlrow} = [$well, $olname, $olseq];
        $len += length($olseq);
        $xlrow++;
        $currcol++;
        $olcount++;
      }
      while ($currcol % ($mpn + 1) != 0)
      {
        $xlrow++;
        $currcol++;
        last if ($currcol > $TEMPL{$p{WC}}->{cols});
      }
    }
    while ($currcol <= $TEMPL{$p{WC}}->{cols})
    {
      $xlrow++;
      $currcol++;
    }
    $currcol = 1;
    $currrow++;
    if ($currrow gt $TEMPL{$p{WC}}->{rows})
    {
      my $newplate = $currplate;
      push @$plates, $newplate;
      $currplate = {};
      $xlrow = 1;
      $currrow = "A";
      $currcol = 1;
      $wellcount = 0;
    }
  }
}
push @$plates, $currplate;
my $platecount = scalar @$plates;
my $totplate = $GD->pad($platecount, 2);
my $platenum = $GD->pad('1', 2);

foreach my $plate (@$plates)
{
  my $barcode = $batch . $platenum . $totplate;
  my $currsheet = $order->add_worksheet($barcode);
  $currsheet->keep_leading_zeros();
  my $pcol = "A";
  my $prow = '1';
  $currsheet->write(0, 0, "Well", $header);
  $currsheet->write(0, 1, "Name", $header);
  $currsheet->write(0, 2, "Sequence", $header);
  
  for (my $xrow = 1; $xrow <= $p{WC}; $xrow++)
  {
    my $frow = $GD->pad($prow, 2);
    $currsheet->write($xrow, 0, $pcol . $frow);
    if (exists $plate->{$xrow})
    {
      my ($well, $olname, $olseq) = @{$plate->{$xrow}};
      print "woops $olname! $pcol$frow $well\n" if ($well ne $pcol . $frow);
      $currsheet->write($xrow, 1, $olname);
      $currsheet->write($xrow, 2, $olseq);
    }
    else
    {
      $currsheet->write($xrow, 1, "empty");
      $currsheet->write($xrow, 2, "empty");
    }
    $prow++;
    if ($prow > $TEMPL{$p{WC}}->{cols})
    {
      $prow = '1';
      $pcol++;
    }
  }
  $platenum = $GD->pad($platenum+1, 2);
}

$order->close();
print "\n\n";

print "Wrote $xlpath\n";

print "\nGenerated $platecount $p{WC} well plates containing $olcount oligos ";
print "and $len bases.\n\n";

print "JGISB_WARNING: Constructs (@cfails) were not ordered.\n\n" if (scalar @cfails);
print "JGISB_WARNING: Building blocks (@bfails) were not ordered.\n\n" if (scalar @bfails);

print $GD->attitude() . " brought to you by $JGISBV\n\n";

exit;

sub wellname
{
  my ($row, $num) = @_;
  $num = "0" . $num while length($num) < 2;
  return $row . $num;
}

__END__

=head1 NAME

  JGISB_IDT_Order.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  Converts a genbank file containing assembly oligos to an excel order sheet

=head1 USAGE

=head1 ARGUMENTS

Required arguments:

  -i,   --input : A genbank file to use as an input

Optional arguments:

  -o,   --output: Where should files be written
  -b,   --batch : The barcode for the oligo plate (defaults to 1)
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