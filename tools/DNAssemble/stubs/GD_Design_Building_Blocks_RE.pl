#!/usr/bin/env perl

**USELIB**
use Bio::GeneDesign::Basic qw(:all);
use Bio::GeneDesign::Codons qw(:all);
use Bio::GeneDesign::RestrictionEnzymes qw(:all);
use Bio::GeneDesign::BuildingBlock;
use Bio::SeqIO;
use Bio::Seq;
use Config::Auto;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;

use strict;
use warnings;

my $VERSION = '4.00';
my $GDV = "GD_Design_Building_Blocks_RE_$VERSION";

local $| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'        => \$p{INPUT},
      'output=s'       => \$p{OUTPUT},
      'length:i'      => \$p{TARBBLEN},
      'lap:i'          => \$p{TARBBLAP},
      'help'          => \$p{HELP}
);
pod2usage(-verbose=>99, -sections=>"DESCRIPTION|ARGUMENTS") if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $GD = configure_GeneDesign("**CONFLOC**");

die "\n ERROR: You must supply an input FASTA file.\n"
  if (! $p{INPUT});

die "\n ERROR: $p{INPUT} does not exist.\n"
  if (! -e $p{INPUT});

die "\n ERROR: You must supply an output path.\n"
  if (! $p{OUTPUT});

$p{TARBBLEN} = $p{TARBBLEN}  ||  740;
die "\n ERROR: building block size is outside of allowable range.\n"
  if ($p{TARBBLEN} < 400 || $p{TARBBLEN} > 1000);

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $RES = define_sites($GD->{enzyme_list});
my @enzes = values %$RES;

$p{TARBBLAP} = $p{TARBBLAP}  ||  40;
$p{BBLEN_MIN} = 350;
$p{BBLEN_MAX} = 800;

my $FIN = Bio::SeqIO->new(-file => $p{INPUT}, -format => 'FASTA');

my $inputfilename = basename($p{INPUT});
$inputfilename =~ s{\.[^.]+$}{};
my $outputfilename = $inputfilename . "_BB.FASTA";
my $outputpath = $p{OUTPUT} . "/" . $outputfilename;
open (my $OUTFH, ">" . $outputpath );
my $FOUT = Bio::SeqIO->new(-fh=>$OUTFH, -format=>'FASTA', -width => 80);
my $UOUT;

################################################################################
################################### CARVING ####################################
################################################################################
while ( my $obj = $FIN->next_seq() )
{
  print "\n";
  warn "\n WARNING: " . $obj->id . " will be a single building block.\n"
    if (length($obj->seq) < $p{TARBBLEN});
  my $chunkname = $obj->id;
  my $chseq = $obj->seq;
  my $chlen = length($chseq);
  my $tarbblen = $p{TARBBLEN};
  $tarbblen    = $chlen if ($chlen < 1.5 * $p{TARBBLEN});
  $tarbblen    = ($chlen / 2) if ($chlen <= .5 * $p{TARBBLEN});
  my @BBs;

  ##-Find the unique sites
  my $lastenz  = undef;
  my $tarmult  = .96 * $tarbblen;
  my $RE_STATUS = define_site_status($chseq, \@enzes);
  my %epos;
  my @set = grep { $RE_STATUS->{$_} == 1 } keys %$RES;
     @set = grep { $RES->{$_}->type !~ /b/ } @set;
     @set = grep { $RES->{$_}->class eq "IIP" } @set;
  foreach my $enzyme (map {$RES->{$_}} @set)
  {
    my $positions = $enzyme->positions($chseq);
    my @temparr = keys %$positions;
    my $pos = shift(@temparr);
    my $site = $pos + $enzyme->len;
    $epos{$enzyme->id} = $site;
  }
  die ("$chunkname does not contain enough unique restriction enzyme "
      . "recognition sites.\n")
      if (scalar(keys %epos) < ($chlen / $tarmult)+1);

  my ($y, $start, $last) = (1, 0, .3 * $tarbblen);
  if ($chlen > ($tarbblen + (0.2 * $tarbblen)))
  {
    my $target = $tarbblen;
    my $limit = int($chlen / $tarmult) + 1;
    while ($target < ($limit * $tarmult))
    {
      my @potents = sort {abs($epos{$a}-$target) <=> abs($epos{$b}-$target)}
                    grep {abs($epos{$_}-$target) <= 2 * $tarbblen}
                    grep {$epos{$_} > $last}
                    @set;
      unless (scalar(@potents))
      {
        warn ("There are not enough unique restriction enzyme recognition "
             . "sites in " . $obj->id . "\n");
        next;
      }
      my ($chosenpos, $enz) = ($epos{$potents[0]}, $potents[0]);
      my $lap = int(.5 * ($p{TARBBLAP} - $RES->{$enz}->len) + 1);
      my $bbSTOP = $chosenpos + $lap;
      my $bbSEQ = substr($chseq, $start, $bbSTOP - $start + 1);
      my $bbSTART = $start + 1;
      my $t = $y;
      $t = "0" . $t while (length($t) < 2);
      $y++;
      $start = $bbSTOP - $p{TARBBLAP} - 1;
      $target = $epos{$potents[0]} + $tarbblen;
      push @BBs, Bio::GeneDesign::BuildingBlock->new(
            -seq => $bbSEQ,
            -start => $bbSTART,
            -stop => $bbSTOP,
            -enzyme3 => $enz,
            -enzyme5 => $lastenz,
            -id => $obj->id . "." . $t,
            -comment => $GDV);
      $last = $chosenpos;
      $lastenz = $enz;
      last if ($chlen - $target <= ($tarbblen-(.5 * $tarbblen)));
    }
  }
  {
    my $t = $y;
    $t = "0" . $t while (length($t) < 2);
    push @BBs, Bio::GeneDesign::BuildingBlock->new(
          -seq => substr($chseq, $start),
          -start => $start + 1,
          -stop => $chlen,
          -enzyme5 => $lastenz,
          -id => $obj->id . "." . $t,
          -comment => $GDV);
  }

  foreach my $bb (@BBs)
  {
    my $newobj = Bio::Seq->new( 
                    -seq  => $bb->seq, 
                    -id   => $bb->id, 
                    -desc => $bb->description);
    $FOUT->write_seq($newobj);
  }
}

exit;

__END__

=head1 NAME

  GD_Design_Building_Blocks_RE.pl

=head1 VERSION

  Version 4.00

=head1 DESCRIPTION

    The Design_Building_Blocks_RE script will break each nucleotide sequence it
    is given into evenly sized Building Blocks, which can be composed of sets of
    overlapping oligos.

    Output will be tagged with the GDbb suffix. Default values are assumed for 
    every parameter not provided.

    Any sequence provided that is less than one and a half times the target
    building block size will not be divided.

    Restriction Enzyme Overlap: Building Blocks will overlap on unique
    restriction enzyme recognition sites. An overlap parameter may be defined to
    determine how much sequence overlaps, including the restriction enzyme 
    recognition site itself. This algorithm makes use of existing sites only and
    does not add or modify sites. If there are not enough evenly spaced, unique
    RE sites the algorithm will fault. This algorithm is not suitable for 
    dividing sequence over 12000 bp long.

=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a FASTA file containing protein sequences.
    -o,  --output : a path in which to deposit building block sequences.

  Optional arguments:

    -le, --length: the length in bp of building blocks, between 400 and 1000
        (default is 740)
    -la, --lap: the target overlap between building blocks.  This parameter is
       default is 40 otherwise
    -h,  --help: display this message

=cut