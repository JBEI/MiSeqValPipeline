#!/usr/bin/env perl

**USELIB**
use Bio::GeneDesign::Basic qw(:all);
use Bio::GeneDesign::Codons qw(:all);
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
my $GDV = "GD_Design_Building_Blocks_JHU_$VERSION";

local $| = 1;

##Get Arguments
my %pa = ();
GetOptions (
      'input=s'        => \$pa{INPUT},
      'output=s'       => \$pa{OUTPUT},
      'length:i'      => \$pa{TARBBLEN},
      'lap:i'          => \$pa{TARBBLAP},
      'temp:i'        => \$pa{TARCHNMEL},
      'tolerance:f'    => \$pa{CHNMELTOL},
      'pin=s'         => \$pa{PIN},
      'help'          => \$pa{HELP}
);
pod2usage(-verbose=>99, -sections=>"DESCRIPTION|ARGUMENTS") if ($pa{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $GD = configure_GeneDesign("**CONFLOC**");

die "\n ERROR: You must supply an input FASTA file.\n"
  if (! $pa{INPUT});

die "\n ERROR: $pa{INPUT} does not exist.\n"
  if (! -e $pa{INPUT});

die "\n ERROR: You must supply an output path.\n"
  if (! $pa{OUTPUT});

die "\n ERROR: $pa{PIN} does not exist.\n"
  if ($pa{PIN} && ! -e $pa{PIN});

$pa{TARBBLEN} = $pa{TARBBLEN}  ||  740;
die "\n ERROR: building block size is outside of allowable range.\n"
  if ($pa{TARBBLEN} < 400 || $pa{TARBBLEN} > 1000);

################################################################################
################################# CONFIGURING ##################################
################################################################################
$pa{TARBBLAP} = $pa{TARBBLAP}  ||  40;
$pa{TARCHNMEL} = $pa{TARCHNMEL}  ||  56;
$pa{CHNMELTOL} = $pa{CHNMELTOL}  ||  2.5;
$pa{BBLEN_MIN} = 350;
$pa{BBLEN_MAX} = 800;

my $FIN = Bio::SeqIO->new(-file => $pa{INPUT}, -format => 'FASTA');

my $inputfilename = basename($pa{INPUT});
$inputfilename =~ s{\.[^.]+$}{};
my $outputfilename = $inputfilename . "_BB.FASTA";
my $outputpath = $pa{OUTPUT} . "/" . $outputfilename;
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
    if (length($obj->seq) < $pa{TARBBLEN});
  my $chunkname = $obj->id;
  my $chseq = $obj->seq;
  my $chlen = length($chseq);
  my $tarbblen = $pa{TARBBLEN};
  $tarbblen    = $chlen if ($chlen < 1.5 * $pa{TARBBLEN});
  $tarbblen    = ($chlen / 2) if ($chlen <= .5 * $pa{TARBBLEN});
  my @BBs;

  if ($pa{PIN})
  {
    my %ssrseqs = ();
    my $PIN = Bio::SeqIO->new(-file => $pa{PIN}, -format => 'FASTA');
    while (my $seq = $PIN->next_seq())
    {
      $ssrseqs{$seq->id} = $seq->seq;
    }

   ## Mask for site specific recognition sites (should be hairpins, later)
    my @MASK = (0) x $chlen;
    foreach my $hpseq (values %ssrseqs)
    {
      my @positions = pattern_finder($chseq, $hpseq, 1);
      foreach my $x (@positions)
      {
        for ($x - $pa{TARBBLAP} .. ($x + length($hpseq)) + $pa{TARBBLAP})
        {
          $MASK[$_]++ unless ($MASK[$_] > 0);
        }
        for ($x .. ($x + length($hpseq)) -1)
        {
          $MASK[$_]++ unless ($MASK[$_] == 2);
        }
      }
    }
    my $mask = join(q{}, @MASK);
  ## introduce breaks between pin sites that are too close
    my @positions = ();
    my $exp = qr/2[1]+2/;
    while ($mask =~ /(?=$exp)/ig)
    {
      push @positions, (pos $mask);
    }
    foreach my $spot (@positions)
    {
      my $stop = $spot+ 1;
      $stop++ while (substr($mask, $stop, 1) != 2);
      if (($stop - $spot) < 40)
      {
        print "\tPINS ARE TOO CLOSE @ $spot: ", substr($mask, $spot, $stop-$spot+1), "\n";
        substr($mask, $spot + 41, 1) = "0";
      }
      else
      {
        substr($mask, $spot + 41, 1) = "0";
      }
    }
    my @coords = substr($mask, 0, 1) eq "0"	?	(0)	:	();
  	while ($mask =~ /(?=[^0]0|0[^0])/ig)
  	{
  	  push @coords, (pos $mask)+1;
  	}
  	push @coords, length($mask) if (substr($mask, -1, 1) eq "0");
    my %filter = @coords;
    my @fkeys = sort {$a <=> $b} keys %filter;
    shift @fkeys;
    push @fkeys, $chlen;
    print $obj->id . ", $chlen bp, @fkeys\n";
    my ($bbcount, $lastpos) = (1, 0);
    while (scalar(@fkeys) && $chlen - $lastpos > $pa{BBLEN_MIN})
    {
      my $closest = shift @fkeys;
      #$closest = $chlen if ($chlen - $closest < $pa{BBLEN_MIN});
      my $distance = $closest - $lastpos;
      #      print "\t last: $lastpos, $closest, $distance\n";
      #The next one is too small - first try moving to the bit before the next site
       if ($distance < $pa{BBLEN_MIN} && $chlen - $closest > $pa{BBLEN_MIN})
      {
        while($distance < $pa{BBLEN_MIN})
        {
          my $next = shift @fkeys;
          $distance = $next - $lastpos;
          #$distance = ($fkeys[0] - 115) - $lastpos;
        }
        #next if ($distance < $pa{BBLEN_MIN} && $chlen - ($fkeys[0] - 115) > $pa{BBLEN_MIN});
      }

      if ($distance < $pa{BBLEN_MAX})# && $distance > $pa{BBLEN_MIN})
      {
        print "\tMaking $bbcount, A $distance\n";
        my $bbMASK = substr($mask, $lastpos, $distance);
        my $te = $bbcount;
        $te = "0" . $te while length($te) < 2;
        $bbcount++;
        $lastpos += $distance - $pa{TARBBLAP};
        push @BBs, Bio::GeneDesign::BuildingBlock->new(
              -seq => substr($chseq, $lastpos, $distance),
              -start => $lastpos,
              -stop => $distance + $lastpos - 1,
              -id => $obj->id . "." . $te,
              -comment => $ALGORITHM{$pa{ALGORITHM}} . " ($GDV)");
        my $maskcopy = $bbMASK;
        my $ssrcount = ($maskcopy =~ tr/2//);
        if ($ssrcount != 0 && $ssrcount != 34)
        {
          print "\t\t$te has wrong number of pin bases!!!\n";
          print $distance, "\n", $BBs[-1]->seq, "\n", $bbMASK, "\n\n";
        }
        if ($distance < $pa{BBLEN_MIN} || $distance > $pa{BBLEN_MAX})
        {
          print "\t\t$te is too short!!!\n";
          print $distance, "\n", $BBs[-1]->seq, "\n", $bbMASK, "\n\n";
        }
      }
      else
      {
        my $tar_num = ceil(($distance / ($pa{TARBBLEN} - $pa{TARBBLAP})));
        my $diff = $distance - (($tar_num * $pa{TARBBLEN}) - ($pa{TARBBLAP} * ($tar_num - 1)));
        my $tarbblen = $pa{TARBBLEN};
        #        print "\t\t$tarbblen, $tar_num\n";
        if (abs($diff) >= $tar_num)
        {
          $tarbblen = $tarbblen + int($diff / $tar_num);
          $diff = $diff - ($tar_num * (int($diff / $tar_num)));
        }
        for my $x (1..$tar_num)
        {
          my $curbblen = $tarbblen;
          $curbblen++ if ( $x <= abs($diff) && $diff > 0);
          $curbblen-- if ( $x <= abs($diff) && $diff < 0);
          print "\tMaking $bbcount, B $curbblen\n";
          my $bbMASK = substr($mask, $lastpos, $curbblen);
          my $te = $bbcount;
          $te = "0" . $te while length($te) < 2;
          push @BBs, Bio::GeneDesign::BuildingBlock->new(
                -seq => substr($chseq, $lastpos, $curbblen),
                -start => $lastpos,
                -stop => $curbblen + $lastpos - 1,
                -id => $obj->id . "." . $te,
                -comment => $ALGORITHM{$pa{ALGORITHM}} . " ($GDV)");

          my $maskcopy = $bbMASK;
          my $ssrcount = ($maskcopy =~ tr/2//);
          if ($ssrcount != 0 && $ssrcount != 34)
          {
            print "\t\t$bbcount has wrong number of pin bases!!\n";
            print $curbblen, "\n", $BBs[-1]->seq, "\n", $bbMASK, "\n\n";
          }
          if ($curbblen < $pa{BBLEN_MIN} || $curbblen > $pa{BBLEN_MAX})
          {
            print "\t\t$bbcount is too short!!\n";
            print $curbblen, "\n", $BBs[-1]->seq, "\n", $bbMASK, "\n\n";
          }
          $lastpos += $curbblen - $pa{TARBBLAP};
          $bbcount++;
        }
      }
    }
  }
  else
  {
    ## Adjust the target building block size so as to avoid outliers.
     my $tarnum = int(($chlen / ($tarbblen - $pa{TARBBLAP})) + 0.5);
     my $diff = $tarnum * $tarbblen - $pa{TARBBLAP} * ($tarnum - 1);
        $diff = $chlen - $diff;
     if (abs($diff) >= $tarnum)
     {
       $tarbblen = $tarbblen + int($diff / $tarnum);
       $diff = $diff - ($tarnum * (int($diff / $tarnum)));
     }

     my ($y, $last) = (1, 0);
     for my $cur (1..$tarnum)
     {
       my $curbblen = $tarbblen;
       $curbblen++ if ( $cur <= abs($diff) && $diff > 0);
       $curbblen-- if ( $cur <= abs($diff) && $diff < 0);
       my $t = $y;
       $t = "0" . $t while (length($t) < 2);
       $y++;
       $last += $curbblen - $pa{TARBBLAP};
       push @BBs, Bio::GeneDesign::BuildingBlock->new(
             -seq => substr($chseq, $last, $curbblen),
             -start => $last + 1,
             -stop => $curbblen + $last,
             -id => $obj->id . "." . $t,
             -comment => $ALGORITHM{$pa{ALGORITHM}} . " ($GDV)");
     }
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

  GD_Design_Building_Blocks_JHU.pl

=head1 VERSION

  Version 4.00

=head1 DESCRIPTION

    The Design_Building_Blocks_JHU script will break each nucleotide sequence it
    is given into evenly sized Building Blocks, which can be composed of sets of
    overlapping oligos.

    Output will be tagged with the GDbb suffix. Default values are assumed for 
    every parameter not provided.

    Any sequence provided that is less than one and a half times the target
    building block size will not be divided.

    Length Overlap: Building Blocks will overlap by a user-defined overlap
    length parameter. Input sequences must be at least 1000 bp long. An optional
    parameter will force building blocks to avoid the presence of more than one
    specified subsequence per building block.

=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a FASTA file containing protein sequences.
    -o,  --output : a path in which to deposit building block sequences.

  Optional arguments:

    -le, --length: the length in bp of building blocks, between 400 and 1000
        (default is 740)
    -la, --lap: the target overlap between building blocks. Default is 40
    -p,  --pin: A FASTA file containing sequences which are not allowed
        to occur more than once per building block; algorithm 1 only.
    -h,  --help: display this message

=cut