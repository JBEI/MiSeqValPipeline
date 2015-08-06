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
my $GDV = "GD_Design_Building_Blocks_USER_$VERSION";

local $| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'        => \$p{INPUT},
      'output=s'       => \$p{OUTPUT},
      'length:i'       => \$p{TARBBLEN},
      'utemp:i'        => \$p{TARUSRMEL},
      'ulength:s'      => \$p{USRUNILEN},
      'help'           => \$p{HELP}
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
my %useroliname = (1 => "LT", 2 => "LU", 3 => "RT", 0 => "RU");

$p{TARUSRMEL} = $p{TARUSRMEL}  ||  56;
$p{USRUNILEN} = $p{USRUNILEN}  || "5,7,9,11";
my @usrlens = split(",", $p{USRUNILEN});
$p{USRUNILEN} = \@usrlens;
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
my $uoutputfilename = $inputfilename . "_BB_USERs.FASTA";
my $uoutputpath = $GD->{output_dir} . "/" . $uoutputfilename;
open (my $UOUTFH, ">" . $uoutputpath );
$UOUT = Bio::SeqIO->new(-fh=>$UOUTFH, -format=>'FASTA', -width => 80);

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

  if ($chlen == $tarbblen)
  {
    push @BBs, Bio::GeneDesign::BuildingBlock->new(
         -seq => $chseq,
         -start => 1,
         -stop => $chlen,
         -id => $obj->id . ".01",
         -comment => $GDV);
  }
  else
  {
    my %rank = ( 9 => 0, 7 => 1, 5 => 2, 3 => 3, 11 => 4, 13 => 5);

    my (%USites, %UCoors) = ((), ());
    my (@UniUsers, @Collection) = ((), ());

   ## Load two hashes with the number of times ANxT is seen and where the last
   ## one was seen. Parse for unique ones, store in array.
    foreach my $tiv (@{$p{USRUNILEN}})
    {
      my @sites;
      $sites[0] = "A" . ("N" x $tiv) . "T";
      my $comp = complement($sites[0], 1);
      push @sites, $comp if ($comp ne $sites[0]);
      foreach my $sit (@sites)
      {
        my $exp = regres($sit, 1);
        while ($chseq =~ /($exp)/ig)
        {
          my $sitsta = (pos $chseq) + 1;
          if ($sitsta ne '')
          {
            $USites{$1}++;
            $UCoors{$1} = $sitsta - length($1);
          }
        }
      }
    }
    my @possibles = grep {$USites{$_} == 1} keys %USites;
    @possibles = grep {! exists $USites{complement($_)}} @possibles;
    foreach my $tiv (@possibles)
    {
      my $user = {};
      $user->{SEQ} = $tiv;
      $user->{LEN} = length($user->{SEQ});
      $user->{NUMBER} = length($user->{SEQ}) - 2;
      $user->{START} = $UCoors{$user->{SEQ}};
      push @UniUsers, $user;
    }

   ## Adjust the target building block size so as to avoid outliers.
    my $tarnum = int(length($chseq) / $tarbblen);
    my $diff = length($chseq) - ($tarnum * $tarbblen);
    $tarbblen = (($diff / $tarnum) * 12.5 >= $tarbblen)
              ?  int(length($chseq) / ($tarnum + 1) + .5)
              :  $tarbblen + int($diff / $tarnum);

   ## Pick the unique sites as close as possible to the requested intervals.
    my ($last, $y, $tar) = (0, 1, $tarbblen);
    while ($tar < length($chseq))
    {
      my ($seen, $door, $int) = (0, 1, 1);
      while ($seen == 0)
      {
        my @g = sort {abs($a->{START} - $tar) <=> abs($b->{START} - $tar)}
                grep {abs($_->{START} - $tar) <= 2 * $tarbblen}
                grep {$_->{START} > $last}
                @UniUsers;
        if (scalar(@g))
        {
          @g = sort {$rank{$a->{LEN}} <=> $rank{$b->{LEN}}} @g;
          my $currchoice = $g[0];
          $last = $tar;
          $tar = $currchoice->{START} + $tarbblen;
          push @Collection, $currchoice;
          $seen = 1;
        }
        $tar = $tar - $int if ($door == 1);
        $tar = $tar + $int if ($door == 0);
        $door = abs($door-1);
        $int++;
        if ($tar <= ($last + (.4 * $tarbblen)))
        {
          $last += ($tarbblen * .5);last;
        }
        $seen = 0;
      }
      $tar = $tarbblen + $last;
    }

   ## Form chunk objects, fill with user primers and oligos
    my ($lastfound, $laststart, $lastlength) = (1, 0, 0);
    foreach my $user (@Collection)
    {
      my @Users;
      my $count = $y-1;
      my $length = $user->{START} - $lastfound + $user->{NUMBER} + 2;
      my $bbSEQ = substr($chseq, $lastfound - 1, $length);
      if (length($bbSEQ) < length(substr($chseq, $lastfound - 1))
          && $count == int((length($chseq) / $tarbblen) + .5) - 1)
      {
        $bbSEQ = substr($chseq, $lastfound-1);
      }
      my $remain = length(substr($chseq, $user->{START}-1));
      if ($remain < $tarbblen / 4)
      {
        $bbSEQ = substr($chseq, $lastfound - 1);
      }
      my $bbLEN = length($bbSEQ);
      my $t = $y;
      $t = "0" . $t while (length($t) < 2);
      my $pri_len = 0;
      if ($count > 0)
      {
        $Users[0] = substr($chseq, $laststart, $lastlength + 2);
        until (melt($Users[0]) >= $p{TARUSRMEL})
        {
          $pri_len++;
          $Users[0] = substr($chseq, $laststart, $pri_len)
        }
        $Users[1] = substr($Users[0], 0, $lastlength+1) . 'U'
                  . substr($Users[0], $lastlength + 2);
      }
      else
      {
        $Users[0] = substr($chseq, 0, 5);
        until (melt($Users[0]) >= $p{TARUSRMEL})
        {
          $pri_len++;
          $Users[0] = substr($chseq, 0, $pri_len)
        }
        $Users[1] = "-";
      }
      $pri_len = 0;
      if ($count != (@Collection-1))
      {
        $Users[2] = substr($chseq, $user->{START} - 1, $user->{NUMBER} + 2);
        until (melt($Users[2])) >= $p{TARUSRMEL})
        {
          $pri_len++;
          $Users[2] = substr($chseq,
                             $user->{START} - 1 - $pri_len,
                             $user->{NUMBER} + 2 + $pri_len);
        }
        $Users[2] = complement($Users[2], 1);
        $Users[3] = substr($Users[2], 0, $user->{NUMBER} + 1) . 'U'
                  . substr($Users[2], $user->{NUMBER} + 2);
        $laststart = $user->{START} - 1;
        $lastlength = $user->{NUMBER};
      }
      else
      {
        $Users[2] = substr($chseq, -10);
        until (melt($Users[2]) >= $p{TARUSRMEL})
        {
          $pri_len++;
          $Users[2] = substr($chseq, -(10 + $pri_len));
        }
        $Users[2] = complement($Users[2], 1);
        $Users[3] = "-";
      }
      $y++;
      $lastfound = $user->{START};#$lastseq = $user->{SEQUENCE};

      push @BBs, Bio::GeneDesign::BuildingBlock->new(
           -seq => $bbSEQ,
           -start => $lastfound,
           -stop => $bbLEN + $lastfound - 1,
           -id => $obj->id . "." . $t,
           -comment => $GDV,
           -userlist => \@Users);
    }
  }

  foreach my $bb (@BBs)
  {
    my $newobj = Bio::Seq->new( 
                    -seq  => $bb->seq, 
                    -id   => $bb->id, 
                    -desc => $bb->description);
    $FOUT->write_seq($newobj);
    my $x = 1;
    foreach my $user (grep {$_ ne "-"} @{$bb->userlist})
    {
      my $desc = length($user) . "bp";
      my $id = $bb->id . ".u" . $useroliname{$x%4};
      my $newuser = Bio::Seq->new( -seq => $user, -id => $id, -desc => $desc);
      $UOUT->write_seq($newuser);
      $x++;
    }
  }
}

exit;

__END__

=head1 NAME

  GD_Design_Building_Blocks_USER.pl

=head1 VERSION

  Version 4.00

=head1 DESCRIPTION

    The Design_Building_Blocks_USER script will break each nucleotide sequence 
    it is given into evenly sized Building Blocks, which can be composed of sets
    of overlapping oligos.

    Output will be tagged with the GDbb suffix. Default values are assumed for 
    every parameter not provided.

    Any sequence provided that is less than one and a half times the target
    building block size will not be divided.

    Building Blocks will overlap on A(N)xT sequences, so as to be compatible 
    with a uracil exicision (USER) assembly protocol. The width of overlap is 
    user definable, and a melting temperature may be defined for USER oligos. 
    Input sequences must be at least 1000 bp long.

=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a FASTA file containing protein sequences.
    -o,  --output : a path in which to deposit building block sequences.

  Optional arguments:

    -le, --length: the length in bp of building blocks, between 400 and 1000
        (default is 740)
    -ut, --utemp: the target melting temperature in degrees C of user oligos
        (default is 56)
    -ul, --ulength: the target length of unique sequence in user oligos. This
        refers to the number of Ns in the sequence A(N)xT.  More than one may be
        provided seperated by commas from the set 5, 7, 9, or 11.
        (default is 5,7,9,11)
    -h,  --help: display this message

=cut