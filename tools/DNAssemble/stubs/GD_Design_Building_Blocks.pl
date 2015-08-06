#!/usr/bin/env perl

use Bio::GeneDesign;
use Bio::GeneDesign::Basic qw(:GD);
use Bio::Annotation::Comment;
use Getopt::Long;
use Pod::Usage;
use POSIX;
use Carp;

use strict;
use warnings;

my $VERSION = '5.00';
my $GDV = "GD_Design_Building_Blocks_$VERSION";
my $GDS = "_BB";

local $| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'         => \$p{INPUT},
      'output=s'        => \$p{OUTPUT},
      'format=s'        => \$p{FORMAT},
      'length:i'        => \$p{TARBBLEN},
      'maxlength:i'     => \$p{MAXBBLEN},
      'minlength:i'     => \$p{MINBBLEN},
      'lap:i'           => \$p{TARBBLAP},
      'list:s'          => \$p{LIST},
      'verbose'         => \$p{VERBOSE},
      'cloningvector=s' => \$p{CLONVEC},
      'destvector=s'    => \$p{DESTVEC},
      'simple'          => \$p{SIMPLE},
      'stitch:i'        => \$p{STITCH},
      'help'            => \$p{HELP}
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();

$p{LIST} = $p{LIST} || "outside";
  
#The input file must exist and be a format we care to read.
die "\n GDERROR: You must supply an input file.\n"
  if (! $p{INPUT});
  
my ($iterator, $filename, $suffix) = $GD->import_seqs($p{INPUT});

$p{FORMAT} = $p{FORMAT} || $suffix || "genbank";

#The output path must exist, and we'll need it to end with a slash
$p{OUTPUT} = $p{OUTPUT} || ".";
$p{OUTPUT} .= "/" if (substr($p{OUTPUT}, -1, 1) !~ /[\/]/);
die "\n GDERROR: $p{OUTPUT} does not exist.\n"
  if ($p{OUTPUT} && ! -e $p{OUTPUT});

$p{TARBBLEN} = $p{TARBBLEN} || 1000;
$p{MAXBBLEN} = $p{MAXBBLEN} || 1023;
$p{MINBBLEN} = $p{MINBBLEN} || 200;
die "\n ERROR: building block size is outside of allowable range.\n"
  if ($p{TARBBLEN} < $p{MINBBLEN} || $p{TARBBLEN} > $p{MAXBBLEN});

################################################################################
################################# CONFIGURING ##################################
################################################################################
my @seqstowrite;
my @seqstoignore;
my @seqstofail;

$p{TARBBLAP} = $p{TARBBLAP}  ||  40;

#Set up two restriction enzyme sets for filtering enzymes through
my $ORE = Bio::GeneDesign->new();
$ORE->set_restriction_enzymes(-enzyme_set => $p{LIST});
my $ARE = Bio::GeneDesign->new();
$ARE->set_restriction_enzymes(-enzyme_set => "blunts");

#Set up vectors, if necessary; note which enzymes should be excluded
my %excludes;
my ($destvec, $clonvec) = (undef, undef);
my (%fpnos, %tpnos) = ((), ());
if ($p{DESTVEC})
{
  $destvec = $GD->load_vector(-name => $p{DESTVEC});
  $excludes{$_}++ foreach $destvec->enzyme_list;

  my $stat5 = $ARE->restriction_status(-sequence => $destvec->chew5);
  %fpnos = map {$_ => 1} grep {$stat5->{$_} != 0} keys %{$stat5};

  my $stat3 = $ARE->restriction_status(-sequence => $destvec->chew3);
  %tpnos = map {$_ => 1} grep {$stat3->{$_} != 0} keys %{$stat3};
}
if ($p{CLONVEC})
{
  $clonvec = $GD->load_vector(-name => $p{CLONVEC});
  $excludes{$_}++ foreach $clonvec->enzyme_list;
}
my @enzes = keys %excludes;

my @outsides = values %{$ORE->enzyme_set};
$ARE->add_to_enzyme_set(-enzymes => \@outsides);


################################################################################
################################### CARVING ####################################
################################################################################
while ( my $obj = $iterator->next_seq() )
{
  my $chunkname = $obj->id;
  my $chseq = uc $obj->seq;
  my $chlen = length($chseq);
  my $sublet = 'A';
  
  #Check that the chunk has not already been carved
  my @bbs = grep {$_->primary_tag eq "building_block"} $obj->get_SeqFeatures;
  if (scalar(@bbs))
  {
    print "\t" . $obj->id . " already has building blocks... skipping\n";
    push @seqstoignore, $obj;
    push @seqstowrite, $obj;
    next;
  }

  #Attempt to carve the chunk
  print "Working on $chunkname ($chlen bp)...\n";
  my $tarbblen = $chlen < $p{TARBBLEN}  ? $chlen  : $p{TARBBLEN};
  my $BBS = $GD->carve_building_blocks(
          -sequence       => $obj,
          -algorithm      => "enzyme",
          -target_length  => $tarbblen,
          -max_length     => $p{MAXBBLEN},
          -min_length     => $p{MINBBLEN},
          -overlap_length => $p{TARBBLAP},
          -excludes       => \@enzes,
          -fpexcludes     => \%fpnos,
          -tpexcludes     => \%tpnos,
          -verbose        => $p{VERBOSE}
  );
  my $num = scalar(@$BBS);
  unless ($num > 0)
  {
    print "\t", $obj->id . " cannot be made into building blocks... skipping\n";
    push @seqstofail, $obj;
    push @seqstowrite, $obj;
    next;
  }
  
  
  #Decide which chunks should have flanking enzymes
  my @FINALS;
  if ($num > 5)
  {
    my @proced;
    my $f = floor $num / 4;
    my $t = (4 - ($num % 4)) % 4;
    $f = $f - ($t - 1) if ($t > 1);
    my $arr = 4 x $f . 3 x $t;
    foreach my $size (split(q{}, $arr))
    {
      my @subset;
      for (my $x = 1; $x <= $size; $x++)
      {
        my $bb = shift @$BBS;
        push @subset, $bb;
      }

      my $first = shift @subset;
      my $last = pop @subset;
      my $name = join(q{}, $first->get_tag_values('name')) . "-"
               . join(q{}, $last->get_tag_values('name'));
      my $subseq = substr($chseq, $first->start -1, $last->end - $first->start + 1);
      
      if (! scalar @proced)
      {
        my $inenzid = join(q{}, $last->get_tag_values("enz_t"));
        my $inenz = $ARE->enzyme_set->{$inenzid};
        my $prepend = $destvec  ? $destvec->chew5  . $subseq  : $subseq;
        my $outenzarr = pick_enzyme($prepend, $inenz, $clonvec, undef);
        carp ("\t$name has no outside enzyme available") unless ($outenzarr);
          
        clean_first($first, $outenzarr, 0, $clonvec, $destvec);
        clean_last($last, undef, 1, $clonvec, undef);
      }
      elsif (! scalar @$BBS)
      {
        my $inenzid = join(q{}, $first->get_tag_values("enz_f"));
        my $inenz = $ARE->enzyme_set->{$inenzid};
        my $postpend = $destvec  ? $subseq . $destvec->chew3  : $subseq;
        my $outenzarr = pick_enzyme($postpend, $inenz, undef, $clonvec);
        carp ("\t$name has no outside enzyme available") unless ($outenzarr);
          
        clean_first($first, undef, 1, $clonvec, undef);
        clean_last($last, $outenzarr, 0, $clonvec, $destvec);
      }
      else
      {
        clean_first($first, undef, 1, $clonvec, undef);
        clean_last($last, undef, 1, $clonvec, undef);
      }
      $first->add_tag_value("subset", $sublet);
      $last->add_tag_value("subset", $sublet);
      push @proced, $first, $last;

      foreach my $bb (@subset)
      {
        $bb->add_tag_value("subset", $sublet);
        clean_mundane($bb, $clonvec);
      }
      push @proced, @subset;
      $sublet++;
    }
    push @FINALS, @proced;
  }

  elsif ($num > 1)
  {
    my $first = shift @$BBS;
    my $last = pop @$BBS;
    my @mundane = @$BBS;
    if ($p{SIMPLE})
    {
      clean_first($first, undef, 0, $clonvec, $destvec);
      clean_last($last, undef, 0, $clonvec, $destvec);
    }
    else
    {
      my $name = join(q{}, $first->get_tag_values('name')) . "-"
               . join(q{}, $last->get_tag_values('name'));
      my $rseq = $destvec ? $destvec->chew5 . $chseq . $destvec->chew3 : $chseq;
      $rseq = $clonvec ? $clonvec->chew5 . $rseq . $clonvec->chew3 : $rseq;
      my $outenzarr = pick_enzyme($rseq);
      carp ("\t$name has no outside enzyme available") unless ($outenzarr);
    
      clean_first($first, $outenzarr, 0, $clonvec, $destvec);
      clean_last($last, $outenzarr, 0, $clonvec, $destvec);
    }
    foreach my $bb (@mundane)
    {
      clean_mundane($bb, $clonvec);
    }
    push @FINALS, $first;
    push @FINALS, @mundane;
    push @FINALS, $last;
  }

  else
  {
    my $bb = shift @$BBS;
    if ($p{SIMPLE})
    {
      clean_only($bb, undef, $clonvec, $destvec);
    }
    else
    {
      my $name = join(q{}, $bb->get_tag_values('name'));
      my $rseq = $destvec ? $destvec->chew5 . $chseq . $destvec->chew3 : $chseq;
      $rseq = $clonvec ? $clonvec->chew5 . $rseq . $clonvec->chew3 : $rseq;
      my $outenzarr = pick_enzyme($rseq);
      carp ("\t$name has no outside enzyme available") unless ($outenzarr);
    
      clean_only($bb, $outenzarr, $clonvec, $destvec);
    }
    push @FINALS, $bb;
  }
  
  if ($p{DESTVEC} || $p{CLONVEC})
  {
    my $vector = $p{DESTVEC}  ? $p{DESTVEC}  : $p{CLONVEC};
    #$obj->id("(" . $vector . ")$chunkname");
    
    my $collec = $obj->annotation;
    my $comment = Bio::Annotation::Comment->new();
    $comment->text("destination_vector = $vector");
    $collec->add_Annotation('comment', $comment);
    $obj->annotation($collec);
  }
  foreach my $bb (sort {$a->start <=> $b->start} @FINALS)
  {
    if ($p{STITCH})
    {
      my $bbseq = substr($chseq, $bb->start - 1, $bb->end - $bb->start + 1);
      my ($lprimer, $rprimer) = $GD->make_amplification_primers(
          -sequence    => $bbseq,
          -temperature => $p{STITCH}
      );
      unless ($bb->has_tag("enz_f"))
      {
        $bb->add_tag_value("stitch_left", $lprimer);
      }
      unless ($bb->has_tag("enz_t"))
      {
        $bb->add_tag_value("stitch_right", $rprimer);
      }
    }
    $bb->add_tag_value("vector", $p{CLONVEC}) if ($p{CLONVEC});
    #$bb->add_tag_value("gdversion", $GDV);
    $obj->add_SeqFeature($bb);
  }
  push @seqstowrite, $obj;
}

#Maybe report

print "\n\n";
my $outputfilename = $filename  . $GDS . "." . $p{FORMAT};
my $reportpath = $p{OUTPUT} . $filename . "_BBreport.txt";
open (my $REP, '>', $reportpath);
if (scalar @seqstoignore)
{
  print $REP $_->id . " : IGNORED FOR BUILDING BLOCKS\n" foreach @seqstoignore;
}
if (scalar @seqstofail)
{
  print $REP $_->id . " : FAILED BUILDING BLOCKS\n" foreach @seqstofail;
}
close $REP;
if (scalar @seqstowrite)
{
  $GD->export_seqs(
          -filename   => $outputfilename,
          -path       => $p{OUTPUT},
          -format     => $p{FORMAT},
          -sequences  => \@seqstowrite
  );
}

print "\n";
print "Wrote " . $p{OUTPUT} . "$outputfilename\n";
print "Wrote $reportpath\n";
print "\n";
print $GD->attitude() . " brought to you by $GDV\n\n";

exit;

sub clean_first
{
  my ($bb, $outenzarr, $inenz, $clonvec, $destvec) = @_;
  my ($fprefix, $fsuffix) = (q{}, q{});
  if ($outenzarr)
  {
    $fprefix = $outenzarr->[1];
    $bb->add_tag_value("enz_f", $outenzarr->[0]);
  }
  if ($destvec)
  {
    $fprefix .= $destvec->chew5;
  }
  if ($inenz)
  {
    $fprefix = join(q{}, $bb->get_tag_values("zne_f"));
    $fprefix = $GD->replace_ambiguous_bases($fprefix);
  }
  if ($clonvec)
  {
    $fprefix = $clonvec->chew5 . $fprefix;
    $fsuffix .= $clonvec->chew3;
  }

  $bb->remove_tag("enz_t") if ($bb->has_tag("enz_t"));
  $bb->remove_tag("zne_t") if ($bb->has_tag("zne_t"));
  $bb->remove_tag("zne_f") if ($bb->has_tag("zne_f"));
  $bb->add_tag_value("BBsuffix", $fsuffix) if ($fsuffix);
  $bb->add_tag_value("BBprefix", $fprefix) if ($fprefix);
  return;
}

sub clean_last
{
  my ($bb, $outenzarr, $inenz, $clonvec, $destvec) = @_;
  my ($lprefix, $lsuffix) = (q{}, q{});
  if ($outenzarr)
  {
    $lsuffix = $GD->complement($outenzarr->[1], 1);
    $bb->add_tag_value("enz_t", $outenzarr->[0]);
  }
  if ($destvec)
  {
    $lsuffix = $destvec->chew3 . $lsuffix;
  }
  if ($inenz)
  {
    my $znet = join(q{}, $bb->get_tag_values("zne_t"));
    $znet = $GD->replace_ambiguous_bases($znet);
    $lsuffix .= $znet;
  }
  if ($clonvec)
  {
    $lprefix = $lprefix . $clonvec->chew5;
    $lsuffix .= $clonvec->chew3;
  }
  $bb->remove_tag("enz_f") if ($bb->has_tag("enz_f"));
  $bb->remove_tag("zne_f") if ($bb->has_tag("zne_f"));
  $bb->remove_tag("zne_t") if ($bb->has_tag("zne_t"));
  $bb->add_tag_value("BBsuffix", $lsuffix) if ($lsuffix);
  $bb->add_tag_value("BBprefix", $lprefix) if ($lprefix);
  return;
}

sub clean_only
{
  my ($bb, $outenzarr, $clonvec, $destvec) = @_;
  
  my ($prefix, $suffix) = (q{}, q{});
  if ($outenzarr)
  {
    $prefix = $outenzarr->[1];
    $suffix = $GD->complement($outenzarr->[1], 1);
    $bb->add_tag_value("enz_t", $outenzarr->[0]);
    $bb->add_tag_value("enz_f", $outenzarr->[0]);
  }
  if ($clonvec)
  {
    $suffix .= $clonvec->chew3;
    $prefix = $clonvec->chew5 . $prefix;
    if ($destvec)
    {
      $prefix .= $destvec->chew5;
      $suffix = $destvec->chew3 . $suffix;
    }
  }
  $bb->add_tag_value("BBsuffix", $suffix) if ($suffix);
  $bb->add_tag_value("BBprefix", $prefix) if ($prefix);
  return;
}

sub clean_mundane
{
  my ($bb, $clonvec) = @_;
  my ($nprefix, $nsuffix) = (q{}, q{});
  $bb->remove_tag("zne_t") if ($bb->has_tag("zne_t"));
  $bb->remove_tag("enz_t") if ($bb->has_tag("enz_t"));
  if ($clonvec)
  {
    $nprefix = $clonvec->chew5;
    $nsuffix .= $clonvec->chew3;
  }
  $bb->remove_tag("zne_f") if ($bb->has_tag("zne_f"));
  $bb->remove_tag("enz_f") if ($bb->has_tag("enz_f"));
  
  $bb->add_tag_value("BBprefix", $nprefix) if ($nprefix);
  $bb->add_tag_value("BBsuffix", $nsuffix) if ($nsuffix);
  return;
}

sub pick_enzyme
{
  my ($seq, $enz, $cprepend, $cpostpend) = @_;
  
  my $substats = $ORE->restriction_status(-sequence => $seq);
  my @candidates =  sort {$a->score <=> $b->score}
                    grep {$substats->{$_->id} == 0}
                    map {$ORE->enzyme_set->{$_}}
                    keys %{$substats};
  my $outenz = $candidates[0];
  
  return q{} unless (scalar @candidates);
  
  if ($enz)
  {
    my @filtered = grep {$enz->temp == $_->temp} @candidates;
    $outenz = scalar @filtered ? $filtered[0] : $outenz;
    
    @filtered = grep {$_->common_buffers($enz, 1)} @filtered;
    $outenz = scalar @filtered ? $filtered[0] : $outenz;
    
    if ($cprepend || $cpostpend)
    {
      my @refiltered = ();
      foreach my $potenz (@filtered)
      {
        if ($cprepend)
        {
          my $tempseq = $cprepend->chew5 . $potenz->seq;
          my $positions = $enz->positions($tempseq);
          push @refiltered, $potenz if (scalar keys %{$positions} == 0);
        }
        if ($cpostpend)
        {
          my $tempseq = $potenz->seq . $cpostpend->chew3;
          my $positions = $enz->positions($tempseq);
          push @refiltered, $potenz if (scalar keys %{$positions} == 0);
        }
      }
      @filtered = @refiltered;
      $outenz = scalar @filtered ? $filtered[0] : $outenz;
    }
  }
  
  my $outenzseq = $ORE->replace_ambiguous_bases($outenz->seq);
  my $rand = $ORE->random_dna(-length => $outenz->inside_cut);
  $outenzseq = $outenzseq . $rand;
  my $finalseq = $seq . $outenzseq . $seq;
  #Check if the randomly generated sequence contains the site twice !
  my $positions = $GD->positions(-sequence => $outenzseq, -query => $finalseq);
  while (scalar keys %{$positions} != 1)
  {
    $outenzseq = $ORE->replace_ambiguous_bases($outenz->seq);
    $rand = $ORE->random_dna(-length => $outenz->inside_cut);
    $outenzseq = $outenzseq . $rand;
    $finalseq = $seq . $outenzseq . $seq;
    $positions = $GD->positions(-sequence => $finalseq, -query => $outenz);
  }
  
  return [$outenz->id, $outenzseq];
}

__END__

=head1 NAME

  GD_Design_Building_Blocks.pl

=head1 VERSION

  Version 5.00

=head1 DESCRIPTION

    The Design_Building_Blocks_JGI script will break each nucleotide sequence it
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

    -i,  --input : a file containing DNA sequences.

  Optional arguments:
  
    -o, --output : a path in which to deposit building block sequences
    -f, --format : default genbank
    -le, --length: the length in bp of building blocks (default is 1000)
    -maxle, --maxlength: the maximum length a building block is allowed to be.
    -minle, --minlength: the minimum length a building block is allowed to be.
    -la, --lap: the target overlap between building blocks. (default is 40)
    -c, --cloningvector: A vector in config/vectors that is to serve as the
        cloning vector
    -d, --destvector: A vector in config/vectors that is to serve as the
        destination vector
    -si, --simple: If a chunk has five or fewer building blocks, don't bother
        assigning enzymes to them
    -st, --stitch: The target melting temperature of assembly primers for
        the assembly of building blocks.
    --verbose : show deliberations happening
    -h,  --help: display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012, GeneDesign developers
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