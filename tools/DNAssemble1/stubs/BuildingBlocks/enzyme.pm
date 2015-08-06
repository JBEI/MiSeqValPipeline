=head1 NAME

Bio::GeneDesign::BuildingBlocks::enzyme

=head1 VERSION

Version 5.00

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::GeneDesign::BuildingBlocks::enzyme;
require Exporter;

use Bio::GeneDesign::Basic qw(:GD);
use Bio::GeneDesign::Exceptions;
use Bio::SeqFeature::Generic;
use POSIX;

use strict;
use warnings;

our $VERSION = 5.00;

use base qw(Exporter);
our @EXPORT_OK = qw(
  _carve_building_blocks
);
our %EXPORT_TAGS =  ( GD => \@EXPORT_OK );
  
=head2 _carve_building_blocks()

=cut

sub _carve_building_blocks
{
  my ($GD, $seqobj, $tarbblen, $maxbblen, $minbblen, $tarbblap, $stitch, $excludes, $fpexcludes, $tpexcludes, $v) = @_;

  my $minbblap = 25;
  my $maxbblap = 50;
  my $chseq = uc $seqobj->seq;
  my $chlen = length($chseq);
  $tarbblen  = $chlen if ($chlen < $tarbblen);
  my $chname = $seqobj->id;

  ## Decide what size bbs to go for and how many
  ## Adjust the target building block size so as to avoid outliers.
  my $tarnum = $chlen > $maxbblen 
              ? ceil($chlen / ($tarbblen - $tarbblap))
              : 1;
  my $diff = $chlen - (($tarnum * $tarbblen) - ($tarbblap * ($tarnum - 1)));
  print "\ttarget: $tarnum bbs size $tarbblen bp rem $diff\n" if ($v);
  if (abs($diff) >= $tarnum)
  {
    my $rem = sprintf "%.0f", $diff / $tarnum;
    $tarbblen = $tarbblen + $rem;
    $diff = $diff - ($tarnum * $rem);
  }
  print "\t final: $tarnum bbs size $tarbblen bp rem $diff\n" if ($v);
 
  my @BBS;
  
  if ($tarnum > 1)
  {
    ## Decide what restriction enzyme half sites to use
    $GD->set_restriction_enzymes(-enzyme_set => "blunts");
    $GD->remove_from_enzyme_set(-enzymes => $excludes) if ($excludes);
    my $RES = $GD->enzyme_set;
    my $objstats = $GD->restriction_status(-sequence => $chseq);
    my @objcandidates = sort {$RES->{$a}->score <=> $RES->{$b}->score}
                        grep {$objstats->{$_} == 0}
                        keys %$objstats;
    print "\tCH CANDIDATES: @objcandidates\n" if ($v);
    Bio::GeneDesign::Exception::UnBBable->throw("no RE candidates!")
      unless scalar(@objcandidates);
    
    ## Create a suffix tree for each end of a bb overlap
    ## The 5' half of an RE site goes in ftree, the 3' half goes in ttree
    my $ftree = Bio::GeneDesign::PrefixTree->new();
    my $ttree = Bio::GeneDesign::PrefixTree->new();
    my %seeks;
    foreach my $cand (@objcandidates)
    {
      my $seq = uc $RES->{$cand}->seq;
      my $qes = _complement($seq, 1);
      my $half = length($seq) / 2;
      $seeks{$cand} = {len => $half};
      my $a = substr($seq, 0, $half);
      my $b = substr($seq, $half);
      my $fnucs = _amb_transcription($b);
      foreach my $fnuc (@$fnucs)
      {
        $ftree->add_prefix($fnuc, $cand);
        $seeks{$cand}->{$fnuc} = $a;
      }
      my $tnucs = _amb_transcription($a);
      foreach my $tnuc (@$tnucs)
      {
        $ttree->add_prefix($tnuc, $cand);
        $seeks{$cand}->{$tnuc} = $b;
      }
      if ($seq ne $qes)
      {
        my $c = substr($qes, 0, $half);
        my $d = substr($qes, $half);
        $fnucs = _amb_transcription($d);
        foreach my $fnuc (@$fnucs)
        {
          $ftree->add_prefix($fnuc, $cand);
          $seeks{$cand}->{$fnuc} = $c;
        }
        $tnucs = _amb_transcription($c);
        foreach my $tnuc (@$tnucs)
        {
          $ttree->add_prefix($tnuc, $cand);
          $seeks{$cand}->{$tnuc} = $d;
        }
      }
    }

    
    ## Search the construct sequence for each half RE with the suffix trees
    ## Then create every sequences that can overlap BBs
    my $laps = [];
    my %anns = ();
    my $id = 1;
    my $fref = $GD->search_prefix_tree(-tree => $ftree, -sequence => $chseq);
    my $tref = $GD->search_prefix_tree(-tree => $ttree, -sequence => $chseq);
    foreach my $fhit (@$fref)
    {
      my $lef = $fhit->[1];
      next if ($lef < $minbblen - $maxbblap);
      my $enz_f = $fhit->[0];
      my $zne_f = $seeks{$enz_f}->{$fhit->[2]};
      foreach my $thit (@$tref)
      {
        my $enz_t = $thit->[0];
        my $rig = $thit->[1] + $seeks{$enz_t}->{len};
        my $dist = $rig - $lef;
        next if ($dist < $minbblap || $dist > $maxbblap);
        next if ($rig > $chlen - ($minbblen - $maxbblap));
        my $zne_t = $seeks{$enz_t}->{$thit->[2]};
        my $ol = substr($chseq, $lef, $dist);
        my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
        $anns{$id} = {seq => $ol, left => $lef + 1, right => $rig,
                      enz_f => $enz_f, zne_f => $zne_f,
                      enz_t => $enz_t, zne_t => $zne_t};
        push @$laps, $lapobj;
        $id++;
      }
    }
    print "\t", scalar(@$laps), " overlaps found\n" if ($v);
    Bio::GeneDesign::Exception::UnBBable->throw("No overlaps passed filtering!")
      unless scalar(@$laps);

    #Filter overlaps for palindromes and mispriming
    #If vmatch or blast are not available, filter for uniqueness
    if ($GD->EMBOSS)
    {
      $laps = $GD->filter_palindromes(-sequences => $laps);
      print "\t", scalar(@$laps), " overlaps palindrome free\n" if ($v);
    }
    if ($GD->vmatch)
    {
      $laps = $GD->filter_vmatch(-sequences => $laps, -parent => $seqobj);
      print "\t", scalar(@$laps), " overlaps informative\n" if ($v);
    }
    elsif ($GD->BLAST)
    {
      $laps = $GD->filter_blast(-sequences => $laps, -parent => $seqobj);
      print "\t", scalar(@$laps), " overlaps informative\n" if ($v);
    }
    else
    {
      $laps = $GD->filter_uniqueness(-sequences => $laps);
      print "\t", scalar(@$laps), " overlaps unique\n" if ($v);
    }
    my @overlaps = map {$_->id} @$laps;
    
    
    my $f = floor $tarnum / 4;
    my $t = (4 - ($tarnum % 4)) % 4;
    $f = $f - ($t - 1) if ($t > 1);
    my $arr = 4 x $f . 3 x $t;
    my $frange = 5;
    my $trange = 4;
    
    
    #Pick overlaps properly spaced so as to obey length parameters.
    #Reroll if length parameters are violated by a choice.
    my $rbound = $maxbblen ? $chlen - $maxbblen : $chlen - $tarbblen;
    my ($lpos, $rpos) = (1, 1);
    my $lenz = undef;
    my %fails;
    my $fail = 0;
    my $redoflag = 0;
    my $counter = 1;
    my $amideadyet = 1;
    my @chosenlaps;
    while ($lpos < $rbound && $amideadyet < 500)
    {
      $amideadyet++;
      my $rtarget = ($lpos + $tarbblen);
      my @laps = grep {! exists $fails{$anns{$_}->{seq}}} @overlaps;
      if ($maxbblen)
      {
        @laps = grep {$anns{$_}->{right} <= $lpos + $maxbblen} @laps;
      }
      if ($minbblen)
      {
        @laps = grep {$anns{$_}->{right} >= $lpos + $minbblen} @laps;
      }
      if ($lenz)
      {
        @laps = grep {$RES->{$anns{$_}->{enz_t}}->temp == $RES->{$lenz}->temp} @laps;
        @laps = grep {$RES->{$anns{$_}->{enz_t}}->common_buffers($RES->{$lenz}, 1)} @laps;
      }
      if (defined $fpexcludes && scalar @chosenlaps < $frange)
      {
        @laps = grep {! exists $fpexcludes->{$anns{$_}->{enz_t}} } @laps;
      }
      if (defined $tpexcludes && scalar @chosenlaps > $trange)
      {
        @laps = grep {! exists $tpexcludes->{$anns{$_}->{enz_f}} } @laps;
      }
      @laps = sort { $RES->{$anns{$a}->{enz_f}}->score + $RES->{$anns{$a}->{enz_t}}->score
                 <=> $RES->{$anns{$b}->{enz_f}}->score + $RES->{$anns{$b}->{enz_t}}->score
                 || abs($anns{$a}->{right} - $rtarget) <=> abs($anns{$b}->{right} - $rtarget)
                 || length($anns{$a}->{left}) <=> length($anns{$b}->{left})}
              @laps;
      unless (scalar(@laps))
      {
        #discarding last choice
        my $disgrace = pop @chosenlaps;
        $fails{$anns{$disgrace}->{seq}}++;
        if (scalar @chosenlaps)
        {
          $rpos = $anns{$chosenlaps[-1]}->{right};
          $lpos = $anns{$chosenlaps[-1]}->{left};
          $lenz = $anns{$chosenlaps[-1]}->{enz_f};
        }
        else
        {
          ($rpos, $lpos) = (1, 1);
          $lenz = undef;
          @chosenlaps = ();
          %fails = ();
          print "\t\tAdjusting from $tarbblen " if ($v);
          $tarbblen = $counter % 2 == 0
                      ? $tarbblen + $counter
                      : $tarbblen - $counter;
          print "to $tarbblen " if ($v);
          if (($maxbblen && $tarbblen > $maxbblen)
           || ($minbblen && $tarbblen < $minbblen))
          {
            print "\n\t\tGDWARNING: giving up... OSCILLATED TOO FAR\n" if ($v);
            $fail++;
            Bio::GeneDesign::Exception::UnBBable->throw("no solution apparent");
          }
          $counter++;
        }
      }
      else
      {
        my $lap = $laps[0];
        my $laplen = length($anns{$lap}->{seq});
        $lenz =  $anns{$lap}->{enz_f};
        print "\tChoosing overlap $laplen bp ", $anns{$lap}->{seq}, if ($v);
        print " at ", $anns{$lap}->{left}, "..", $anns{$lap}->{right} if ($v);
        print " (", $anns{$lap}->{enz_f}, "..", $anns{$lap}->{enz_t}, ")" if ($v);
        print " (", $anns{$lap}->{zne_f}, "..", $anns{$lap}->{zne_t}, ")\n" if ($v);
        push @chosenlaps, $lap;
        $rpos = $anns{$chosenlaps[-1]}->{right};
        $lpos = $anns{$chosenlaps[-1]}->{left};
        $redoflag = 0;
        if ($lpos >= $rbound)
        {
          if (($minbblen && ($chlen - $lpos) + 1 < $minbblen)
           || ($maxbblen && ($chlen - $lpos) + 1 > $maxbblen))
          {
            print "\t\t\tdiscarding last choice\n" if ($v);
            my $disgrace = pop @chosenlaps;
            $fails{$anns{$disgrace}->{seq}}++;
            if (scalar @chosenlaps)
            {
              $rpos = $anns{$chosenlaps[-1]}->{right};
              $lpos = $anns{$chosenlaps[-1]}->{left};
              $lenz = $anns{$chosenlaps[-1]}->{enz_f};
            }
          }
        }
      }
      $lpos = $chlen if ($fail);
      if ($amideadyet >= 999)
      {
        $fail++;
        print "\n\t\tGDWARNING: Repeated search over 500 times!\n" if ($v);
        next;
      }
    }
    Bio::GeneDesign::Exception::UnBBable->throw("stuck in a rut!") if ($fail);
    
    
    #Make the building blocks
    $lpos = 1;
    my $bbnum = 1;
    my $lastenz = undef;
    my $lastmatch = undef;
    while (scalar @chosenlaps)
    {
      my $lap = shift @chosenlaps;
      my $bbseq = substr($chseq, $lpos - 1, $anns{$lap}->{right} - $lpos + 1);
      my $count = $GD->count($bbseq);
      my $bbname = $chname . "." . $GD->pad($bbnum, 2);
      my $atts = {name => $bbname, gcp => $count->{GCp}};
      $atts->{enz_f} = $lastenz if ($lastenz);
      $atts->{zne_f} = $lastmatch if ($lastmatch);
      $atts->{enz_t} = $anns{$lap}->{enz_t};
      $atts->{zne_t} = $anns{$lap}->{zne_t};
      $lastenz = $anns{$lap}->{enz_f};
      $lastmatch = $anns{$lap}->{zne_f};
      push @BBS, Bio::SeqFeature::Generic->new(
        -primary    => "building_block",
        -start      => $lpos,
        -end        => $anns{$lap}->{right},
        -tag        => $atts
      );
      $bbnum++;
      print "\t\t$bbname ", length($bbseq), " bp\n" if ($v);
      $lpos = $anns{$lap}->{left};
    }
    my $bbseq = substr($chseq, $lpos - 1);
    my $count = $GD->count($bbseq);
    my $bbname = $chname . "." . $GD->pad($bbnum, 2);
    my $atts = {name => $bbname, gcp => $count->{GCp},
                enz_f => $lastenz, zne_f => $lastmatch};
    push @BBS, Bio::SeqFeature::Generic->new(
      -primary    => "building_block",
      -start      => $lpos,
      -end        => $chlen,
      -tag        => $atts
    );
    $bbnum++;
    print "\t\t$bbname ", length($bbseq), " bp\n" if ($v);
  }
  else
  {
    my $bbseq = $chseq;
    my $count = $GD->count($bbseq);
    my $bbname = $chname . ".01";
    my $atts = {name => $bbname, gcp => $count->{GCp}};
    push @BBS, Bio::SeqFeature::Generic->new(
      -primary    => "building_block",
      -start      => 1,
      -end        => $chlen,
      -tag        => $atts
    );
  }
  
  print "\n\n" if ($v);
  return \@BBS;
}

1;

__END__

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
