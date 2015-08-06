=head1 NAME

Bio::GeneDesign::BuildingBlocks::overlap

=head1 VERSION

Version 5.00

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::GeneDesign::BuildingBlocks::overlap;
require Exporter;

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
  my ($GD, $seqobj, $tarbblen, $maxbblen, $bblenmin, $tarbblap, $stitch, $excludes, $fpexcludes, $tpexcludes, $v)
    = @_;

  my $chseq = $seqobj->seq;
  my $chlen = length($chseq);
  $tarbblen  = $chlen if ($chlen < $tarbblen);
  my @BBS;
  
  my $chname = $seqobj->id;
  my $len = length($seqobj->seq);
  
  
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


  my $laps = [];
  if ($tarnum > 1)
  {
    ## Search the construct sequence for sequences that can overlap BBs
    my $lef = $tarbblen - (3 * $tarbblap);
    my $rig = $lef + $tarbblap - 1;
    my $id = 1;
    while ($rig <= ($chlen - ($tarbblap + 1)))
    {
      my $ol = substr($chseq, $lef, $rig - $lef + 1);
      my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
      my $lefta   = Bio::Annotation::SimpleValue->new(-value => $lef + 1);
      $lapobj->add_Annotation('left', $lefta);
      my $righta  = Bio::Annotation::SimpleValue->new(-value => $rig + 1);
      $lapobj->add_Annotation('right', $righta);
      push @$laps, $lapobj;
      
      $id++;
      $lef++;
      $rig = $lef + $tarbblap - 1;
    }
    print "\t", scalar(@$laps), " overlaps found\n" if ($v);
    
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
    
    
    #Flatten overlap annotations for speed, figure out the average length of the
    #overlaps that passed filtering
    my %LAPS;
    foreach my $lap (@$laps)
    {
      my $lefann  = join(q{}, map {$_->value} $lap->get_Annotations('left'));
      my $rigann = join(q{}, map {$_->value} $lap->get_Annotations('right'));
    
      $LAPS{$lap->seq} = [$lap->seq, $lefann, $rigann];
    }
    my @overlaps = values %LAPS;
    
    
    
    #Pick overlaps properly spaced so as to obey length parameters.
    #Reroll if length parameters are violated by a choice.
    my $rbound = $maxbblen ? $chlen - $maxbblen : $chlen - $tarbblen;
    my ($lpos, $rpos) = (1, 1);
    my %fails;
    my $fail = 0;
    my $redoflag = 0;
    my $counter = 1;
    my $amideadyet = 1;
    my @chosenlaps;
    while ($lpos < $rbound && $amideadyet < 1000)
    {
      $amideadyet++;
      my $rtarget = ($lpos + $tarbblen);
      my @laps = grep {! exists $fails{$_->[0]}} @overlaps;
      if ($maxbblen)
      {
        @laps = grep {$_->[2] <= $lpos + $maxbblen} @laps;
      }
      if ($bblenmin)
      {
        @laps = grep {$_->[2] >= $lpos + $bblenmin} @laps;
      }
      @laps = sort {abs($a->[2] - $rtarget) <=> abs($b->[2] - $rtarget)
                || length($a->[0]) <=> length($b->[0])}
              @laps;
      unless (scalar(@laps))
      {
        #discarding last choice
        my $disgrace = pop @chosenlaps;
        $fails{$disgrace->[0]}++;
        if (scalar @chosenlaps)
        {
          $rpos = $chosenlaps[-1]->[2];
          $lpos = $chosenlaps[-1]->[1];
        }
        else
        {
          ($rpos, $lpos) = (1, 1);
          @chosenlaps = ();
          %fails = ();
          print "\t\tAdjusting from $tarbblen " if ($v);
          $tarbblen = $counter % 2 == 0
                      ? $tarbblen + $counter
                      : $tarbblen - $counter;
          print "to $tarbblen " if ($v);
          if (($maxbblen && $tarbblen > $maxbblen)
           || ($bblenmin && $tarbblen < $bblenmin))
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
        print "\tChoosing overlap $lap->[0] at $lap->[1], $lap->[2]\n" if ($v);
        push @chosenlaps, $lap;
        $rpos = $chosenlaps[-1]->[2];
        $lpos = $chosenlaps[-1]->[1];
        $redoflag = 0;
        if ($lpos >= $rbound)
        {
          if (($bblenmin && ($chlen - $lpos) + 1 < $bblenmin)
           || ($maxbblen && ($chlen - $lpos) + 1 > $maxbblen))
          {
            print "\t\t\tdiscarding last choice\n" if ($v);
            my $disgrace = pop @chosenlaps;
            $fails{$disgrace->[0]}++;
            if (scalar @chosenlaps)
            {
              $rpos = $chosenlaps[-1]->[2];
              $lpos = $chosenlaps[-1]->[1];
            }
          }
        }
      }
      $lpos = $chlen if ($fail);
      if ($amideadyet >= 999)
      {
        $fail++;
        print "GDWARNING: Repeated search over 1000 times!\n" if ($v);
        next;
      }
    }
    Bio::GeneDesign::Exception::UnBBable->throw("stuck in a rut!") if ($fail);
    
    
    #Make the building blocks
    $lpos = 1;
    my $bbnum = 1;
    while (scalar @chosenlaps)
    {
      my $lap = shift @chosenlaps;
      my $bbseq = substr($chseq, $lpos - 1, $lap->[2] - $lpos + 1);
      my $count = $GD->count($bbseq);
      my $bbname = $chname . "." . $GD->pad($bbnum, 2);
      push @BBS, Bio::SeqFeature::Generic->new(
        -primary    => "building_block",
        -start      => $lpos,
        -end        => $lap->[2],
        -tag        => {
          #sequence   => $bbseq,
          name       => $bbname,
          gcp        => $count->{'GCp'}
        }
      );
      $bbnum++;
      print "\t\t$bbname ", length($bbseq), " bp\n" if ($v);
      $lpos = $lap->[1];
    }
    my $bbseq = substr($chseq, $lpos - 1);
    my $count = $GD->count($bbseq);
    my $bbname = $chname . "." . $GD->pad($bbnum, 2);
    push @BBS, Bio::SeqFeature::Generic->new(
      -primary    => "building_block",
      -start      => $lpos,
      -end        => $chlen,
      -tag        => {
        #sequence   => $bbseq,
        name       => $bbname,
        gcp        => $count->{'GCp'}
      }
    );
    $bbnum++;
    print "\t\t$bbname ", length($bbseq), " bp\n" if ($v);
  }
  else
  {
    my $bbseq = $chseq;
    my $count = $GD->count($bbseq);
    my $bbname = $chname . ".01";
    push @BBS, Bio::SeqFeature::Generic->new(
      -primary    => "building_block",
      -start      => 1,
      -end        => $chlen,
      -tag        => {
        #sequence   => $bbseq,
        name       => $bbname,
        gcp        => $count->{'GCp'}
      }
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
