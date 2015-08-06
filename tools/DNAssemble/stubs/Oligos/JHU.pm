
=head1 NAME

Bio::GeneDesign::Oligos::JHU

=head1 VERSION

Version 5.00

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::GeneDesign::Oligos::JHU;
require Exporter;

use Bio::GeneDesign::Basic qw(_complement _melt);
use Bio::GeneDesign::Exceptions;
use Bio::Seq;
use Bio::SeqFeature::Generic;

use strict;
use warnings;

our $VERSION = 5.00;

use base qw(Exporter);
our @EXPORT_OK = qw(
  _chop_oligos
);
our %EXPORT_TAGS =  ( GD => \@EXPORT_OK);
  
=head2 _carve_building_blocks()

  
=cut

sub _chop_oligos
{
  my ($GD, $bb, $olilenmax, $olilenmin, $olilen, $laptemp, $laplenmin,
      $tmtol, $poolsize, $poolnummax, $v) = @_;

  my @OLIGOS;
  my $fail = 0;
  my $univF = 0;
  my $univR = 0;
    
  #Add prefixes and/or suffixes, if necessary; note whether universals are
  #indicated. create a seqobj for BB seq for filtering overlaps later
  my $bbname = join(q{}, $bb->get_tag_values("name"));
  my $bbseq = $bb->seq->seq;
  if ($bb->has_tag("prefix"))
  {
    $bbseq = join(q{}, $bb->get_tag_values("prefix")) . $bbseq;
    $univF = 1;
  }
  if ($bb->has_tag("suffix"))
  {
    $bbseq = $bbseq . join(q{}, $bb->get_tag_values("suffix"));
    $univR = 1;
  }
  $bbseq        =~ s/\s//xg;
  $bbseq        = uc $bbseq;
  my $bblen     = length($bbseq);
  my $bbstart   = $bb->start;
  my $bbend     = $bb->end;
  my $bbseqobj  = Bio::Seq->new(-id => $bbname, -seq => $bbseq);
  print "Chopping $bbname ($bblen bp)\n" if ($v);

    
  #Find all of the overlaps in the building block that fit the Tm profile
  #Store the Tm and position with the overlap sequence as a Bio::Seq object
  my $laps = [];
  my $lef = 0;
  my $id = 0;
  my $rig = $laplenmin || 7;
  while ($rig <= $bblen)
  {
    my $ol = substr($bbseq, $lef, $rig - $lef + 1);
    my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
    my $Tm = $GD->melt(-sequence => $lapobj);
    
    if ($Tm < $laptemp - $tmtol)
    {
      $rig++;
      next;
    }
    elsif ($Tm > $laptemp + $tmtol)
    {
      $lef++;
      $rig = $laplenmin  ? $lef + $laplenmin  : $lef + 7;
      next;
    }
    
    my $TmA     = Bio::Annotation::SimpleValue->new(-value => $Tm);
    $lapobj->add_Annotation('Tm', $TmA);
    my $lefta   = Bio::Annotation::SimpleValue->new(-value => $lef + 1);
    $lapobj->add_Annotation('left', $lefta);
    my $righta  = Bio::Annotation::SimpleValue->new(-value => $rig + 1);
    $lapobj->add_Annotation('right', $righta);
    push @$laps, $lapobj;

    $lef++;
    $id++;
    $rig = $laplenmin  ? $lef + $laplenmin  : $lef + 7;
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
    $laps = $GD->filter_vmatch(-sequences => $laps, -parent => $bbseqobj);
    print "\t", scalar(@$laps), " overlaps informative\n" if ($v);
  }
  elsif ($GD->BLAST)
  {
    $laps = $GD->filter_blast(-sequences => $laps, -parent => $bbseqobj);
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
    my $Tm    = join(q{}, map {$_->value} $lap->get_Annotations('Tm'));
    my $lefann  = join(q{}, map {$_->value} $lap->get_Annotations('left'));
    my $rigann = join(q{}, map {$_->value} $lap->get_Annotations('right'));
    
    $LAPS{$lap->seq} = [$Tm, $lap->seq, $lefann, $rigann];
  }
  my @overlaps = values %LAPS;
  my $avg = 0;
  $avg += length($_->[1]) foreach @overlaps;
  $avg = sprintf "%.0f", $avg / scalar(@overlaps);


  #Decide what size oligos to go for and how many.
  #Adjust the target size to avoid length outliers.
  my $tarlen = $olilen || sprintf "%.0f", ($olilenmin + $olilenmax) / 2;
  my $tarnum = sprintf "%.0f", $bblen / ($tarlen - $avg);
  $tarnum++ unless ($tarnum % 2 == 0);
  my $tarolilen = sprintf "%.0f", $bblen / $tarnum;
  my $diff = $bblen - (($tarnum * $tarolilen) - ($avg * ($tarnum - 1)));
  print "\ttarget: $tarnum oligos size $tarolilen bp rem $diff\n" if ($v);
  if (abs($diff) >= $tarnum)
  {
    my $rem = sprintf "%.0f", $diff / $tarnum;
    $tarolilen = $tarolilen + $rem;
    $diff = $diff - ($tarnum * $rem);
  }
  print "\t final: $tarnum oligos size $tarolilen bp rem $diff\n" if ($v);
    
  #Pick overlaps properly spaced so as to obey length parameters.
  #Reroll if length parameters are violated by a choice.
  my $rbound = $olilenmax ? $bblen - $olilenmax : $bblen - $tarolilen;
  my ($lpos, $rpos) = (1, 1);
  my %fails;
  my $redoflag = 0;
  my $counter = 1;
  my $amideadyet = 1;
  my $pstrand = 1;
  my @chosenlaps;
  while ($lpos < $rbound && $amideadyet < 1000)
  {
    $amideadyet++;
    my $rtarget = ($lpos + $tarolilen);
    my @laps = grep {! exists $fails{$_->[1]}} @overlaps;
    if ($olilenmax)
    {
      @laps = grep {$_->[3] <= $lpos + $olilenmax} @laps;
    }
    if ($olilenmin)
    {
      @laps = grep {$_->[3] >= $lpos + $olilenmin} @laps;
    }
    @laps = sort {abs($a->[3] - $rtarget) <=> abs($b->[3] - $rtarget)
              || length($a->[1]) <=> length($b->[1])}
            @laps;
    unless (scalar(@laps))
    {
      #discarding last choice
      my $disgrace = pop @chosenlaps;
      $pstrand = $pstrand == 1  ? -1  : 1;
      $fails{$disgrace->[1]}++;
      if (scalar @chosenlaps)
      {
        $rpos = $chosenlaps[-1]->[3];
        $lpos = $chosenlaps[-1]->[2];
      }
      else
      {
        ($rpos, $lpos) = (1, 1);
        @chosenlaps = ();
        %fails = ();
        print "\t\tAdjusting from $tarolilen " if ($v);
        $tarolilen = $counter % 2 == 0
                    ? $tarolilen + $counter
                    : $tarolilen - $counter;
        print "to $tarolilen " if ($v);
        if (($olilenmax && $tarolilen > $olilenmax)
         || ($olilenmin && $tarolilen < $olilenmin))
        {
          print "GDWARNING: giving up... OSCILLATED TOO FAR\n" if ($v);
          $fail++;
          Bio::GeneDesign::Exception::UnOLable->throw("no solution apparent");
        }
        $counter++;
      }
    }
    else
    {
      my $lap = $laps[0];
      my $lastl = scalar @chosenlaps  ? $chosenlaps[-1]->[2]  : 1;
      print "\tChoosing $lap->[1] at $lap->[2], $lap->[3], $pstrand, $olilen\n" if ($v);
      push @chosenlaps, $lap;
      $pstrand = $pstrand == 1  ? -1  : 1;
      $rpos = $chosenlaps[-1]->[3];
      $lpos = $chosenlaps[-1]->[2];
      $redoflag = 0;
      if ($lpos >= $rbound)
      {
        print "\t\tChecking if this is the end\n" if ($v);
        my $olilen = $rpos - $lastl + 1;
        my $dist = $bblen - $lpos + 1;
        if (($olilenmin && $dist < $olilenmin)
         || ($olilenmax && $dist > $olilenmax)
         || ($pstrand == 1 && $dist + $olilen > $olilenmax))
        {
          print "\t\t\tdiscarding last choice\n" if ($v);
          $pstrand = $pstrand == 1  ? -1  : 1;
          my $disgrace = pop @chosenlaps;
          $fails{$disgrace->[1]}++;
          if (scalar @chosenlaps)
          {
            $rpos = $chosenlaps[-1]->[3];
            $lpos = $chosenlaps[-1]->[2];
          }
        }
      }
    }
    $lpos = $bblen if ($fail);
    if ($amideadyet >= 999)
    {
      $fail++;
      print "GDWARNING: Repeated search over 1000 times!\n" if ($v);
      next;
    }
  }
  Bio::GeneDesign::Exception::UnOLable->throw("stuck in a rut!") if ($fail);


  #Make the oligos
  #All the oligos
  $lpos = 1;
  my $olnum = 1;
  my $strand = 1;
  while (scalar @chosenlaps)
  {
    my $lap = shift @chosenlaps;
    my $olseq = substr($bbseq, $lpos-1, $lap->[3] - $lpos + 1);
    $olseq = _complement($olseq, 1) if ($strand == -1);
    my $olname = $bbname . "." . $GD->pad($olnum, 2);
    push @OLIGOS, Bio::SeqFeature::Generic->new(
      -primary    => "assembly_oligo",
      -start      => $lpos,
      -end        => $lap->[3],
      -strand     => $strand,
      -tag        => {
        sequence   => $olseq,
        name       => $olname,
        bb         => $bbname,
      }
    );
    $olnum++;
    print "\t\t$olname $strand ", length($olseq), " bp\n" if ($v);
    $strand = $strand == 1 ? -1  : 1;
    $lpos = $lap->[2];
  }
  my $olseq = substr($bbseq, $lpos-1);
  $olseq = _complement($olseq, 1) if ($strand == -1);
  my $olname = $bbname . "." . $GD->pad($olnum, 2);
  push @OLIGOS, Bio::SeqFeature::Generic->new(
    -primary    => "assembly_oligo",
    -start      => $lpos,
    -end        => $bblen,
    -strand     => $strand,
    -tag        => {
      sequence   => $olseq,
      name       => $olname,
      bb         => $bbname,
    }
  );
  $olnum++;
  print "\t\t$olname $strand ", length($olseq), " bp\n" if ($v);
  
    
  #Adjust the start and stop coordinates of the oligos to reflect the building
  #blocks placement in its parent sequence and the non-parental prefix/suffixes
  foreach my $ol (@OLIGOS)
  {
    my ($start, $end) = ($ol->start, $ol->end);
    if ($bbstart != 1)
    {
      $start  += $bbstart;
      $end    += $bbstart;
    }
    if ($bb->has_tag("prefix"))
    {
      my $offset = length(join(q{}, $bb->get_tag_values("prefix")));
      $start  -= $offset;
      $end    -= $offset;
    }
    $ol->start($start);
    $ol->end($end);
  }
  
  print "\n\n" if ($v);
  return \@OLIGOS;
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
