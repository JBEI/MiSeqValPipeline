=head1 NAME

Bio::GeneDesign::Oligos::JGI

=head1 VERSION

Version 5.10

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::GeneDesign::Oligos::JGI;
require Exporter;

use Bio::GeneDesign::Basic qw(_complement _melt);
use Bio::GeneDesign::Exceptions;
use Bio::Seq;
use Bio::SeqFeature::Generic;

use strict;
use warnings;

our $VERSION = 5.10;

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

  # Add prefixes and/or suffixes, if necessary; create a seqobj for BB seq for
  # filtering overlaps later
  #
  my $bbname = join q{}, $bb->get_tag_values('name');
  my $bbseq = $bb->seq->seq;
  if ($bb->has_tag('prefix'))
  {
    $bbseq = join(q{}, $bb->get_tag_values('prefix')) . $bbseq;
  }
  if ($bb->has_tag('suffix'))
  {
    $bbseq = $bbseq . join q{}, $bb->get_tag_values('suffix');
  }
  $bbseq        =~ s/\s//msixg;
  $bbseq        = uc $bbseq;
  my $bblen     = length $bbseq;
  my @bbmask    = (0) x $bblen;
  my $bbstart   = $bb->start;
  my $bbend     = $bb->end;
  my $bbseqobj  = Bio::Seq->new(-id => $bbname, -seq => $bbseq);
  print "Chopping $bbname ($bblen bp)\n" if ($v);


  # Find all of the overlaps in the building block that fit the Tm profile
  # Store the Tm and position with the overlap sequence as a Bio::Seq object
  #
  my $laps = [];
  my $lef = 0;
  my $id = 0;
  my $rig = $laplenmin || 7;
  while ($rig <= $bblen)
  {
    my $ol = substr $bbseq, $lef, $rig - $lef + 1;
    my $Tm = $GD->melt(-sequence => $ol);

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

    my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
    my $TmA     = Bio::Annotation::SimpleValue->new(-value => $Tm);
    $lapobj->add_Annotation('Tm', $TmA);
    my $lefta   = Bio::Annotation::SimpleValue->new(-value => $lef + 1);
    $lapobj->add_Annotation('left', $lefta);
    my $righta  = Bio::Annotation::SimpleValue->new(-value => $rig + 1);
    $lapobj->add_Annotation('right', $righta);
    push @{$laps}, $lapobj;

    $lef++;
    $id++;
    $rig = $laplenmin  ? $lef + $laplenmin  : $lef + 7;
  }
  print "\t", scalar @{$laps}, " overlaps found\n" if ($v);


  #Filter overlaps for palindromes and mispriming
  #If vmatch or blast are not available, filter for uniqueness
  if ($GD->EMBOSS)
  {
    $laps = $GD->filter_palindromes(-sequences => $laps);
    print "\t", scalar @{$laps}, " overlaps palindrome free\n" if ($v);
  }
  if ($GD->vmatch)
  {
    $laps = $GD->filter_vmatch(-sequences => $laps, -parent => $bbseqobj);
    print "\t", scalar @{$laps}, " overlaps informative\n" if ($v);
  }
  elsif ($GD->BLAST)
  {
    $laps = $GD->filter_blast(-sequences => $laps, -parent => $bbseqobj);
    print "\t", scalar @{$laps}, " overlaps informative\n" if ($v);
  }
  else
  {
    $laps = $GD->filter_uniqueness(-sequences => $laps);
    print "\t", scalar @{$laps}, " overlaps unique\n" if ($v);
  }


  #Flatten overlap annotations for speed, figure out the average length of the
  #overlaps that passed filtering
  my %LAPS;
  foreach my $lap (@{$laps})
  {
    my $Tm     = join q{}, map {$_->value} $lap->get_Annotations('Tm');
    my $lefann = join q{}, map {$_->value} $lap->get_Annotations('left');
    my $rigann = join q{}, map {$_->value} $lap->get_Annotations('right');

    $LAPS{$lap->seq} = [$Tm, $lap->seq, $lefann, $rigann];
    for my $index ($lefann .. $rigann)
    {
      $bbmask[$index - 1]++;
    }
  }
  my @overlaps = values %LAPS;
  my $avg = 0;
  $avg += length $_->[1] foreach @overlaps;
  $avg = sprintf '%.0f', $avg / (scalar @overlaps);


  #Assess spread and decide if this is possible
  #Bio::GeneDesign::Exception::UnOLable->throw('no solution apparent');

  #Decide what size oligos to go for and how many.
  #Adjust the target size to avoid length outliers.
  my $tarlen = $olilen || sprintf '%.0f', ($olilenmin + $olilenmax) / 2;
  my $tarnum = sprintf '%.0f', $bblen / ($tarlen - $avg);
  my $tarolilen = sprintf '%.0f', $bblen / $tarnum;
  my $diff = $bblen - (($tarnum * $tarolilen) - ($avg * ($tarnum - 1)));
  print "\ttarget: $tarnum oligos size $tarolilen bp rem $diff\n" if ($v);
  if (abs($diff) >= $tarnum)
  {
    my $rem = sprintf '%.0f', $diff / $tarnum;
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
  my @chosenlaps;
  while ($lpos < $rbound && $amideadyet < 500)
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
              || length $a->[1] <=> length $b->[1] }
            @laps;
    if (! scalar @laps)
    {
      #discarding last choice
      my $disgrace = pop @chosenlaps;
      $fails{$disgrace->[1]}++ if (defined $disgrace);
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
          Bio::GeneDesign::Exception::UnOLable->throw('no solution apparent');
        }
        $counter++;
      }
    }
    else
    {
      my $lap = $laps[0];
      print "\tChoosing overlap $lap->[1] at $lap->[2], $lap->[3]\n" if ($v);
      push @chosenlaps, $lap;
      $rpos = $chosenlaps[-1]->[3];
      $lpos = $chosenlaps[-1]->[2];
      $redoflag = 0;
      if ($lpos >= $rbound)
      {
        if (($olilenmin && ($bblen - $lpos) + 1 < $olilenmin)
         || ($olilenmax && ($bblen - $lpos) + 1 > $olilenmax))
        {
          print "\t\t\tdiscarding last choice\n" if ($v);
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
      print "GDWARNING: Repeated search over 500 times!\n" if ($v);
      next;
    }
  }
  Bio::GeneDesign::Exception::UnOLable->throw('stuck in a rut!') if ($fail);


  #Make the oligos
  my $y = 1;
  $lpos = 1;
  my $divisor = $poolnummax || 2;
  my $count = scalar @chosenlaps;
  my $maxol = $poolnummax  ? $poolnummax * ($poolsize - 1) : undef;
  if ($maxol && ( $count + 1 ) > $maxol)
  {
    print "GDWARNING: $bbname has EXCEEDS MAX POOL NUMBER.\n";
    $divisor++;
  }
  my $flag = $poolsize  ? 1 : 0;
  my $olcount = $count + 1;
  $olcount++ if ($olcount % 2 != 0);
  my $div;
  if ($flag)
  {
    $div = ($count + 1) <= $poolsize - 1
          ?  $poolsize
          : sprintf '%.0f', $olcount / $divisor;
  }
  my $pcanum = 1;
  my $olnum = 1;

  ##Make the small forward
  my $frig = sprintf '%.0f', $avg / 2;
  $frig++ while (_melt((substr $bbseq, 0, $frig)) < $laptemp - $tmtol);
  my $olseq = substr $bbseq, 0, $frig;
  my $olname = $bbname . q{.} . $GD->pad($olnum, 2);
  push @OLIGOS, Bio::SeqFeature::Generic->new(
    -primary    => 'assembly_oligo',
    -start      => 1,
    -end        => $frig + 1,
    -tag        => {
      sequence   => $olseq,
      name       => $olname,
      parent     => $bbname,
      pool       => $pcanum,
    }
  );
  $olnum++;
  print "\t\t$olname ", length $olseq, " bp\n" if ($v);

  while (scalar @chosenlaps)
  {
    my $lap = shift @chosenlaps;
    $olseq = substr $bbseq, $lpos-1, $lap->[3] - $lpos + 1;
    $olname = $bbname . q{.} . $GD->pad($olnum, 2);
    push @OLIGOS, Bio::SeqFeature::Generic->new(
      -primary    => 'assembly_oligo',
      -start      => $lpos,
      -end        => $lap->[3],
      -tag        => {
        sequence   => $olseq,
        name       => $olname,
        parent     => $bbname,
        pool       => $pcanum
      }
    );
    $olnum++;
    print "\t\t$olname ", length $olseq, " bp\n" if ($v);

    #Make the pool spanners if necessary
    #Increment the poolnumber
    if ($flag && $y % $div == 0)
    {
      $olname = $bbname . q{.} . $GD->pad($olnum, 2);
      push @OLIGOS, Bio::SeqFeature::Generic->new(
        -primary    => 'assembly_oligo',
        -start      => $lap->[2],
        -end        => $lap->[3],
        -strand     => -1,
        -tag        => {
          sequence   => _complement($lap->[1], 1),
          name       => $olname,
          parent     => $bbname,
          pool       => $pcanum
        }
      );
      print "\t\t$olname ", length $lap->[1], " bp\n" if ($v);
      $pcanum++;
      $olnum++;
      $olname = $bbname . q{.} . $GD->pad($olnum, 2);
      push @OLIGOS, Bio::SeqFeature::Generic->new(
        -primary    => 'assembly_oligo',
        -start      => $lap->[2],
        -end        => $lap->[3],
        -tag        => {
          sequence   => $lap->[1],
          name       => $olname,
          parent     => $bbname,
          pool       => $pcanum
        }
      );
      print "\t\t$olname ", length $lap->[1], " bp\n" if ($v);
      $olnum++;
    }

    $lpos = $lap->[2];
    $y++;
  }

  #Make the lagging oligo
  $olseq = substr $bbseq, $lpos - 1;
  $olname = $bbname . q{.} . $GD->pad($olnum, 2);
  push @OLIGOS, Bio::SeqFeature::Generic->new(
    -primary    => 'assembly_oligo',
    -start      => $lpos,
    -end        => $bblen,
    -tag        => {
      sequence   => $olseq,
      name       => $olname,
      parent     => $bbname,
      pool       => $pcanum
    }
  );
  $olnum++;
  $bb->add_tag_value('pca_pool_count', $pcanum);
  print "\t\t$olname ", length $olseq, " bp\n" if ($v);

  ##Make the small reverse
  my $flef = sprintf '%.0f', $avg / 2;
  $flef++ while (_melt((substr $bbseq, -$flef)) < $laptemp - $tmtol);
  $olseq = substr $bbseq, -$flef;
  $olseq = $GD->complement($olseq, 1);
  $olname = $bbname . q{.} . $GD->pad($olnum, 2);
  push @OLIGOS, Bio::SeqFeature::Generic->new(
    -primary    => 'assembly_oligo',
    -start      => $bblen - $flef,
    -strand     => -1,
    -end        => $bblen,
    -tag        => {
      sequence   => $olseq,
      name       => $olname,
      parent     => $bbname,
      pool       => $pcanum,
    }
  );
  $olnum++;
  print "\t\t$olname ", length $olseq, " bp\n" if ($v);

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
    if ($bb->has_tag('prefix'))
    {
      my $offset = length join q{}, $bb->get_tag_values('prefix');
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

Copyright (c) 2013, GeneDesign developers
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
