#!/usr/bin/env perl

use Bio::DNAssemble;
use Bio::BioStudio::Mask;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use POSIX;
use Readonly;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_155_CorrectGBlocks_' . $VERSION;
my $DNAS = '_155';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
  'input=s'              => \$p{INPUT},
  'output:s'             => \$p{OUTPUT},
  'logfile:s'            => \$p{LOGPATH},
  'SKIP_WARNINGS'        => \$p{SKIP_WARNINGS},
  'GC_GLOBAL_MIN=i'      => \$p{GC_GLOBAL_MIN},
  'GC_GLOBAL_MAX=i'      => \$p{GC_GLOBAL_MAX},
  'GC_REGIONAL_WINDOW'   => \$p{GC_REGIONAL_WINDOW},
  'GC_REGIONAL_MIN=i'    => \$p{GC_REGIONAL_MIN},
  'GC_REGIONAL_MAX=i'    => \$p{GC_REGIONAL_MAX},
  'GC_LOCAL_WINDOW'      => \$p{GC_LOCAL_WINDOW},
  'GC_LOCAL_MIN=i'       => \$p{GC_LOCAL_MIN},
  'GC_LOCAL_MAX=i'       => \$p{GC_LOCAL_MAX},
  'GC_TERMINAL_WINDOW'   => \$p{GC_TERMINAL_WINDOW},
  'GC_TERMINAL_MIN=i'    => \$p{GC_TERMINAL_MIN},
  'GC_TERMINAL_MAX=i'    => \$p{GC_TERMINAL_MAX},
  'HP_A_LENGTH=i'        => \$p{HP_A_LENGTH},
  'HP_T_LENGTH=i'        => \$p{HP_T_LENGTH},
  'HP_C_LENGTH=i'        => \$p{HP_C_LENGTH},
  'HP_G_LENGTH=i'        => \$p{HP_G_LENGTH},
  'DIMER_LIMIT=i'        => \$p{DIMER_LIMIT},
  'TRIMER_LIMIT=i'       => \$p{TRIMER_LIMIT},
  'DR_LOCAL_WINDOW=i'    => \$p{DR_LOCAL_WINDOW},
  'DR_GLOBAL_RATIO=i'    => \$p{DR_GLOBAL_RATIO},
  'DR_SOLO_RATIO=i'      => \$p{DR_SOLO_RATIO},
  'DR_LENGTH=i'          => \$p{DR_LENGTH},
  'IR_LENGTH_MIN=i'      => \$p{IR_LENGTH_MIN},
  'IR_LENGTH_MAX=i'      => \$p{IR_LENGTH_MAX},
  'IR_WIDTH=i'           => \$p{IR_WIDTH},
  'help'                 => \$p{HELP}
);


################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>'NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE')
  if ($p{HELP});

my $DNA = Bio::DNAssemble->new();

my $OUTH;
if ($p{LOGPATH})
{
  open $OUTH, '>>', $p{LOGPATH};
}
else
{
  $OUTH = *STDOUT;
}
print {$OUTH} "\n\n******* $DNAV WORKING\n\n";

#The input file must exist and be a format we care to read.
my $design = $DNA->load_design_from_file($p{INPUT});
my $GD = $DNA->GD();
my $CT = $GD->codontable();
my $RCT = $GD->reversecodontable();
my $filename = $design->filename();
$filename = $1 if ($filename =~ m{(.+)\_}msix);
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.xml};

##DEBUG
$p{ERROR_NUMBER} = 0;
$p{SKIP_WARNINGS} = $p{SKIP_WARNINGS} || 1;

Readonly my $D_GBLOCK_MAX => 1000;
$p{GBLOCK_MAX} = $p{GBLOCK_MAX} || $D_GBLOCK_MAX;

## GLOBAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_GLOBAL_MIN => 25;
Readonly my $D_GC_GLOBAL_MAX => 69;
$p{GC_GLOBAL_MIN} = $p{GC_GLOBAL_MIN} || $D_GC_GLOBAL_MIN;
$p{GC_GLOBAL_MAX} = $p{GC_GLOBAL_MAX} || $D_GC_GLOBAL_MAX;
die 'BS_ERROR: Global GC percentage arguments do not parse.'
  if ($p{GC_GLOBAL_MIN} >= $p{GC_GLOBAL_MAX});

## REGIONAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_REGIONAL_WINDOW => 100;
$p{GC_REGIONAL_WINDOW} = $p{GC_REGIONAL_WINDOW} || $D_GC_REGIONAL_WINDOW;
Readonly my $D_GC_REGIONAL_MIN => 28;
Readonly my $D_GC_REGIONAL_MAX => 77;
$p{GC_REGIONAL_MIN} = $p{GC_REGIONAL_MIN} || $D_GC_REGIONAL_MIN;
$p{GC_REGIONAL_MAX} = $p{GC_REGIONAL_MAX} || $D_GC_REGIONAL_MAX;
die 'BS_ERROR: Regional GC percentage arguments do not parse.'
  if ($p{GC_REGIONAL_MIN} >= $p{GC_REGIONAL_MAX});

## LOCAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_LOCAL_WINDOW => 20;
$p{GC_LOCAL_WINDOW} = $p{GC_LOCAL_WINDOW} || $D_GC_LOCAL_WINDOW;
Readonly my $D_GC_LOCAL_MIN => 15;
Readonly my $D_GC_LOCAL_MAX => 90;
$p{GC_LOCAL_MIN} = $p{GC_LOCAL_MIN} || $D_GC_LOCAL_MIN;
$p{GC_LOCAL_MAX} = $p{GC_LOCAL_MAX} || $D_GC_LOCAL_MAX;
die 'BS_ERROR: Local GC percentage arguments do not parse.'
  if ($p{GC_LOCAL_MIN} >= $p{GC_LOCAL_MAX});
die 'BS_ERROR: Regional and local window sizes do not parse.'
  if ($p{GC_LOCAL_WINDOW} >= $p{GC_REGIONAL_WINDOW});

## TERMINAL GC PERCENTAGE ARGUMENTS
Readonly my $D_GC_TERMINAL_WINDOW => 30;
$p{GC_TERMINAL_WINDOW} = $p{GC_TERMINAL_WINDOW} || $D_GC_TERMINAL_WINDOW;
Readonly my $D_GC_TERMINAL_MIN => 24;
Readonly my $D_GC_TERMINAL_MAX => 76;
$p{GC_TERMINAL_MIN} = $p{GC_TERMINAL_MIN} || $D_GC_TERMINAL_MIN;
$p{GC_TERMINAL_MAX} = $p{GC_TERMINAL_MAX} || $D_GC_TERMINAL_MAX;
die 'BS_ERROR: Terminal GC percentage arguments do not parse.'
  if ($p{GC_TERMINAL_MIN} >= $p{GC_TERMINAL_MAX});

## HOMOPOLYMER LENGTHS ARGUMENTS
Readonly my $HP_MIN_LENGTH => 5;
Readonly my $D_HP_A_LENGTH => 12;
$p{HP_A_LENGTH} = $p{HP_A_LENGTH} || $D_HP_A_LENGTH;
die 'BS_ERROR: The homopolymer limit for A is too low.'
  if ($p{HP_A_LENGTH} < $HP_MIN_LENGTH);
Readonly my $D_HP_T_LENGTH => 12;
$p{HP_T_LENGTH} = $p{HP_T_LENGTH} || $D_HP_T_LENGTH;
die 'BS_ERROR: The homopolymer limit for T is too low.'
  if ($p{HP_T_LENGTH} < $HP_MIN_LENGTH);
Readonly my $D_HP_C_LENGTH => 8;
$p{HP_C_LENGTH} = $p{HP_C_LENGTH} || $D_HP_C_LENGTH;
die 'BS_ERROR: The homopolymer limit for C is too low.'
  if ($p{HP_C_LENGTH} < $HP_MIN_LENGTH);
Readonly my $D_HP_G_LENGTH => 8;
$p{HP_G_LENGTH} = $p{HP_G_LENGTH} || $D_HP_G_LENGTH;
die 'BS_ERROR: The homopolymer limit for G is too low.'
  if ($p{HP_G_LENGTH} < $HP_MIN_LENGTH);

## REPEAT ARGUMENTS
Readonly my $D_DIMER_LIMIT => 10;
$p{DIMER_LIMIT} = $p{DIMER_LIMIT} || $D_DIMER_LIMIT;
Readonly my $D_TRIMER_LIMIT => 6;
$p{TRIMER_LIMIT} = $p{TRIMER_LIMIT} || $D_TRIMER_LIMIT;

##LOCAL REPEAT ARGUMENTS
Readonly my $D_DR_LOCAL_WINDOW => 70;
$p{DR_LOCAL_WINDOW} = $p{DR_LOCAL_WINDOW} || $D_DR_LOCAL_WINDOW;
Readonly my $D_DR_LOCAL_RATIO => 90;
$p{DR_LOCAL_RATIO} = $p{DR_LOCAL_RATIO} || $D_DR_LOCAL_RATIO;

Readonly my $D_DR_GLOBAL_RATIO => 69;
$p{DR_GLOBAL_RATIO} = $p{DR_GLOBAL_RATIO} || $D_DR_GLOBAL_RATIO;
Readonly my $D_DR_SOLO_RATIO => .4;
$p{DR_SOLO_RATIO} = $p{DR_SOLO_RATIO} || $D_DR_SOLO_RATIO;
Readonly my $D_DR_LENGTH => 8;
$p{DR_LENGTH} = $p{DR_LENGTH} || $D_DR_LENGTH;
Readonly my $D_IR_LENGTH_MIN => 16;
$p{IR_LENGTH_MIN} = $p{IR_LENGTH_MIN} || $D_IR_LENGTH_MIN;
Readonly my $D_IR_LENGTH_MAX => 19;
$p{IR_LENGTH_MAX} = $p{IR_LENGTH_MAX} || $D_IR_LENGTH_MAX;
Readonly my $D_IR_WIDTH => 100;
$p{IR_WIDTH} = $p{IR_WIDTH} || $D_IR_WIDTH;

################################################################################
################################## INSPECTING ##################################
################################################################################
my @constructs = $design->get_constructs(-kind => 'unknown');
my @gblocks = $design->get_constructs(-kind => 'building_block');
my %ghsh = map {$_->id => $_} @gblocks;
my @noissues;
my @someissues;
my $badgblockcount = 0;
my $goodgblockcount = 0;
my $goodconstructcount = 0;
my $badconstructcount = 0;

my ($make, $nomake) = (0, 0);
foreach my $construct (@constructs)
{
  my @features = $construct->get_features();
  print "Working on ", $construct->id, "\n";

  my @subids = $construct->get_subconstruct_ids();
  my @theseblocks = map {$ghsh{$_}} grep {exists $ghsh{$_}} @subids;
  my @myissues = ();
  foreach my $gblock (@theseblocks)
  {
    push @myissues, $design->get_issues(-sid => $gblock->id);
  }
  print "@theseblocks\n@myissues\n";
}
print {$OUTH} "Make: $make; Nomake: $nomake\n\n";

################################################################################
################################## REPORTING  ##################################
################################################################################
$design->dump_xml($p{OUTPUT});

print {$OUTH} "\n\n";
print {$OUTH} "Wrote $p{OUTPUT}\n\n";
print {$OUTH} $DNA->attitude() . " brought to you by $DNAV\n\n";

print {$OUTH} "\n\n******* $DNAV FINISHED\n\n";
close $OUTH;
exit;


################################################################################
################################# SUBROUTINES  #################################
################################################################################

sub repair_tandem_repeat
{
  my ($GD, $obj, $gene, $problem) = @_;
  my $seq = $obj->seq;
  my ($gstart, $gstop) = ($gene->start, $gene->end);
  my $geneseq = $gene->seq->seq;
  
  #Get in frame
  my $change_start = $problem->start < $gstart ? $gstart  : $problem->start;
  while (($change_start - $gstart) % 3 != 0)
  {
    $change_start--;
  }
  my $change_end = $problem->end > $gstop ? $gstop  : $problem->end;
  while (($change_end - $gstop) % 3 != 0)
  {
    $change_end++;
  }
  my $changeseq = substr $seq, $change_start - 1, $change_end - $change_start + 1;
  my $newseq = $changeseq;

  my $repeat = join q{}, $problem->get_tag_values('repeat');
  #print "\t\tFOUND A TRIMER OF $repeat in $changeseq\n";
  #Test to see if trimer extends behind or ahead an in frame repeat
  my $pep = $CT->{$repeat};
  my $firstcod = substr $seq, $change_start - 1, 3;
  if ($CT->{$firstcod} eq $pep)
  {
    my $prevcod = substr $seq, $change_start - 4, 3;
    my $nextcod = substr $seq, $change_end, 3;
    #print "\t\tEXAMINING $prevcod\t$firstcod\t$nextcod\n";
    while ($CT->{$prevcod} eq $CT->{$firstcod})
    {
      $change_start -= 3;
      $prevcod = substr $seq, $change_start - 4, 3;
    }
    while ($CT->{$nextcod} eq $CT->{$firstcod})
    {
      $change_end += 3;
      $nextcod = substr $seq, $change_end, 3;
    }
    $newseq = substr $seq, $change_start - 1, $change_end - $change_start + 1;
    #print "\t\texpanding from $changeseq to $newseq\n" if ($newseq ne $changeseq);
  }

  #Rotate through possible codons for each position to disrupt the trimer
  $changeseq = $newseq;
  my $changepep = $GD->translate($changeseq);
  my %pephsh = map {$_ => []} split q{}, $changepep;
  foreach my $aa (keys %pephsh)
  {
    my @codons = @{$RCT->{$aa}};
    $pephsh{$aa} = \@codons;
  }
  my $tcv = 0;
  while ($tcv < length $changeseq)
  {
    my $curcod = substr $newseq, $tcv, 3;
    my @codons = @{$pephsh{$CT->{$curcod}}};
    my $newcod = shift @codons;
    substr $newseq, $tcv, 3, $newcod;
    push @codons, $newcod;
    $pephsh{$CT->{$curcod}} = \@codons;
    $tcv += 3;
  }
  
  #Make the repair feature
  my $label = 'repair_' . join q{}, $problem->get_tag_values('label');
  my $type  = join q{}, $problem->get_tag_values('kind');
  my $repair = Bio::SeqFeature::Generic->new(
    -start    => $change_start,
    -end      => $change_end,
    -display_name => $label,
    -primary  => 'sequence_alteration',
    -tag      => {
      oldseq => $changeseq,
      newseq => $newseq,
      label => $label . q{(} . $type . q{_} . $repeat . q{)},
    },
  );
  return $repair;
}

sub repair_direct_repeat
{
  my ($GD, $obj, $gene, $problem) = @_;
  my $seq = $obj->seq;
  my ($gstart, $gstop) = ($gene->start, $gene->end);
  my $geneseq = $gene->seq->seq;
  
  #Get in frame
  my $change_start = $problem->start < $gstart ? $gstart  : $problem->start;
  while (($change_start - $gstart) % 3 != 0)
  {
    $change_start--;
  }
  my $change_end = $problem->end > $gstop ? $gstop  : $problem->end;
  while (($change_end - $gstop) % 3 != 0)
  {
    $change_end++;
  }
  my $changeseq = substr $seq, $change_start - 1, $change_end - $change_start + 1;
  my $newseq = $changeseq;

  my $repeat = join q{}, $problem->get_tag_values('repeat');
  print "\t\tFOUND A DIRECT REPEAT OF $repeat in $changeseq\n";
  #Test to see if trimer extends behind or ahead an in frame repeat


  #Check to see if this overlaps any tandem repeats; if so don't do it
  #my @drlaps = $mask->feature_objects_in_range($change_start, length $changeseq);
  #my @lcrs = grep {$_->primary_tag eq 'low_complexity_region'} @drlaps;
  #my $skipflag = 0;
  #foreach my $lcr (@lcrs)
  #{
  #  my $type = join q{}, $lcr->get_tag_values('note');
  #  my $label = join q{}, $lcr->get_tag_values('label');
  #  if ($type eq 'DIMER' || $type eq 'TRIMER' || $type =~ 'HP')
  #  {
  #    $skipflag++;
  #    $notproblems{$label}++;
  #  }
  #}
  #next if ($skipflag > 0);
  my $replen = length $repeat;
  my $index = index $changeseq, $repeat;
  my $pos = $index;
  next if ($pos < 0);
  while ($pos % 3 != 0)
  {
    $pos--;
  }
  my $end = $index + $replen;
  while ($end % 3 != 0)
  {
    $end++;
  }
  $change_end = $change_start + $end - 1;
  $change_start += $pos;
  #print "\t\t\tgrabbed $pos .. $end\n";
  my $subseq = substr $changeseq, $pos, $end - $pos;
  $newseq = $subseq;
  #my $chcounter = 0;
  #my $flag = 0;
  ##print "\t\t$repeat\t$subseq, $newseq\n";
  while ($newseq eq $subseq)
  {
    $newseq = $GD->codon_juggle(-algorithm => 'balanced', -sequence => $subseq);
  #  $chcounter++;
  #  if ($chcounter == 100 && $newseq eq $subseq)
  #  {
  #    $flag++;
  #    last;
  #  }
  }
  #next if ($flag);
    
  #Make the repair feature
  my $label = 'repair_' . join q{}, $problem->get_tag_values('label');
  my $type  = join q{}, $problem->get_tag_values('kind');
  my $repair = Bio::SeqFeature::Generic->new(
    -start    => $change_start,
    -end      => $change_end,
    -display_name => 'repair_' . $label,
    -primary  => 'sequence_alteration',
    -tag      => {
      oldseq => $subseq,
      newseq => $newseq,
      label => 'repair_' . $label . q{(DR_} . $repeat . q{)},
      type => $type,
    },
  );
  return $repair;
}

sub repair_repeat_load
{
  my ($GD, $p, $obj, $gene, $problem, $gblock) = @_;
  my $seq = $obj->seq;
  my $gseq = $gblock->seq->seq;
  print "working on gblock ", $gblock->start, "\n";
  my ($gstart, $gstop) = ($gene->start, $gene->end);
  my $geneseq = $gene->seq->seq;
  
  #Get in frame
  my $change_start = $problem->start < $gstart ? $gstart  : $problem->start;
  while (($change_start - $gstart) % 3 != 0)
  {
    $change_start--;
  }
  my $change_end = $problem->end > $gstop ? $gstop  : $problem->end;
  while (($change_end - $gstop) % 3 != 0)
  {
    $change_end++;
  }
  my $changeseq = substr $seq, $change_start - 1, $change_end - $change_start + 1;
  my $seqlen = length $changeseq;
  my $gindex = index $gseq, $changeseq;
  die 'GAWD' if ($gindex == -1);
  my $newseq = $GD->codon_juggle(
    -algorithm => 'balanced',
    -sequence => $changeseq
  );
  print "$gseq\n";
  substr $gseq, $gindex, $seqlen, $newseq;
  my $gmask = repeat_mask($GD, $p, $gseq);
  my $times = 0;
  while (scan_for_load($GD, $p, $gmask) || check_for_tandem_repeats($GD, $p, $newseq))
  {
    $newseq = $GD->codon_juggle(
      -algorithm => 'balanced',
      -sequence => $newseq
    );
    substr $gseq, $gindex, $seqlen, $newseq;
    $gmask = repeat_mask($GD, $p, $gseq);
        $times++;
        last if ($times > 1000);
  }

  #print "last $gseq\n";
  #print "mask $gmask\n\n";
  #Make the repair feature
  my $label = 'repair_' . join q{}, $problem->get_tag_values('label');
  my $type  = join q{}, $problem->get_tag_values('kind');
  my $repair = Bio::SeqFeature::Generic->new(
    -start    => $change_start,
    -end      => $change_end,
    -display_name => $label,
    -primary  => 'sequence_alteration',
    -tag      => {
      oldseq => $changeseq,
      newseq => $newseq,
      label =>  $label,
      type => $type,
    },
  );
  return $repair;  
}

sub check_for_tandem_repeats
{
  my ($GD, $p, $seq) = @_;
  my @troubles;
  
  #homopolymers
  my @BASES = qw(A T C G);
  foreach my $base (@BASES)
  {
    my $coords = $GD->find_runs(
      -sequence => $seq,
      -pattern => $base,
      -minimum => $p->{'HP_' . $base . '_LENGTH'},
    );
    push @troubles, @{$coords};
  }
  
  my @dimers;
  my @trimers;
  my %complexes;
  foreach my $a (@BASES)
  {
    foreach my $b (@BASES)
    {
      if ($a ne $b)
      {
        my $dimer = $a . $b;
        my $dcoords = $GD->find_runs(
          -sequence => $seq,
          -pattern  => $dimer,
          -minimum  => $p->{DIMER_LIMIT},
        );
        push @troubles, @{$dcoords};
      }
      foreach my $c (@BASES)
      {
        if (! ($a eq $b && $b eq $c))
        {
          my $trimer = $a . $b . $c;
          my $trcoords = $GD->find_runs(
            -sequence => $seq,
            -pattern  => $trimer,
            -minimum  => $p->{TRIMER_LIMIT},
          );
          push @troubles, @{$trcoords};
        }
      }
    }
  }
  return scalar @troubles;
}

sub repeat_mask
{
  my ($GD, $p, $seq) = @_;
  my $seqlen = length $seq;  
  my %hsh = ();
  my $x = 0;
  my %pins;
  while ($x < $seqlen - $p->{DR_LENGTH})
  {
    my $stem = substr $seq, $x, $p->{DR_LENGTH};
    $pins{$stem}++;
    $x++;
  }
  my $mask = '0' x $seqlen;
  foreach my $stem (sort {$pins{$a} <=> $pins{$b}} keys %pins)
  {
    my $stemlen = length $stem;
    my $rep = '1' x $stemlen;
    my $positions = $GD->positions(
      -sequence => $seq,
      -query => $stem,
      -reverse_complement => 1
    );
    my @poses = sort {$a <=> $b} keys %{$positions};
    next if (scalar @poses < 2);
    #print "\t found $stem; @poses";
    foreach my $pos (@poses)
    {
      substr $mask, $pos, $stemlen, $rep;
    }
  }
  return $mask;
}

sub scan_for_load
{
  my ($GD, $p, $mask) = @_;
  my $len = length $mask;
  my $start = 0;
  my $flag = 0;
  while ($start < $len - $p->{DR_LOCAL_WINDOW})
  {
    my $rawscore = substr $mask, $start, $p->{DR_LOCAL_WINDOW};
    my @score = grep {$_ == 1} split q{}, $rawscore;
    my $ratio = 100 * (scalar @score) / $p->{DR_LOCAL_WINDOW};
    my $repratio = sprintf "%.1f", $ratio;
    $flag++ if ($repratio >= $p->{DR_LOCAL_RATIO});
    $start++;
  }
  return $flag;

}

__END__


=head1 NAME

  DNAssemble_155_CorrectGBlocks.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

  Optional arguments:

    -o,  --output : a filepath for dumping output

    -h,  --help : display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, GeneDesign developers
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