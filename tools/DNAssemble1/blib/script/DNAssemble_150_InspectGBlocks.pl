#!/usr/bin/env perl

use Bio::DNAssemble;
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
my $DNAV = 'DNAssemble_150_InspectGBlocks_' . $VERSION;
my $DNAS = '_150';

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

##ENDDEBUG

################################################################################
################################## INSPECTING ##################################
################################################################################
my @gblocks = $design->get_constructs(-kind => 'building_block');
my @noissues;
my @someissues;
my $badgblockcount = 0;
my $goodgblockcount = 0;
my $goodconstructcount = 0;
my $badconstructcount = 0;
foreach my $gblock (@gblocks)
{
  print "Examining ", $gblock->id, "\n";
  my @issues = ();
  my @troubles = ();
  ################################ GC PERCENTAGE
  ## Global
  push @troubles, markup_global_gc($GD, \%p, $gblock);
  ## Regional
  push @troubles, markup_regional_gc($GD, \%p, $gblock);
  ## Local
  push @troubles, markup_local_gc($GD, \%p, $gblock);
  ## Terminal
  push @troubles, markup_terminal_gc($GD, \%p, $gblock);

  ################################ REPEATS
  ## Tandem Repeats
  push @troubles, markup_tandem_repeats($GD, \%p, $gblock);
  ## Direct Repeats
  push @troubles, markup_direct_repeats($GD, \%p, $gblock);
  ## Inverted Repeats
  push @troubles, markup_inverted_repeats($GD, \%p, $gblock);
  
  
  if (scalar @troubles)
  {
    $badgblockcount++;
    print "\t gblock ", $gblock->start, " has ", scalar @troubles, " issues\n";
    foreach my $issue (@troubles)
    {
      push @issues, $issue;
      my $attr = {
        #kind => $issue->primary_tag,
        start => $issue->start,
        end => $issue->end,
      };
      if ($issue->has_tag('repeat'))
      {
        $attr->{repeat} = join q{}, $issue->get_tag_values('repeat');
      }
      if ($issue->has_tag('kind'))
      {
        $attr->{kind} = join q{}, $issue->get_tag_values('kind');
      }
      $design->critique(-construct => $gblock, -attributes => $attr);
    }
  }
  else
  {
    $goodgblockcount++;
  }
  if (scalar @issues)
  {
    push @someissues, $gblock;
    $badconstructcount++;
  }
  else
  {
    push @noissues, $gblock;
    $goodconstructcount++;
  }
}
print {$OUTH} "$goodgblockcount good gblocks; $badgblockcount bad gblocks\n\n";
print {$OUTH} "$goodconstructcount good constructs; $badconstructcount bad constructs\n\n";

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

sub markup_global_gc
{
  my ($GD, $p, $gblock) = @_;
  my $seq = $gblock->sequence;
  my $seqlen = length $seq;
  my $globalstats = $GD->count($seq);
  my $globalGC = $globalstats->{GCp};
  my @troubles;
  if ($globalGC > $p->{GC_GLOBAL_MAX})
  {
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start,
      -end      => $gblock->end,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        kind => 'gc_global',
        note => $globalGC,
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;
  }
  elsif ($globalGC < $p->{GC_GLOBAL_MIN})
  {
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => 1,
      -end      => $seqlen,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        kind => 'gc_global',
        note => $globalGC,
        label => 'error_' . $p->{ERROR_NUMBER},
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;
  }
  return @troubles;
}

sub markup_regional_gc
{
  my ($GD, $p, $gblock) = @_;
  my @troubles;
  my $region_problems = $GD->GC_windows(
    -sequence => $gblock->sequence,
    -window   => $p->{GC_REGIONAL_WINDOW},
    -minimum  => $p->{GC_REGIONAL_MIN},
    -maximum  => $p->{GC_REGIONAL_MAX},
  );
  my ($regionlo, $regionhi) = @{$region_problems};
  my @lotroubles;
  my @hitroubles;
  if (! $p->{SKIP_WARNINGS})
  {
    foreach my $coord (@{$regionlo})
    {
      my ($lef, $rig) = @{$coord};
      my $trouble = Bio::SeqFeature::Generic->new(
        -start    => $gblock->start + $lef - 1,
        -end      => $gblock->start + $rig - 1,
        -display_name => 'error_' . $p->{ERROR_NUMBER},
        -primary  => 'low_complexity_region',
        -tag      => {
          kind => 'gc_regional',
          label => 'error_' . $p->{ERROR_NUMBER},
          gblock => $gblock->start,
        },
      );
      $p->{ERROR_NUMBER}++;
      push @lotroubles, $trouble;
    }
  }
  foreach my $coord (@{$regionhi})
  {
    my ($lef, $rig) = @{$coord};
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start + $lef - 1,
      -end      => $gblock->start + $rig - 1,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        kind => 'gc_regional',
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @hitroubles, $trouble;
  }
  my @lofeatures = @{coalesce(\@lotroubles)};
  my @hifeatures = @{coalesce(\@hitroubles)};
  foreach my $feature (@lofeatures, @hifeatures)
  {
    push @troubles, $feature;
  }
  return @troubles;
}

sub markup_local_gc
{
  my ($GD, $p, $gblock) = @_;
  my @troubles;
  my $region_problems = $GD->GC_windows(
    -sequence => $gblock->sequence,
    -window   => $p->{GC_LOCAL_WINDOW},
    -minimum  => $p->{GC_LOCAL_MIN},
    -maximum  => $p->{GC_LOCAL_MAX},
  );
  my ($regionlo, $regionhi) = @{$region_problems};
  my @lotroubles;
  my @hitroubles;
  if (! $p->{SKIP_WARNINGS})
  {
    foreach my $coord (@{$regionlo})
    {
      my ($lef, $rig) = @{$coord};
      my $trouble = Bio::SeqFeature::Generic->new(
        -start    => $gblock->start + $lef - 1,
        -end      => $gblock->start + $rig - 1,
        -display_name => 'error_' . $p->{ERROR_NUMBER},
        -primary  => 'low_complexity_region',
        -tag      => {
          kind => 'gc_local',
          label => 'error_' . $p->{ERROR_NUMBER},
          gblock => $gblock->start,
        },
      );
      $p->{ERROR_NUMBER}++;
      push @lotroubles, $trouble;
    }
  }
  foreach my $coord (@{$regionhi})
  {
    my ($lef, $rig) = @{$coord};
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start + $lef - 1,
      -end      => $gblock->start + $rig - 1,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        kind => 'gc_local',
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @hitroubles, $trouble;
  }
  my @lofeatures = @{coalesce(\@lotroubles)};
  my @hifeatures = @{coalesce(\@hitroubles)};
  foreach my $feature (@lofeatures, @hifeatures)
  {
    push @troubles, $feature;
  }
  return @troubles;
}

sub markup_terminal_gc
{
  my ($GD, $p, $gblock) = @_;
  my @troubles;
  my $seq = $gblock->sequence;
  my $seqlen = length $seq;
  my $termfivestats = $GD->count(substr $seq, 0, $p->{GC_TERMINAL_WINDOW});
  my $termthreestats = $GD->count(substr $seq, -$p->{GC_TERMINAL_WINDOW});
  my $termfiveGC = $termfivestats->{GCp};
  if ($termfiveGC > $p->{GC_TERMINAL_MAX})
  {
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start + 1 - 1,
      -end      => $gblock->start + $p->{GC_TERMINAL_WINDOW} - 1,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'T_GC>' . $p->{GC_TERMINAL_MAX},
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;
  }
  elsif ($termfiveGC < $p->{GC_TERMINAL_MIN})
  {
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start + 1 - 1,
      -end      => $gblock->start + $p->{GC_TERMINAL_WINDOW} - 1,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'T_GC<' . $p->{GC_TERMINAL_MIN},
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;
  }
  my $termthreeGC = $termthreestats->{GCp};
  if ($termthreeGC > $p->{GC_TERMINAL_MAX})
  {
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->end - $p->{GC_TERMINAL_WINDOW},
      -end      => $gblock->end,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'T_GC>' . $p->{GC_TERMINAL_MAX},
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;
  }
  elsif ($termthreeGC < $p->{GC_TERMINAL_MIN})
  {
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->end - $p->{GC_TERMINAL_WINDOW},
      -end      => $gblock->end,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'T_GC<' . $p->{GC_TERMINAL_MIN},
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start,
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;
  }
  return @troubles;
}

sub markup_tandem_repeats
{
  my ($GD, $p, $gblock) = @_;
  my $seq = $gblock->sequence;
  my @troubles;
  my @reptroubles;
  
  #homopolymers
  my @BASES = qw(A T C G);
  foreach my $base (@BASES)
  {
    my $coords = $GD->find_runs(
      -sequence => $seq,
      -pattern => $base,
      -minimum => $p->{'HP_' . $base . '_LENGTH'},
    );
    foreach my $coord (@{$coords})
    {
      my ($lef, $rig) = @{$coord};
      my $trouble = Bio::SeqFeature::Generic->new(
        -start    => $gblock->start + $lef - 1,
        -end      => $gblock->start + $rig - 1,
        -display_name => 'error_' . $p->{ERROR_NUMBER},
        -primary  => 'low_complexity_region',
        -tag      => {
          repeat => $base,
          kind => 'homopolymer',
          label => 'error_' . $p->{ERROR_NUMBER},
          gblock => $gblock->start,
        },
      );
      $p->{ERROR_NUMBER}++;
      push @troubles, $trouble;
    }
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
        foreach my $dcoord (@{$dcoords})
        {
          my ($lef, $rig) = @{$dcoord};
          my $trouble = Bio::SeqFeature::Generic->new(
            -start    => $gblock->start + $lef - 1,
            -end      => $gblock->start + $rig - 1,
            -primary  => 'low_complexity_region',
            -display_name => 'error_' . $p->{ERROR_NUMBER},
            -tag      => {
              repeat => $dimer,
              kind => 'dimer',
              label => 'error_' . $p->{ERROR_NUMBER},
              gblock => $gblock->start, 
            },
          );
          $p->{ERROR_NUMBER}++;
          push @reptroubles, $trouble;
        }
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
          foreach my $trcoord (@{$trcoords})
          {
            my ($lef, $rig) = @{$trcoord};
            my $trouble = Bio::SeqFeature::Generic->new(
              -start    => $gblock->start + $lef - 1,
              -end      => $gblock->start + $rig - 1,
              -display_name => 'error_' . $p->{ERROR_NUMBER},
              -primary  => 'low_complexity_region',
              -tag      => {
                kind => 'trimer',
                repeat => $trimer,
                number => ($rig - $lef + 1) / 3,
                label => 'error_' . $p->{ERROR_NUMBER},
                gblock => $gblock->start, 
              },
            );
            $p->{ERROR_NUMBER}++;
            push @reptroubles, $trouble;
          }
        }
      }
    }
  }

  my @features = @{coalesce(\@reptroubles)};
  foreach my $feature (@features)
  {
    push @troubles, $feature;
  }
  return @troubles;
}

sub markup_direct_repeats
{
  my ($GD, $p, $gblock) = @_;
  my @troubles;
  my $seq = $gblock->sequence;
  my $seqlen = length $seq;  
  my %hsh = ();
  my $x = 0;
  my %pins;
  while ($x <= $seqlen - ($p->{DR_LENGTH}))
  {
    my $stem = substr $seq, $x, $p->{DR_LENGTH};
    $pins{$stem}++;
    $x++;
  }
  my $mask = '0' x $seqlen;
  my %seens = ();
  my $rep = '1' x $p->{DR_LENGTH};
  foreach my $stem (sort keys %pins)
  {
    my $positions = $GD->positions(
      -sequence => $seq,
      -query => $stem,
      -reverse_complement => 1 
    );
    my @poses = sort {$a <=> $b} keys %{$positions};
    next if (scalar @poses < 2);
    next if (exists $seens{$GD->rcomplement($stem)});
    $seens{$stem}++;
    foreach my $index (@poses)
    {
      substr $mask, $index, $p->{DR_LENGTH}, $rep;
    }
  }
  my $pos = -1;
  while ($pos <= $seqlen - $p->{DR_LOCAL_WINDOW})
  {
    $pos++;
    my $lrawscore = substr $mask, $pos, $p->{DR_LOCAL_WINDOW};
    my @lscore = grep {$_ == 1} split q{}, $lrawscore;
    my $lratio = 100 * (scalar @lscore) / $p->{DR_LOCAL_WINDOW};
    my $lrepratio = sprintf "%.1f", $lratio;
    if ($lrepratio >= $p->{DR_LOCAL_RATIO}) 
    {
      my $start = $pos + 1;
      my $end = $pos + $p->{DR_LOCAL_WINDOW};
      #print "$lrepratio% (", scalar @lscore, "); $start..$end ";
      #print join q{,}, split q{}, $lrawscore;
      #print "\n";
      my $trouble = Bio::SeqFeature::Generic->new(
        -start    => $gblock->start + $start - 1,
        -end      => $gblock->start + $end - 1,
        -display_name => 'error_' . $p->{ERROR_NUMBER},
        -primary  => 'low_complexity_region',
        -tag      => {
          note => 'local',
          kind => 'repeat_load',
          label => 'error_' . $p->{ERROR_NUMBER},
          gblock => $gblock->start, 
        },
      );
      $p->{ERROR_NUMBER}++;
      push @troubles, $trouble;        
    }
  } 
  #Total repeats
  my @score = grep {$_ == 1} split q{}, $mask;
  my $ratio = 100 * (scalar @score) / $seqlen;
  my $repratio = sprintf "%.1f", $ratio;
  if ($repratio >= $p->{DR_GLOBAL_RATIO}) 
  {
    print "$repratio!\n";
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start,
      -end      => $gblock->end,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'global',
        kind => 'repeat_load',
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start, 
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;        
  }

  #Terminal repeats
  my $failbit = '1' x $p->{DR_LENGTH};
  my $fterbit = substr $mask, 0, $p->{DR_LENGTH};
  if ($fterbit eq $failbit)
  {
    print "gotta5terminal repeat: $fterbit\n";
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->start,
      -end      => $gblock->start + $p->{DR_LENGTH} - 1,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'local',
        kind => 'terminal_repeat',
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start, 
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble; 
  }
  my $tterbit = substr $mask, -$p->{DR_LENGTH};
  if ($tterbit eq $failbit)
  {
    print "gotta3terminal repeat: $tterbit\n";
    my $trouble = Bio::SeqFeature::Generic->new(
      -start    => $gblock->end - $p->{DR_LENGTH} + 1,
      -end      => $gblock->end,
      -display_name => 'error_' . $p->{ERROR_NUMBER},
      -primary  => 'low_complexity_region',
      -tag      => {
        note => 'local',
        kind => 'terminal_repeat',
        label => 'error_' . $p->{ERROR_NUMBER},
        gblock => $gblock->start, 
      },
    );
    $p->{ERROR_NUMBER}++;
    push @troubles, $trouble;      
  }
  #my $printmask = join q{,}, split q{}, substr $mask, 0;
  #print "\"$printmask\"\n";
  #print "$seq\n";
  #print "$mask\n";
  #print "\n\n";
  @troubles = @{coalesce(\@troubles)};
  return @troubles;
}

sub markup_inverted_repeats
{
  my ($GD, $p, $gblock) = @_;
  my $seq = $gblock->sequence;
  my $seqlen = length $seq;
  my @troubles;
  my %hsh = ();
  my %pins;
  my $x = 0;
  my $onnastem = 0;
  my $length = $p->{IR_LENGTH_MIN};
  my $limit = $x + $length;
  my $stem = substr $seq, $x, $length;
  my ($start, $end) = (undef, undef);
  while ($x < $seqlen - $p->{IR_LENGTH_MIN})
  {
    my $mets = $GD->rcomplement($stem);
    #replace this with mismatch detection if desired
    my $position = index $seq, $mets, $limit;
    if ($position == -1 && $onnastem == 0)
    {
      $x++;
    }
    elsif ($position == -1 && $onnastem == 1)
    {
      $stem = substr $stem, 0, $length - 1;
      $pins{$stem} = [$start, $end];
      $length = $p->{IR_LENGTH_MIN};
      $onnastem = 0;
      $x = $x + length $stem;
    }
    else
    {
      $onnastem = 1;
      $start = $x + 1;
      $end = $position + $length;
      $length++;
    }
    $limit = $x + $length;
    $stem = substr $seq, $x, $length;
  }
  foreach my $stem (keys %pins)
  {
    my ($start, $stop) = @{$pins{$stem}};
    if ($stop - $start + 1 < $p->{IR_WIDTH} || length $stem > $p->{IR_LENGTH_MAX})
    {
      my $trouble = Bio::SeqFeature::Generic->new(
        -start    => $gblock->start + $start - 1,
        -end      => $gblock->start + $stop - 1,
        -display_name => 'error_' . $p->{ERROR_NUMBER},
        -primary  => 'hairpin',
        -tag      => {
          kind => 'hairpin',
          repeat => $stem,
          label => 'error_' . $p->{ERROR_NUMBER},
          gblock => $gblock->start, 
        },
      );
      $p->{ERROR_NUMBER}++;
      push @troubles, $trouble
    }
  }
  return @troubles;
}

sub coalesce
{
  my ($arr) = @_;
  my @thisarr = sort {$a->start <=> $b->start} @{$arr};
  my @features = ();
  while (scalar @thisarr)
  {
    my $feature = shift @thisarr;
    my $thisend = $feature->end;
    my $thisstart = $feature->start;
    while (scalar @thisarr && $thisarr[0]->start <= $feature->end)
    {
      my $nextend = $thisarr[0]->end;
      my $nextstart = $thisarr[0]->start;
      my $label = join q{}, $feature->get_tag_values('kind');
      if ($label =~ 'repeat_load' && $feature->has_tag('repeat'))
      {
        my $ediff = $nextend - $thisend;
        my $sdiff = $nextstart - $thisstart;
        if ($nextstart > $thisstart && $nextstart < $thisend)
        {
          my $repeat = join q{}, $feature->get_tag_values('repeat');
          my $extend = join q{}, $thisarr[0]->get_tag_values('repeat');
          $feature->remove_tag('repeat');
          $repeat .= substr $extend, -1;
          $feature->add_tag_value('repeat', $repeat);
          $feature->end($nextend);
          shift @thisarr;
        }
        else
        {
          last;
        }
      }
      elsif ($label =~ 'trimer')
      {
        shift @thisarr;
      }
      else
      {
        $feature->end($nextend);
        shift @thisarr;
      }
    }
    push @features, $feature;
  }
  return \@features;
}



__END__


=head1 NAME

  DNAssemble_150_InspectGBlocks.pl

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