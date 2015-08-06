#!/usr/bin/env perl

use Bio::GeneDesign;
use Bio::DNAssemble;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use Readonly;
use POSIX;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_122_Carve_Simple_' . $VERSION;
my $DNAS = '_122';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
  'input=s'              => \$p{INPUT},
  'output:s'             => \$p{OUTPUT},
  'gblock_length=i'      => \$p{GBLOCK_LENGTH},
  'gblock_length_max=i'  => \$p{GBLOCK_LENGTH_MAX},
  'gblock_length_min=i'  => \$p{GBLOCK_LENGTH_MIN},
  'gblock_overlap=i'     => \$p{GBLOCK_OVERLAP},
  'gblock_overlap_min=i' => \$p{GBLOCK_OVERLAP_MIN},
  'gc_terminal_min=i'    => \$p{GC_TERMINAL_MIN},
  'gc_terminal_max=i'    => \$p{GC_TERMINAL_MAX},
  'gc_terminal_window=i' => \$p{GC_TERMINAL_WINDOW},
  'assemblyTm=i'         => \$p{ASSEMBLYTM},
  'dr_length=i'          => \$p{DR_LENGTH},
  'logfile:s'            => \$p{LOGPATH},
  'help'                 => \$p{HELP}
);


################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>'NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE')
  if ($p{HELP});

my $DNA = Bio::DNAssemble->new();
my $GD = $DNA->GD();

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
my $filename = $design->filename();
$filename = $1 if ($filename =~ m{(.+)\_}msix);
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.xml};

Readonly my $SEARCHLIMIT  =>  500;

################################################################################
################################# CONFIGURING ##################################
################################################################################


################################################################################
################################### CARVING ####################################
################################################################################
my @chunks = $design->get_constructs(-kind => 'unknown');
foreach my $chunk ( @chunks )
{
  # Gather data about the chunk
  #
  #
  my $chunkname = $chunk->id();
  my $chseq = uc $chunk->sequence();
  my $chlen = length $chseq;
  print {$OUTH} "Gblocking $chunkname ($chlen bp)...\n";
  my $gblocknum = 1;
  my $sublet = 'A';
  my $params = {
      stitchtemp => $p{ASSEMBLYTM}
  };
  
  my @gblockables = $chunk->get_features(-primary_tag => 'gblockable');  
  foreach my $gblockable (@gblockables)
  {
    my ($gstart, $gend) = ($gblockable->start, $gblockable->end);
    my $offset = $gstart - 1;
    my $gblseq = substr $chseq, $gstart - 1, $gend - $gstart + 1;
    my $gbllen = length $gblseq;
    my $ps = {
      chunkseq   => $chseq,
      maxbblen   => $p{GBLOCK_LENGTH_MAX},
      minbblen   => $p{GBLOCK_LENGTH_MIN},
      tarbblen   => $p{GBLOCK_LENGTH},
      overlap    => $p{GBLOCK_OVERLAP},
      minbblap   => $p{GBLOCK_OVERLAP_MIN},
      maxgcterm  => $p{GC_TERMINAL_MAX},
      mingcterm  => $p{GC_TERMINAL_MIN},
      wingcterm  => $p{GC_TERMINAL_WINDOW},
    };

    my $return = carve_building_blocks($ps);
    my $BBS = $return->{constructs};
    my $num = scalar @{$BBS};

    # If it cannot be carved, abandon ship
    #
    #
    if ($num < 1)
    {
      print {$OUTH} "\t$chunkname cannot be made into building blocks... skipping\n";
      $design->weep(-construct => $chunk, -note=> 'unfragmentable');
      $chunk->deliverable('true');
      $chunk->method('magic');
      last;
    }

    my $mask = find_repeats($GD, \%p, $gblseq);
    my $gblockcount = 0;
    my $glimit = scalar @{$BBS};
    my $failbit = '1' x $p{DR_LENGTH};

    my $fragname = $chunkname . q{.} . $sublet;
    my $patts = {number => 1, method => 'amplification'};
    
    foreach my $bb (@{$BBS})
    {
      $gblockcount++;
      my ($start, $end) = ($bb->{start}, $bb->{end});
      my $gblockseq = substr $gblseq, $start - 1, $end - $start + 1;
      my $gmask = find_repeats($GD, \%p, $gblockseq);
      if ($gblockcount != 1)
      {
        my $termfivestats = $GD->count(substr $chseq, $start - 1, $p{GC_TERMINAL_WINDOW});
        my $termfiveGC = $termfivestats->{GCp};
        my $gmaskbit = substr $gmask, 0, $p{DR_LENGTH};
        #print "\tbeginning adjustment of start ($termfiveGC) from $start to ";
        while ($gmaskbit eq $failbit || $termfiveGC > $p{GC_TERMINAL_MAX} || $termfiveGC < $p{GC_TERMINAL_MIN})
        {
         #print "(th $termfiveGC)" if ($termfiveGC > $p{GC_TERMINAL_MAX});
         #print "(tl $termfiveGC)" if ($termfiveGC < $p{GC_TERMINAL_MIN});
         #print "(mb $gmaskbit)" if ($gmaskbit eq $failbit);
          $start--;
          $gblockseq = substr $chseq, $start - 1, $end - $start + 1;
          $gmask = find_repeats($GD, \%p, $gblockseq);
          $gmaskbit = substr $gmask, 0, $p{DR_LENGTH};
          $termfivestats = $GD->count(substr $chseq, $start - 1, $p{GC_TERMINAL_WINDOW});
          $termfiveGC = $termfivestats->{GCp};
        }
        #print " $start\n";
      }
      my $termthreestats = $GD->count(substr $chseq, $end - $p{GC_TERMINAL_WINDOW}, $p{GC_TERMINAL_WINDOW});
      my $termthreeGC = $termthreestats->{GCp};  
      my $gmaskbit = substr $gmask, -$p{DR_LENGTH};
      my $oldend = $end;
      #print "beginning adjustment of end ($termthreeGC) from $end to ";
      while ($gmaskbit eq $failbit || $termthreeGC > $p{GC_TERMINAL_MAX} || $termthreeGC < $p{GC_TERMINAL_MIN})
      {
       #print "(th $termthreeGC)" if ($termthreeGC > $p{GC_TERMINAL_MAX});
       #print "(tl $termthreeGC)" if ($termthreeGC < $p{GC_TERMINAL_MIN});
       #print "(mb $gmaskbit)" if ($gmaskbit eq $failbit);
        $end++;
        $gblockseq = substr $chseq, $start - 1, $end - $start + 1;
        $gmask = find_repeats($GD, \%p, $gblockseq);
        $gmaskbit = substr $gmask, -$p{DR_LENGTH};
        $termthreestats = $GD->count(substr $chseq, $end - $p{GC_TERMINAL_WINDOW}, $p{GC_TERMINAL_WINDOW});
        $termthreeGC = $termthreestats->{GCp};
        if ($end == ($oldend + 20))
        {
          print "BAILING\n";
          $end = $oldend;
          last;
        }
      }
      #print " $end\n";                              
      #$gblockseq = substr $chseq, $start - 1, $end - $start + 1;
      #$gmask = find_repeats($GD, \%p, $gblockseq);
      #print "gblock: $gblockseq\n";
      #print "\ngmask: $gmask\n\n";
      #print "\n\n\n";
            
      my $bbname = $fragname . q{.} . $DNA->pad($patts->{number}, 2);
      my %mparams = %{$params};
      $mparams{bbhsh} = $bb;
      $mparams{bbseq} = bb_seq($bb, $chseq);
      $mparams{name} = $bbname;
      $mparams{offset} = $offset;
      my ($bbcon, $slcon, $srcon) = make_bbcon(\%mparams);
      my $pool = $chunk->add_pool($patts);
      $pool->add_subconstructs([$slcon, $srcon, $bbcon]);
      $patts->{number}++;
    }
    $sublet++;
  }
}

################################################################################
################################## REPORTING ###################################
################################################################################
$design->dump_xml($p{OUTPUT});

print {$OUTH} "\n\n";
print {$OUTH} "Wrote $p{OUTPUT}\n\n";
print {$OUTH} $DNA->attitude() . " brought to you by $DNAV\n\n";

print {$OUTH} "\n\n******* $DNAV FINISHED\n\n";
close $OUTH;
exit;


################################################################################
################################# SUBROUTINES ##################################
################################################################################
sub carve_building_blocks
{
  my ($ps) = @_;
  my $v = 1;
  my $SGD = Bio::GeneDesign->new();
  my $return = {error => undef, warning => undef, constructs => []};

  my $chseq = $ps->{chunkseq};
  my $chlen = length $chseq;
  my $chseqobj  = Bio::Seq->new(-id => 'chunk', -seq => $chseq);

  my $tarbblap = $ps->{overlap};
  my $minbblap = $ps->{minbblap};
  my $tarbblen = $ps->{tarbblen};
  $tarbblen = $chlen < $tarbblen  ? $chlen  : $tarbblen;
  my $maxbblen = $ps->{maxbblen};
  my $minbblen = $ps->{minbblen};
  my $maxgcterm = $ps->{maxgcterm};
  my $mingcterm = $ps->{mingcterm};
  my $wingcterm = $ps->{wingcterm};

  ## Decide what size bbs to go for and how many
  ## Adjust the target building block size so as to avoid outliers.
  my $tarnum = $chlen > $maxbblen
              ? ceil(($chlen) / ($tarbblen - $tarbblap))
              : 1;
  my $diff = $chlen - (($tarnum * $tarbblen) - ($tarbblap * ($tarnum - 1)));
  print {$OUTH} "\ttarget: $tarnum bbs size $tarbblen bp rem $diff\n" if ($v);
  if (abs($diff) >= $tarnum)
  {
    my $rem = sprintf '%.0f', $diff / $tarnum;
    $tarbblen = $tarbblen + $rem;
    $diff = $diff - ($tarnum * $rem);
  }
  print {$OUTH} "\t final: $tarnum bbs size $tarbblen bp rem $diff\n" if ($v);
  my @BBS;

  if ($tarnum > 1)
  {
    # Find all of the overlaps in the construct that fit the length profile
    # the overlap sequence as a Bio::Seq object
    #
    my $laps = [];
    my $lef = 0;
    my $id = 0;
    my $rig = $tarbblap;
    while ($rig <= $chlen)
    {
      my $ol = substr $chseq, $lef, $rig - $lef + 1;
      my $fstats = $SGD->count(substr $ol, 0, $wingcterm);
      my $tstats = $SGD->count(substr $ol, -$wingcterm);
      my $fGC = $fstats->{GCp};
      my $tGC = $tstats->{GCp};
      if ($fGC < $mingcterm || $fGC > $maxgcterm
       || $tGC < $mingcterm || $tGC > $maxgcterm)
      {
        $lef++;
        $rig = $tarbblap  ? $lef + $tarbblap  : $lef + 7;
        next;
      }
      my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
      my $lefta   = Bio::Annotation::SimpleValue->new(-value => $lef + 1);
      $lapobj->add_Annotation('left', $lefta);
      my $righta  = Bio::Annotation::SimpleValue->new(-value => $rig + 1);
      $lapobj->add_Annotation('right', $righta);
      push @{$laps}, $lapobj;
      $lef++;
      $id++;
      $rig = $tarbblap  ? $lef + $tarbblap  : $lef + 7;
    }
    my $subreturn = $DNA->choose_overlaps(
      -construct_sequence => $chseqobj,
      -overlap_list       => $laps,
      -overlap_len_min    => $minbblap,
      -span_len           => $tarbblen,
      -span_len_min       => $minbblen,
      -span_len_max       => $maxbblen,
    );
    my @chosenlaps = @{$subreturn->{overlaps}};
    if (! scalar @chosenlaps)
    {
      $return->{error} = $subreturn->{error};
      print "NO OVERLAPS RETURNED! ", $subreturn->{error}, "\n";
      return $return;
    }

    #Make the building blocks
    my $lpos = 1;
    while (scalar @chosenlaps)
    {
      my $lap = shift @chosenlaps;
      my $atts = {start => $lpos, end => $lap->[3]};
      push @BBS, $atts;
      $lpos = $lap->[2];
    }
    my $atts = {start => $lpos, end => $chlen};
    push @BBS, $atts;
  }
  else
  {
    my $bbseq = $chseq;
    my $atts = {start => 1, end => $chlen};
    push @BBS, $atts;
  }

  print {$OUTH} "\n\n" if ($v);
  $return->{constructs} = \@BBS;
  return $return;
}

sub make_bbcon
{
  my ($ps) = @_;
  my $stemp = $ps->{stitchtemp};
  my ($bbhsh, $bbseq, $name) = ($ps->{bbhsh}, $ps->{bbseq}, $ps->{name});
  my $offset = $ps->{offset};

  # Begin processing sequence; make insert stitches
  #
  my ($lprimer, $rprimer) = $GD->make_amplification_primers(
    -sequence    => $bbseq,
    -temperature => $stemp,
  );

  # Create and annotate the constructs
  #
  my $buildingb = $design->create_construct({
    id        => $name,
    start     => $bbhsh->{start} + $offset,
    end       => $bbhsh->{end} + $offset,
    kind      => 'building_block',
    method    => 'gblock',
    sequence  => $bbseq,
  });

  my $soligol = $design->create_construct({
    id       => $name . q{.F},
    kind     => 'stitching_oligo',
    method   => 'order',
    orient   => 'F',
    sequence => $lprimer,
  });

  my $soligor = $design->create_construct({
    id       => $name . q{.R},
    kind     => 'stitching_oligo',
    method   => 'order',
    orient   => 'R',
    sequence => $rprimer,
  });

  return ($buildingb, $soligol, $soligor);
}

sub bb_seq
{
  my ($bb, $chseq) = @_;
  my $start = $bb->{start};
  my $seq = substr $chseq, $start - 1, $bb->{end} - $start + 1;
  return $seq;
}

sub find_repeats
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
    #no need to worry if they're more than a gblock's width away from each other
    for my $index (0 .. (scalar @poses) - 2)
    {
      my $pos = $poses[$index];
      my $npos = $poses[$index+1];
      my $range = $npos + $stemlen - $pos;
      if ($range <= $p->{GBLOCK_LENGTH_MAX})
      {
        substr $mask, $pos, $stemlen, $rep;
        substr $mask, $npos, $stemlen, $rep;
      }
    }
  }
  return $mask;
}

__END__

=head1 NAME

  DNAssemble_121_Carve_Simple.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

    -c,  --cloningvector: A vector to serve as the cloning plasmid.

  Optional arguments:

    -o,  --output : a filepath for dumping output

    -g, --gblock_length : the length in bp of building blocks
            Default is 1000 bp

    -ma, --maxlength : the maximum length a building block is allowed to be.
            Default is 1023 bp

    -mi, --minlength : the minimum length a building block is allowed to be.
            Default is 200 bp

    -la, --lap : the target overlap between building blocks.
            Defaults is 40 bp

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