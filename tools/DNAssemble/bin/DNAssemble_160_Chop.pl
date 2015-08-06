#!/usr/bin/env perl

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
my $DNAV = 'DNAssemble_160_Chop_' . $VERSION;
my $DNAS = '_160';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'       => \$p{INPUT},
      'output:s'      => \$p{OUTPUT},
      'olilen:i'      => \$p{OLILEN},
      'laptemp:i'     => \$p{LAPTEMP},
      'poolsize:i'    => \$p{POOLSIZE},
      'maxpoolnum:i'  => \$p{MAXPOOLNUM},
      'maxolilen:i'   => \$p{MAXOLILEN},
      'minolilen:i'   => \$p{MINOLILEN},
      'minlaplen:i'   => \$p{MINLAPLEN},
      'tolerance:f'   => \$p{TMTOL},
      'verbose:s'     => \$p{VERBOSE},
      'logfile:s'     => \$p{LOGPATH},
      'help'          => \$p{HELP}
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

Readonly my $TAROLLEN     =>  150;
Readonly my $MINOLLEN     =>   45;
Readonly my $MAXOLLEN     =>  200;
Readonly my $POOLSIZE     =>    5;
Readonly my $MAXPOOLNUM   =>    2;
Readonly my $STITCH       =>   70;
Readonly my $MINOLLAP     =>   10;
Readonly my $TMTOL        =>    0.5;
Readonly my $SEARCHLIMIT  =>  500;

$p{OLILEN}     = $p{OLILEN}     || $TAROLLEN;
$p{MINOLILEN}  = $p{MINOLILEN}  || $MINOLLEN;
$p{MAXOLILEN}  = $p{MAXOLILEN}  || $MAXOLLEN;
$p{LAPTEMP}    = $p{LAPTEMP}    || $STITCH;
$p{MINLAPLEN}  = $p{MINLAPLEN}  || $MINOLLAP;
$p{TMTOL}      = $p{TMTOL}      || $TMTOL;
$p{POOLSIZE}   = $p{POOLSIZE}   || $POOLSIZE;
$p{MAXPOOLNUM} = $p{MAXPOOLNUM} || $MAXPOOLNUM;

croak "\n DNA_ERROR: oligo size is outside of allowable range.\n"
  if ($p{OLILEN} < $p{MINOLILEN} || $p{OLILEN} > $p{MAXOLILEN});


################################################################################
################################# CONFIGURING ##################################
################################################################################
my @bbs = $design->get_constructs(-kind => 'building_block');


################################################################################
################################## CHOPPING  ###################################
################################################################################
foreach my $bb (@bbs)
{
  my $ps = {
    building_block   => $bb,
    oligo_len_min    => $p{MINOLILEN},
    oligo_len_max    => $p{MAXOLILEN},
    oligo_len        => $p{OLILEN},
    overlap_tm       => $p{LAPTEMP},
    overlap_len_min  => $p{MINLAPLEN},
    tm_tolerance     => $p{TMTOL},
    pool_size        => $p{POOLSIZE},
    pool_num_max     => $p{MAXPOOLNUM},
  };
  my $return = chop_oligos($ps);
  my $num = $return->{constructs_made};
  if (! $num || $num < 1)
  {
    $bb->method('magic');
    $design->weep(-construct => $bb, -note => $return->{error});
    next;
  }
  if (defined $return->{warning})
  {
    $design->warn(-construct => $bb, -note => $return->{warning});
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
################################# SUBROUTINES  #################################
################################################################################

sub chop_oligos
{
  my ($ps) = @_;
  my $v = 1;
  my $bbcon = $ps->{building_block};
  my $laptemp = $ps->{overlap_tm};
  my $laplenmin = $ps->{overlap_len_min};
  my $tmtol = $ps->{tm_tolerance};
  my $poolsize = $ps->{pool_size};
  my $poolnummax = $ps->{pool_num_max};
  my $return = {error => undef, warning => undef, constructs => []};
  my @constructs = ();

  # Add prefixes and/or suffixes, if necessary; create a seqobj for BB seq for
  # filtering overlaps later
  #
  my $bbname    = $bbcon->id();
  my $bbseq     = uc $bbcon->sequence();
  my $bblen     = length $bbseq;
  my $bbseqobj  = Bio::Seq->new(-id => $bbname, -seq => $bbseq);
  print {$OUTH} "Chopping $bbname ($bblen bp)\n" if ($v);

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
    my $tm = $GD->melt(-sequence => $ol);

    if ($tm < $laptemp - $tmtol)
    {
      $rig++;
      next;
    }
    elsif ($tm > $laptemp + $tmtol)
    {
      $lef++;
      $rig = $laplenmin  ? $lef + $laplenmin  : $lef + 7;
      next;
    }

    my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
    my $tma     = Bio::Annotation::SimpleValue->new(-value => $tm);
    $lapobj->add_Annotation('Tm', $tma);
    my $lefta   = Bio::Annotation::SimpleValue->new(-value => $lef + 1);
    $lapobj->add_Annotation('left', $lefta);
    my $righta  = Bio::Annotation::SimpleValue->new(-value => $rig + 1);
    $lapobj->add_Annotation('right', $righta);
    push @{$laps}, $lapobj;
    $lef++;
    $id++;
    $rig = $laplenmin  ? $lef + $laplenmin  : $lef + 7;
  }
  my $subreturn = $DNA->choose_overlaps(
    -construct_sequence => $bbseqobj,
    -overlap_list       => $laps,
    -overlap_len_min    => $laplenmin,
    -span_len           => $ps->{oligo_len},
    -span_len_min       => $ps->{oligo_len_min},
    -span_len_max       => $ps->{oligo_len_max},
  );
  my @chosenlaps = @{$subreturn->{overlaps}};
  if (scalar @chosenlaps < 1)
  {
    $return->{error} = $subreturn->{error};
    return $return;
  }
  
  #Make the oligos
  my $y = 1;
  my $lpos = 1;
  my $divisor = $poolnummax || 2;
  my $count = scalar @chosenlaps;
  my $maxol = $poolnummax  ? $poolnummax * ($poolsize - 1) : undef;
  if ($maxol && ( $count + 1 ) > $maxol)
  {
    print {$OUTH} "DNAWARNING: $bbname has more than $poolnummax pools.\n";
    $return->{warning} .= "More than $poolnummax PCA pools";
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
  my $olnum = 1;
  my $poolatts = {number => 1, method => 'PCA'};
  $bbcon->method('fusion');
  my @universals = ();

  ##Make the small forward
  my $frig = $laplenmin;
  $frig++ while ($GD->melt((substr $bbseq, 0, $frig)) < $laptemp - $tmtol);
  my $olseq = substr $bbseq, 0, $frig;
  my $olname = $bbname . q{.} . $DNA->pad($olnum, 2);
  my $oligo = make_oligo($olname, $olseq, q{F});
  push @constructs, $oligo;
  push @universals, $oligo;
  $olnum++;

  while (scalar @chosenlaps)
  {
    # Make a forward oligo
    my $lap = shift @chosenlaps;
    $olseq = substr $bbseq, $lpos-1, $lap->[3] - $lpos + 1;
    $olname = $bbname . q{.} . $DNA->pad($olnum, 2);
    $oligo = make_oligo($olname, $olseq, q{F});
    push @constructs, $oligo;
    $olnum++;

    # Make the pool spanners if necessary, also incrementing poolnumber
    if ($flag && $y % $div == 0)
    {
      $olseq = $GD->complement(-sequence => $lap->[1], -reverse => 1);
      $olname = $bbname . q{.} . $DNA->pad($olnum, 2);
      $oligo = make_oligo($olname, $olseq, q{R});
      push @constructs, $oligo;
      $olnum++;

      my $pool = $bbcon->add_pool($poolatts);
      $pool->add_subconstructs(\@constructs);
      @constructs = ();
      $poolatts->{number}++;

      $olseq = $lap->[1];
      $olname = $bbname . q{.} . $DNA->pad($olnum, 2);
      $oligo = make_oligo($olname, $olseq, q{F});
      push @constructs, $oligo;
      $olnum++;
    }
    $lpos = $lap->[2];
    $y++;
  }

  #Make the lagging oligo
  $olseq = substr $bbseq, $lpos - 1;
  $olname = $bbname . q{.} . $DNA->pad($olnum, 2);
  $oligo = make_oligo($olname, $olseq, q{F});
  push @constructs, $oligo;
  $olnum++;

  ##Make the small reverse
  my $flef = $laplenmin;
  $flef++ while ($GD->melt((substr $bbseq, -$flef)) < $laptemp - $tmtol);
  $olseq = substr $bbseq, -$flef;
  $olseq = $GD->complement(-sequence => $olseq, -reverse => 1);
  $olname = $bbname . q{.} . $DNA->pad($olnum, 2);
  $oligo = make_oligo($olname, $olseq, q{R});
  push @constructs, $oligo;
  push @universals, $oligo;
  $olnum++;
  my $pool = $bbcon->add_pool($poolatts);
  $pool->add_subconstructs(\@constructs);
  $bbcon->add_subconstructs(\@universals);
  $return->{constructs_made} = $olnum;
  return $return;
}

sub make_oligo
{
  my ($name, $sequence, $orient) = @_;
  my $oligocon = $design->create_construct({
    id       => $name,
    kind     => 'assembly_oligo',
    method   => 'order',
    orient   => $orient,
    sequence => $sequence,
  });
  return $oligocon;
}


__END__

=head1 NAME

  DNAssemble_160_Chop.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

    This script will break each nucleotide sequence it is given into a set of
    overlapping assembly oligonucleotides. It uses EMBOSS palindrome to check
    for potential secondary structures.

    Output will be named by the identification of the part, and will be tagged
    with the _OL suffix. Default values are assumed for every parameter not
    provided.

    You must provide either a target oligo length or a minimum maximum range.
    You can provide both.

=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing building block sequences.


  Optional arguments:
    --output : a filepath for dumping output

    --minlaplen : minimum length of overlap allowed (default 10)

    --laptemp : Tm in temperature degrees C of oligo overlaps (default 70)

    --olilen : the length in bp of assembly oligos (default 150)

    --minolilen : minimum length of assembly oligos permitted (default 45)

    --maxolilen : maximum length of assembly oligos permitted (default 200)

    -t,  --tolerance : amount of +/- variation allowed in Tm (default 2.5)

    -p,  --poolsize : oligos will be pooled; GD will make bridging primers
            between pools of the specified size

    --maxpoolnum : the maximum number of pools to create.

    -f,  --format : default genbank

    -v,  --verbose : show deliberations happening

    -h,  --help: display this message

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