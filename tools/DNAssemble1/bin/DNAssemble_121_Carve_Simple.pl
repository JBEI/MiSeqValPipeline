#!/usr/bin/env perl

use Bio::DNAssemble;
use Bio::GeneDesign;
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
my $DNAV = 'DNAssemble_121_Carve_Simple_' . $VERSION;
my $DNAS = '_121';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'             => \$p{INPUT},
      'cloningvector=s'     => \$p{CLONVEC},
      'output:s'            => \$p{OUTPUT},
      'length:i'            => \$p{TARBBLEN},
      'maxlength:i'         => \$p{MAXBBLEN},
      'minlength:i'         => \$p{MINBBLEN},
      'lap:i'               => \$p{TARBBLAP},
      'stitch:i'            => \$p{STITCH},
      'verbose:s'           => \$p{VERBOSE},
      'logfile:s'           => \$p{LOGPATH},
      'help'                => \$p{HELP}
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

#The cloning vector must be provided.
croak "\n DNA_ERROR: You must name a cloning vector.\n" if (! $p{CLONVEC});
my $clonvec = $DNA->load_vector(-name => $p{CLONVEC});

Readonly my $TARBBLEN      => 1000;
Readonly my $MAXBBLEN      => 1023;
Readonly my $MINBBLEN      =>  140;
Readonly my $TARBBLAP      =>   40;
Readonly my $MINBBLAP      =>   25;
Readonly my $MAXBBLAP      =>   50;
Readonly my $STITCH        =>   70;
Readonly my $INTBBLIMIT    =>    5;
Readonly my $BBGROUPSIZE   => $INTBBLIMIT - 1;
Readonly my $BANDSIZERANGE =>  500;
Readonly my $BUFFERLIMIT   =>   50;
Readonly my $SEARCHLIMIT   =>  500;

$p{TARBBLEN} = $p{TARBBLEN} || $TARBBLEN;
$p{MAXBBLEN} = $p{MAXBBLEN} || $MAXBBLEN;
$p{MINBBLEN} = $p{MINBBLEN} || $MINBBLEN;
$p{TARBBLAP} = $p{TARBBLAP} || $TARBBLAP;
$p{STITCH}   = $p{STITCH}   || $STITCH;

croak "\n DNA_ERROR: building block size is outside of allowable range.\n"
  if ($p{TARBBLEN} < $p{MINBBLEN} || $p{TARBBLEN} > $p{MAXBBLEN});

croak "\n DNA_ERROR: chewback overlap is too small.\n"
  if ($p{TARBBLAP} < $MINBBLAP);


################################################################################
################################# CONFIGURING ##################################
################################################################################
#Set up vectors, if necessary
my $cvname = $clonvec->name;
my (%vseqs, %vstats) = ((), ());

my $newstart = $clonvec->chew3loc();
my $tplen = length $clonvec->chew3;
my $fplen = length $clonvec->chew5;
my $addlen = $tplen + $fplen;
my $oldseq = $clonvec->seq();
my $newseq = substr $oldseq, $newstart - 2 + $tplen;
$newseq .= substr $oldseq, 0, $newstart - 1 - $fplen;
$vseqs{$cvname} = $newseq;


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
  my $username = $chunk->uid();
  $username = $username ne q{}  ? $username : $chunkname;

  # Attempt to carve the chunk
  #
  #
  print {$OUTH} "Working on $chunkname ($chlen bp)...\n";
  my $ps = {
    chunkseq   => $chseq,
    maxbblen   => $p{MAXBBLEN},
    minbblen   => $p{MINBBLEN},
    tarbblen   => $p{TARBBLEN},
    overlap    => $p{TARBBLAP},
    addlen     => $addlen,
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
    next;
  }

  # Address how building blocks will be assembled
  #
  #
  my $sublet = 'A';
  my $params = {
      clonvec    => $clonvec,
      stitchtemp => $p{STITCH}
  };

  # If there are at least 2 building blocks for this construct,
  # it can be made as one pool.
  if ($num > 1)
  {
    my $fragname = $chunkname . q{.} . $sublet;

    my $firstbb = shift @{$BBS};
    my $lastbb = pop @{$BBS};

    $design->delete_construct($chunk);
    my $iatts = {
      id          => $fragname,
      deliverable => 'true',
      kind        => 'intermediate',
      method      => 'chewback',
      vector      => $cvname,
      sequence    => $clonvec->chew5 . $chseq . $clonvec->chew3,
      uid         => $username
    };
    my $intermediate = $design->create_construct($iatts);
    my $patts = {number => 1, method => 'amplification'};
    my $firstname = $fragname . q{.} . $DNA->pad($patts->{number}, 2);
    my %fparams = %{$params};
    $fparams{bbhsh} = $firstbb;
    $fparams{bbseq} = bb_seq($firstbb, $chseq);
    $fparams{isfirst} = 1;
    $fparams{islast} = 0;
    $fparams{isedge} = 1;
    $fparams{name} = $firstname;
    my ($firstcon, $fslcon, $fsrcon) = make_bbcon(\%fparams);
    my $pool = $intermediate->add_pool($patts);
    $pool->add_subconstructs([$fslcon, $fsrcon, $firstcon]);
    $patts->{number}++;

    foreach my $bb (@{$BBS})
    {
      my $bbname = $fragname . q{.} . $DNA->pad($patts->{number}, 2);
      my %mparams = %{$params};
      $mparams{bbhsh} = $bb;
      $mparams{bbseq} = bb_seq($bb, $chseq);
      $mparams{isfirst} = 0;
      $mparams{islast} = 0;
      $mparams{isedge} = 0;
      $mparams{name} = $bbname;
      my ($bbcon, $slcon, $srcon) = make_bbcon(\%mparams);
      $pool = $intermediate->add_pool($patts);
      $pool->add_subconstructs([$slcon, $srcon, $bbcon]);
      $patts->{number}++;
    }

    my $lastname = $fragname . q{.} . $DNA->pad($patts->{number}, 2);
    my %lparams = %{$params};
    $lparams{bbhsh} = $lastbb;
    $lparams{bbseq} = bb_seq($lastbb, $chseq);
    $lparams{isfirst} = 0;
    $lparams{islast} = 1;
    $lparams{isedge} = 1;
    $lparams{name} = $lastname;
    my ($lastcon, $lslcon, $lsrcon) = make_bbcon(\%lparams);
    $pool = $intermediate->add_pool($patts);
    $pool->add_subconstructs([$lslcon, $lsrcon, $lastcon]);
    $patts->{number}++;
  }


  # If there is only 1 building block for this construct, its a lot simpler.
  # Just prepend and postpend the vector sequences.
  else
  {
    $design->delete_construct($chunk);
    my $fragname = $chunkname . q{.} . $sublet;
    my $iatts = {
      id          => $fragname,
      deliverable => 'true',
      kind        => 'intermediate',
      method      => 'chewback',
      uid         => $username,
      sequence    => $clonvec->chew5 . $chseq . $clonvec->chew3,
      vector      => $cvname
    };
    my $intermediate = $design->create_construct($iatts);
    my $bb = shift @{$BBS};
    my $bbname = $fragname . q{.} . $DNA->pad(1, 2);
    my %oparams = %{$params};
    $oparams{bbhsh} = $bb;
    $oparams{bbseq} = bb_seq($bb, $chseq);
    $oparams{isfirst} = 1;
    $oparams{islast} = 1;
    $oparams{isedge} = 1;
    $oparams{name} = $bbname;
    my ($bbcon, $slcon, $srcon) = make_bbcon(\%oparams);
    $intermediate->add_pool({
      number => 1,
      method => 'gblock',
      subconstruct => [$bbcon->id, $slcon->id, $srcon->id],
    });
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

  my $minbblap = $MINBBLAP;
  my $maxbblap = $MAXBBLAP;
  my $tarbblap = $ps->{overlap};
  my $tarbblen = $ps->{tarbblen};
  $tarbblen = $chlen < $tarbblen  ? $chlen  : $tarbblen;
  my $maxbblen = $ps->{maxbblen};
  my $minbblen = $ps->{minbblen};
  my $addlen   = $ps->{addlen};

  ## Decide what size bbs to go for and how many
  ## Adjust the target building block size so as to avoid outliers.
  my $tarnum = ($chlen + $addlen) > $maxbblen
              ? ceil(($chlen + $addlen) / ($tarbblen - $tarbblap))
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
      -span_len_min       => $ps->{minbblen},
      -span_len_max       => $ps->{maxbblen},
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
  my $cvector = $ps->{clonvec};
  my ($first, $lasty, $edge) = ($ps->{isfirst}, $ps->{islast}, $ps->{isedge});
  my $stemp = $ps->{stitchtemp};
  my ($bbhsh, $bbseq, $name) = ($ps->{bbhsh}, $ps->{bbseq}, $ps->{name});
  my $bblen = length $bbseq;
  my $offset = 0;
  my $toss = undef;

  # Begin processing sequence; make insert stitches
  #
  my ($lprimer, $rprimer) = $GD->make_amplification_primers(
      -sequence    => $bbseq,
      -temperature => $stemp,
  );

  # Add cloning vector sequence if this is an edge construct
  #
  if ($edge)
  {
    if ($first)
    {
      my $prefix = $cvector->chew5();
      $bbseq = $prefix . $bbseq;
      $offset += length $prefix;
      ($lprimer, $toss) = $GD->make_amplification_primers(
        -sequence    => $bbseq,
        -temperature => $stemp,
      );
    }
    if ($lasty)
    {
      $bbseq .= $cvector->chew3;
      ($toss, $rprimer) = $GD->make_amplification_primers(
        -sequence    => $bbseq,
        -temperature => $stemp,
      );
    }
  }

  # Create and annotate the constructs
  #
  my $buildingb = $design->create_construct({
    id        => $name,
    kind      => 'building_block',
    method    => 'gblock',
    istart    => $offset,
    ilen      => $bblen,
    sequence  => $bbseq,
    vector    => $cvector->{name},
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

    -le, --length : the length in bp of building blocks
            Default is 1000 bp

    -ma, --maxlength : the maximum length a building block is allowed to be.
            Default is 1023 bp

    -mi, --minlength : the minimum length a building block is allowed to be.
            Default is 200 bp

    -la, --lap : the target overlap between building blocks.
            Defaults is 40 bp

    -st, --stitch : The target melting temperature of assembly primers for
            the assembly of building blocks.
            Default is 70 degrees C

    -v,  --verbose : show deliberations happening
            Default is 0

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