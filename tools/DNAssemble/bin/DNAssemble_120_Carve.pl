#!/usr/bin/env perl

use Bio::GeneDesign;
use Bio::DNAssemble;
use Bio::DNAssemble::Digestion;
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
my $DNAV = 'DNAssemble_120_Carve_' . $VERSION;
my $DNAS = '_120';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'             => \$p{INPUT},
      'cloningvector=s'     => \$p{CLONVEC},
      'destinationvector=s' => \$p{DESTVEC},
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
my $RES = $GD->set_restriction_enzymes();

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

#The cloning vector and destination vector must be provided and different
croak "\n DNA_ERROR: You must name a cloning vector.\n" if (! $p{CLONVEC});
croak "\n DNA_ERROR: You must name a destination vector.\n" if (! $p{DESTVEC});
my $clonvec = $DNA->load_vector(-name => $p{CLONVEC});
my $destvec = $DNA->load_vector(-name => $p{DESTVEC});

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
#Set up two restriction enzyme sets for filtering enzymes through
my $ORE = Bio::GeneDesign->new();
$ORE->set_restriction_enzymes(-enzyme_set => 'outside');
my $ARE = Bio::GeneDesign->new();
$ARE->set_restriction_enzymes(-enzyme_set => 'blunts');

#Set up vectors, if necessary; note which enzymes should be excluded
my $cvname = $clonvec->name;
my $dvname = $destvec->name;
my (%excludes, %fpnos, %tpnos, %vseqs, %vstats) = ((), (), (), (), ());
$excludes{$_}++ foreach ($clonvec->enzyme_list);
$excludes{$_}++ foreach ($destvec->enzyme_list);
my @enzes = keys %excludes;

my $newstart = $clonvec->chew3loc();
my $tplen = length $clonvec->chew3;
my $fplen = length $clonvec->chew5;
my $oldseq = $clonvec->seq();
my $newseq = substr $oldseq, $newstart - 2 + $tplen;
$newseq .= substr $oldseq, 0, $newstart - 1 - $fplen;
$vseqs{$cvname} = $newseq;
$vstats{$cvname} = $ORE->restriction_status(-sequence => $newseq);
my $stat5 = $ARE->restriction_status(-sequence => $destvec->chew5);
%fpnos = map {$_ => 1} grep {$stat5->{$_} != 0} keys %{$stat5};
my $stat3 = $ARE->restriction_status(-sequence => $destvec->chew3);
%tpnos = map {$_ => 1} grep {$stat3->{$_} != 0} keys %{$stat3};
my @outsides = values %{$ORE->enzyme_set};
$ARE->add_to_enzyme_set(-enzymes => \@outsides);


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
  my $chseq = $chunk->sequence();
  my $chlen = length $chseq;
  my $username = $chunk->uid();
  $username = defined $username ? $username : $chunkname;

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
    excludes   => \@enzes,
    fpexcludes => \%fpnos,
    tpexcludes => \%tpnos,
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
      destvec    => $destvec,
      stitchtemp => $p{STITCH}
  };

  # If there are more than $INTBBLIMIT building blocks for this construct, it
  # will have to be made into intermediates that will then be put together. So,
  # all building blocks will be given stitching primers that amplify cloning
  # vector sequence on the outside of groups, and intermediates will be put
  # back into the same cloning vector. The intermediates will be cut out of
  # their cloning vector by restriction enzymes (with care taken that the bands
  # can be cut out without ambiguity) and then cloned directly into the
  # destination vector. Two vectors are required for this to work properly.
  if ($num > $INTBBLIMIT)
  {
    if ($dvname eq $cvname)
    {
      print {$OUTH} "\t$chunkname will require two vectors and is being skipped.\n";
      $design->weep(-construct=> $chunk, -note => 'two vectors required');
      next;
    }
    my @proced;
    my $f = floor $num / $BBGROUPSIZE;
    my $t = ($BBGROUPSIZE - ($num % $BBGROUPSIZE)) % $BBGROUPSIZE;
    $f = $f - ($t - 1) if ($t > 1);
    my $forder = $BBGROUPSIZE x $f;
    my $remainder = $BBGROUPSIZE - 1;
    my $torder = $remainder x $t;
    my $order = $forder . $torder;
    $chunk->method('chewback');
    my $chatts = {number => 1, method => 'enzymatic_digestion'};
    foreach my $size (split q{}, $order)
    {
      my @subset;
      for (my $x = 1; $x <= $size; $x++)
      {
        my $bb = shift @{$BBS};
        push @subset, $bb;
      }
      my $fragname = $chunkname . q{.} . $sublet;
      $sublet++;
      my $iatts = {
        id     => $fragname,
        kind   => 'intermediate',
        method => 'chewback'
      };
      my $intermediate = $design->create_construct($iatts);

      my $firstbb = shift @subset;
      my $lastbb = pop @subset;
      my $start = $firstbb->{start};
      my $end = $lastbb->{end};
      my ($fpenz, $tpenz, $openz) = (undef, undef, undef);
      my $subseq = substr $chseq, $start - 1, $end - $start + 1;
      my $patts = {number => 1, method => 'fragamp'};
      my $firstname = $fragname . q{.} . $DNA->pad($patts->{number}, 2);
      my %fparams = %{$params};
      $fparams{bbhsh} = $firstbb;
      $fparams{bbseq} = bb_seq($firstbb, $chseq);
      $fparams{isfirst} = 1;
      $fparams{islast} = 0;
      $fparams{isedge} = 0;
      $fparams{name} = $firstname;

      if (! scalar @proced)
      {
        $fparams{isedge} = 1;
        my $inenz = $ARE->enzyme_set->{$lastbb->{enz_t}};
        my $prepend = $destvec  ? $destvec->chew5  . $subseq  : $subseq;
        my $outenzarr = pick_enzyme($prepend, $inenz, $clonvec, undef);
        if (! $outenzarr)
        {
          my $reason = 'missing enzyme for ' . $fragname;
          $design->warn(-construct => $chunk, -note => $reason);
          $chunk->method('magic');
        }
        else
        {
          $fparams{outenzseq} = $outenzarr->[1];
          $fpenz = $outenzarr->[0];
        }
      }
      else
      {
        $fparams{doinenz} = 1;
        $fpenz = $firstbb->{enz_f};
      }
      my ($firstcon, $fslcon, $fsrcon, $fpref, $fsuff) = make_bbcon(\%fparams);
      my $pool = $intermediate->add_pool($patts);
      $pool->add_subconstructs([$fslcon, $fsrcon, $firstcon]);
      $patts->{number}++;
      push @proced, $firstcon;

      foreach my $bb (@subset)
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
        push @proced, $bbcon;
      }

      my $lastname = $fragname . q{.} . $DNA->pad($patts->{number}, 2);
      my %lparams = %{$params};
      $lparams{bbhsh} = $lastbb;
      $lparams{bbseq} = bb_seq($lastbb, $chseq);
      $lparams{isfirst} = 0;
      $lparams{islast} = 1;
      $lparams{isedge} = 0;
      $lparams{name} = $lastname;
      if (! scalar @{$BBS})
      {
        $lparams{isedge} = 1;
        my $inenz = $ARE->enzyme_set->{$firstbb->{enz_f}};
        my $postpend = $destvec  ? $subseq . $destvec->chew3  : $subseq;
        my $outenzarr = pick_enzyme($postpend, $inenz, undef, $clonvec);
        if (! $outenzarr)
        {
          my $reason = 'missing enzyme for ' . $fragname;
          $design->warn(-construct=> $chunk, -note => $reason);
          $chunk->method('magic');
        }
        else
        {
          $lparams{outenzseq} = $outenzarr->[1];
          $tpenz = $outenzarr->[0];
        }
      }
      else
      {
        $lparams{doinenz} = 1;
        $tpenz = $lastbb->{enz_t};
      }
      my ($lastcon, $lslcon, $lsrcon, $lpref, $lsuff) = make_bbcon(\%lparams);
      $pool = $intermediate->add_pool($patts);
      $pool->add_subconstructs([$lslcon, $lsrcon, $lastcon]);
      $patts->{number}++;
      push @proced, $lastcon;

      $intermediate->vector($cvname);
      my $intseq = $fpref . $subseq . $lsuff;
      $intermediate->sequence($intseq);
      $chatts->{subconstruct} = [$intermediate->id];
      my $cconpool = $chunk->add_pool($chatts);
      $chatts->{number}++;

      # Define enzyme reaction conditions
      # Check the cloning vector sequence to make sure the cleaving enzymes
      # will distinguish intermediate band from the vector band
      # if not, pick a third compatible enzyme
      my $seqlen = length $intseq;
      my $wholeseq = $intseq . $vseqs{$cvname};

      my $renzf = $RES->{$fpenz};
      my $renzt = $RES->{$tpenz};
      my $renzttemp = $renzt->temp();
      my $renzftemp = $renzf->temp();
      my $fbuffers = $renzf->common_buffers($renzf);
      my $tbuffers = $renzt->common_buffers($renzt);
      my ($rxna, $rxnb) = (undef, undef);
      my $bestbuffer = $renzf->acceptable_buffer($renzt, $BUFFERLIMIT);
      if (defined $bestbuffer && $renzttemp == $renzftemp)
      {
        $rxna = $cconpool->add_reaction({
          temperature => $renzttemp,
          buffer      => $bestbuffer
        });
        my $funits = $renzf->units($bestbuffer, $wholeseq);
        if ($fpenz ne $tpenz)
        {
          my $tunits = $renzt->units($bestbuffer, $wholeseq);
          $rxna->add_enzyme({name => $fpenz, units => $funits, type => '5p'});
          $rxna->add_enzyme({name => $tpenz, units => $tunits, type => '3p'});
        }
        else
        {
          $rxna->add_enzyme({name => $fpenz, units => $funits, type => '53p'});
        }
      }
      else
      {
        my $fusebuff = $fbuffers->[0];
        $rxna = $cconpool->add_reaction({
          temperature => $renzftemp,
          buffer      => $fusebuff
        });
        my $funits = $renzf->units($fusebuff, $wholeseq);
        $rxna->add_enzyme({name => $fpenz, units => $funits, type => '5p'});
        
        my $tusebuff = $tbuffers->[0];
        $rxnb = $cconpool->add_reaction({
          temperature  =>  $renzttemp,
          buffer       =>  $tusebuff,
        });
        my $tunits = $renzt->units($tusebuff, $wholeseq);
        $rxnb->add_enzyme({name => $tpenz, units => $tunits, type => '3p'});
      }

      my $bands = bands([$renzf, $renzt], $vseqs{$cvname});
      my $vprob = find_obscuring_band($seqlen, $bands, $BANDSIZERANGE);
      if (defined $vprob)
      {
        my $buffa = $rxna->buffer;
        my $tempa = $rxna->temperature;
        my $istats = $GD->restriction_status(-sequence => $intseq);
        my $pstats = $GD->restriction_status(-sequence => $vprob);
        my @cands = grep {$istats->{$_} == 0 } keys %{$RES};
        @cands = grep {$pstats->{$_} > 0} @cands;
        @cands = sort {$RES->{$a}->score <=> $RES->{$b}->score} @cands;
        my @lastditch;
        foreach my $enz (@cands)
        {
          my $renzo = $RES->{$enz};
          my $obands = bands([$renzo], $vprob);
          my $sameflag = find_obscuring_band($seqlen, $obands, $BANDSIZERANGE);
          next if defined $sameflag;

          my $otemp = $renzo->temp();
          my $obuffs = $renzo->buffers();
          if (exists $obuffs->{$buffa}
           && $obuffs->{$buffa} >= $BUFFERLIMIT
           && $otemp == $tempa)
          {
            $openz = $enz;
            my $ounits = $renzo->units($buffa, $wholeseq);
            $rxna->add_enzyme({name => $openz, units => $ounits, type => 'vprob'});
            last;
          }
          elsif (defined $rxnb)
          {
            my $buffb = $rxnb->buffer;
            my $tempb = $rxnb->temperature;
            if (exists $obuffs->{$buffb}
             && $obuffs->{$buffb} >= $BUFFERLIMIT
             && $otemp == $tempb)
            {
              $openz = $enz;
              my $ounits = $renzo->units($buffb, $wholeseq);
              $rxnb->add_enzyme({name => $openz, units => $ounits, type => 'vprob'});
              last;
            }
          }
          push @lastditch, $enz;
        }
        if (! defined $openz && ! defined $rxnb && scalar @lastditch)
        {
          $openz = $lastditch[0];
          my $renzo = $RES->{$openz};
          my $obuff = $renzo->acceptable_buffer($renzo);
          my $ounits = $renzo->units($obuff, $wholeseq);
          $rxnb = Bio::DNAssemble::Digestion->new(
            -temp   => $renzo->temp(),
            -buffer => $obuff,
          );
          $rxnb->add_enzyme({name => $openz, units => $ounits, type => 'vprob'});
        }
        if (! defined $openz)
        {
          my $warnmsg = 'obscured digestion band on ' . $fragname;
          $design->warn(-construct=> $chunk, -note => $warnmsg);
        }
      }
    }
    $chunk->vector($dvname);
    $chunk->deliverable('true');
    $chunk->kind('chunk');
    
    my $newchseq = $destvec->chew5 . $chseq . $destvec->chew3;
    $chunk->sequence($newchseq);
  }


  # If there are between 2 and $INTBBLIMIT building blocks for this construct,
  # it can be made as one pool. Each building block will be created in the
  # cloning vector, and stitching primers will amplify destination vector
  # sequence on the outside of the group. The intermediate can therefore be
  # cloned directly into the destination vector.
  elsif ($num > 1)
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
      vector      => $dvname,
      sequence    => $destvec->chew5 . $chseq . $destvec->chew3,
      uid         => $username
    };
    my $intermediate = $design->create_construct($iatts);
    my $patts = {number => 1, method => 'fragamp'};
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
      sequence    => $destvec->chew5 . $chseq . $destvec->chew3,
      vector      => $dvname
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
    my ($bbcon, $slcon, $srcon, $pref, $suff) = make_bbcon(\%oparams);
    $intermediate->add_pool({
      number => 1,
      method => 'fragamp',
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
  my $excludes = $ps->{excludes};
  my $fpexcludes = $ps->{fpexcludes};
  my $tpexcludes = $ps->{tpexcludes};

  ## Decide what size bbs to go for and how many
  ## Adjust the target building block size so as to avoid outliers.
  my $tarnum = $chlen > $maxbblen
              ? ceil($chlen / ($tarbblen - $tarbblap))
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
    ## Decide what restriction enzyme half sites to use
    $SGD->set_restriction_enzymes(-enzyme_set => 'blunts');
    $SGD->remove_from_enzyme_set(-enzymes => $excludes) if ($excludes);
    my $SRES = $SGD->enzyme_set;
    my $objstats = $SGD->restriction_status(-sequence => $chseq);
    my @objcandidates = sort {$SRES->{$a}->score <=> $SRES->{$b}->score}
                        grep {$objstats->{$_} == 0}
                        keys %{$objstats};
    print {$OUTH} "\tCH CANDIDATES: @objcandidates\n" if ($v);
    if (! scalar @objcandidates)
    {
      $return->{error} = 'no blunt enzyme candidates for partitioning';
      return $return;
    }

    ## Create a prefix tree for each end of a bb overlap
    ## The 5' half of an RE site goes in ftree, the 3' half goes in ttree
    my $ftree = Bio::GeneDesign::PrefixTree->new();
    my $ttree = Bio::GeneDesign::PrefixTree->new();
    my %seeks;
    foreach my $cand (@objcandidates)
    {
      my $seq = uc $SRES->{$cand}->seq;
      my $qes = $GD->complement(-sequence => $seq, -reverse => 1);
      my $half = length($seq) / 2;
      $seeks{$cand} = {len => $half};
      my $a = substr $seq, 0, $half;
      my $b = substr $seq, $half;
      my $fnucs = $GD->ambiguous_transcription($b);
      foreach my $fnuc (@{$fnucs})
      {
        $ftree->add_prefix($fnuc, $cand);
        $seeks{$cand}->{$fnuc} = $a;
      }
      my $tnucs = $GD->ambiguous_transcription($a);
      foreach my $tnuc (@{$tnucs})
      {
        $ttree->add_prefix($tnuc, $cand);
        $seeks{$cand}->{$tnuc} = $b;
      }
      if ($seq ne $qes)
      {
        my $c = substr $qes, 0, $half;
        my $d = substr $qes, $half;
        $fnucs = $GD->ambiguous_transcription($d);
        foreach my $fnuc (@{$fnucs})
        {
          $ftree->add_prefix($fnuc, $cand);
          $seeks{$cand}->{$fnuc} = $c;
        }
        $tnucs = $GD->ambiguous_transcription($c);
        foreach my $tnuc (@{$tnucs})
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
    my $mindist = $minbblen - $maxbblap;
    my %ranges;
    foreach my $fhit (@{$fref})
    {
      my $lef = $fhit->[1];
      next if ($lef < $mindist);
      my $enz_f = $fhit->[0];
      my $zne_f = $seeks{$enz_f}->{$fhit->[2]};
      foreach my $thit (@{$tref})
      {
        my $enz_t = $thit->[0];
        my $rig = $thit->[1] + $seeks{$enz_t}->{len};
        my $dist = $rig - $lef;
        next if ($dist < $minbblap || $dist > $maxbblap);
        next if ($rig > $chlen - $mindist);
        my $zne_t = $seeks{$enz_t}->{$thit->[2]};
        my $ol = substr $chseq, $lef, $dist;
        my $lapobj = Bio::Seq->new( -seq => $ol, -id => $id);
        $anns{$id} = {seq => $ol, left => $lef + 1, right => $rig,
                      enz_f => $enz_f, zne_f => $zne_f,
                      enz_t => $enz_t, zne_t => $zne_t};
        push @{$laps}, $lapobj;
        $ranges{$lef} = $rig if (! $ranges{$lef} || $ranges{$lef} < $rig);
        $id++;
      }
    }
    print {$OUTH} "\t", scalar(@{$laps}), " overlaps found\n" if ($v);
    if (! scalar @{$laps})
    {
      $return->{error} = 'no usable overlaps for partitioning';
      return $return;
    }

    #Filter overlaps for palindromes and mispriming
    #If vmatch or blast are not available, filter for uniqueness
    if ($DNA->EMBOSS)
    {
      $laps = $DNA->filter_palindromes(-sequences => $laps);
      print {$OUTH} "\t", scalar(@{$laps}), " overlaps palindrome free\n" if ($v);
    }
    if ($DNA->vmatch)
    {
      $laps = $DNA->filter_vmatch(-sequences => $laps, -parent => $chseqobj);
      print {$OUTH} "\t", scalar(@{$laps}), " overlaps informative\n" if ($v);
    }
    elsif ($DNA->BLAST)
    {
      $laps = $DNA->filter_blast(-sequences => $laps, -parent => $chseqobj);
      print {$OUTH} "\t", scalar(@{$laps}), " overlaps informative\n" if ($v);
    }
    else
    {
      $laps = $DNA->filter_uniqueness(-sequences => $laps);
      print {$OUTH} "\t", scalar(@{$laps}), " overlaps unique\n" if ($v);
    }
    my @overlaps = map {$_->id} @{$laps};

    #Here we should assess spread and decide if this is possible
    my $mask = 'o' x $chlen;
    foreach my $lef (keys %ranges)
    {
      my $rig = $ranges{$lef};
      my $laplen = $rig - $lef + 1;
      substr $mask, $lef - 1, $laplen, 'i' x $laplen;
    }
    my $gapsize = $maxbblen - ($minbblap * 2);
    my $gap = 'o{' . $gapsize . ',}';
    my $deadzones = {};
    while ($mask =~ /($gap)/g)
    {
      my $size = length $1;
      my $pos = pos $mask;
      $deadzones->{$pos - $size} = length $1;
    }
    my @errors;
    foreach my $pos (sort keys %{$deadzones})
    {
      my $len = $deadzones->{$pos};
      my $end = $pos + $len;
      push @errors, 'no coverage from ' . $pos . ' to ' . $end;
    }
    if (scalar @errors)
    {
      $return->{error} = join q{, }, @errors;
      return $return;
    }

    my $f = floor $tarnum / $BBGROUPSIZE;
    my $t = ($BBGROUPSIZE - ($tarnum % $BBGROUPSIZE)) % $BBGROUPSIZE;
    $f = $f - ($t - 1) if ($t > 1);
    my $remainder = $BBGROUPSIZE - 1;
    my $arr = $BBGROUPSIZE x $f . $remainder x $t;
    my $frange = $INTBBLIMIT;
    my $trange = $BBGROUPSIZE;


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
    while ($lpos < $rbound && $amideadyet < $SEARCHLIMIT)
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
        @laps = grep {$SRES->{$anns{$_}->{enz_t}}->temp == $SRES->{$lenz}->temp} @laps;
        @laps = grep {$SRES->{$anns{$_}->{enz_t}}->common_buffers($SRES->{$lenz}, 1)} @laps;
      }
      if (defined $fpexcludes && scalar @chosenlaps < $frange)
      {
        @laps = grep {! exists $fpexcludes->{$anns{$_}->{enz_t}} } @laps;
      }
      if (defined $tpexcludes && scalar @chosenlaps > $trange)
      {
        @laps = grep {! exists $tpexcludes->{$anns{$_}->{enz_f}} } @laps;
      }
      @laps = sort { $SRES->{$anns{$a}->{enz_f}}->score + $SRES->{$anns{$a}->{enz_t}}->score
                 <=> $SRES->{$anns{$b}->{enz_f}}->score + $SRES->{$anns{$b}->{enz_t}}->score
                 || abs($anns{$a}->{right} - $rtarget) <=> abs($anns{$b}->{right} - $rtarget)
                 || length($anns{$a}->{left}) <=> length($anns{$b}->{left})}
              @laps;
      if (! scalar @laps)
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
          print {$OUTH} "\t\tAdjusting from $tarbblen " if ($v);
          $tarbblen = $counter % 2 == 0
                      ? $tarbblen + $counter
                      : $tarbblen - $counter;
          print {$OUTH} "to $tarbblen " if ($v);
          if (($maxbblen && $tarbblen > $maxbblen)
           || ($minbblen && $tarbblen < $minbblen))
          {
            print {$OUTH} "\n\t\tDNAWARNING: giving up... OSCILLATED TOO FAR\n" if ($v);
            $fail++;
            $return->{error} = 'no building blocks can be made';
            return $return;
          }
          $counter++;
        }
      }
      else
      {
        my $lap = $laps[0];
        my $laplen = length($anns{$lap}->{seq});
        $lenz =  $anns{$lap}->{enz_f};
        print {$OUTH} "\tChoosing overlap $laplen bp ", $anns{$lap}->{seq}, if ($v);
        print {$OUTH} ' at ', $anns{$lap}->{left}, q{..}, $anns{$lap}->{right} if ($v);
        print {$OUTH} q{ (}, $anns{$lap}->{enz_f}, q{..}, $anns{$lap}->{enz_t}, q{)} if ($v);
        print {$OUTH} q{ (}, $anns{$lap}->{zne_f}, q{..}, $anns{$lap}->{zne_t}, ")\n" if ($v);
        push @chosenlaps, $lap;
        $rpos = $anns{$chosenlaps[-1]}->{right};
        $lpos = $anns{$chosenlaps[-1]}->{left};
        $redoflag = 0;
        if ($lpos >= $rbound)
        {
          if (($minbblen && ($chlen - $lpos) + 1 < $minbblen)
           || ($maxbblen && ($chlen - $lpos) + 1 > $maxbblen))
          {
            print {$OUTH} "\t\t\tdiscarding last choice\n" if ($v);
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
      if ($amideadyet >= $SEARCHLIMIT)
      {
        $fail++;
        print {$OUTH} "\n\t\tDNAWARNING: Repeated search over $SEARCHLIMIT times!\n" if ($v);
        next;
      }
    }
    if ($fail)
    {
      $return->{error} = 'no building blocks can be made';
      return $return;
    }


    #Make the building blocks
    $lpos = 1;
    my $lastenz = undef;
    my $lastmatch = undef;
    while (scalar @chosenlaps)
    {
      my $lap = shift @chosenlaps;
      my $atts = {start => $lpos, end => $anns{$lap}->{right}};
      $atts->{enz_f} = $lastenz if ($lastenz);
      $atts->{zne_f} = $lastmatch if ($lastmatch);
      $atts->{enz_t} = $anns{$lap}->{enz_t};
      $atts->{zne_t} = $anns{$lap}->{zne_t};
      $lastenz = $anns{$lap}->{enz_f};
      $lastmatch = $anns{$lap}->{zne_f};
      push @BBS, $atts;
      $lpos = $anns{$lap}->{left};
    }
    my $atts = {start => $lpos, end => $chlen,
      enz_f => $lastenz, zne_f => $lastmatch};
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
  my ($prefix, $suffix) = (q{}, q{});
  my ($cvector, $dvector) = ($ps->{clonvec}, $ps->{destvec});
  my ($cname, $dname) = ($cvector->name, $dvector->name);
  my ($first, $lasty, $edge) = ($ps->{isfirst}, $ps->{islast}, $ps->{isedge});
  my $stemp = $ps->{stitchtemp};
  my ($bbhsh, $bbseq, $name) = ($ps->{bbhsh}, $ps->{bbseq}, $ps->{name});
  my ($inenz, $outenzseq) = ($ps->{doinenz}, $ps->{outenzseq});
  my $bblen = length $bbseq;
  my $offset = 0;

  # Begin processing sequence; make insert stitches
  #
  my ($liprimer, $riprimer) = $GD->make_amplification_primers(
      -sequence    => $bbseq,
      -temperature => $stemp,
  );

  # Add destination vector sequence if this is an edge construct, make stitches
  #
  my ($ldprimer, $rdprimer) = (undef, undef);
  if ($edge && $cname ne $dname)
  {
    if ($first)
    {
      my $dvecseq = $dvector->chew5;
      $prefix = $dvecseq . $prefix;
      $bbseq = $dvecseq . $bbseq;
    }
    if ($lasty)
    {
      $suffix .= $dvector->chew3;
      $bbseq .= $dvector->chew3;
    }
    ($ldprimer, $rdprimer) = $GD->make_amplification_primers(
        -sequence    => $bbseq,
        -temperature => $stemp,
    );
  }

  # Add an outside cutter if this is an edge construct that requires one
  #
  if (defined $outenzseq && $edge)
  {
    if ($first)
    {
      $prefix = $outenzseq . $prefix;
      $bbseq = $outenzseq . $bbseq;
    }
    if ($lasty)
    {
      my $touseq = $GD->complement(-sequence => $outenzseq, -reverse => 1);
      $suffix .= $touseq;
      $bbseq .= $touseq;
    }
  }

  # Add the blunt finishers if this is a construct that requires them
  #
  if (defined $inenz)
  {
    if ($first)
    {
      my $enzfix = $bbhsh->{zne_f};
      $enzfix = $GD->replace_ambiguous_bases($enzfix);
      $prefix = $enzfix . $prefix;
      $bbseq = $enzfix . $bbseq;
    }
    if ($lasty)
    {
      my $enzfix = $bbhsh->{zne_t};
      $enzfix = $GD->replace_ambiguous_bases($enzfix);
      $suffix .= $enzfix;
      $bbseq .= $enzfix;
    }
  }

  # Add the cloning vector sequence and make the cloning vector stitches
  #
  $prefix = $cvector->chew5 . $prefix;
  $suffix .= $cvector->chew3;
  $offset = length $prefix;
  $bbseq = $cvector->chew5 . $bbseq . $cvector->chew3;
  my ($lcprimer, $rcprimer) = $GD->make_amplification_primers(
      -sequence    => $bbseq,
      -temperature => $stemp,
  );

  # Pick the appropriate stitching primers
  #
  my ($lprimer, $rprimer) = (undef, undef);
  if ($first)
  {
    if ($edge && $cname ne $dname)
    {
      $lprimer = $ldprimer;
    }
    else
    {
      $lprimer = $lcprimer;
    }
  }
  else
  {
    $lprimer = $liprimer;
  }
  if ($lasty)
  {
    if ($edge)
    {
      $rprimer = $rdprimer;
    }
    else
    {
      $rprimer = $rcprimer;
    }
  }
  else
  {
    $rprimer = $riprimer;
  }

  # Create and annotate the constructs
  #
  my $buildingb = $design->create_construct({
    id        => $name,
    kind      => 'building_block',
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

  return ($buildingb, $soligol, $soligor, $prefix, $suffix);
}

sub bb_seq
{
  my ($bb, $chseq) = @_;
  my $start = $bb->{start};
  my $seq = substr $chseq, $start - 1, $bb->{end} - $start + 1;
  return $seq;
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

  return q{} if (! scalar @candidates);

  my @filtered = grep {$enz->temp == $_->temp} @candidates;
  $outenz = scalar @filtered ? $filtered[0] : $outenz;

  @filtered = grep {$_->common_buffers($enz, 1)} @filtered;
  $outenz = scalar @filtered ? $filtered[0] : $outenz;

  my @refiltered = ();
  foreach my $potenz (@filtered)
  {
    if (defined $cprepend)
    {
      my $tempseq = $cprepend->chew5 . $potenz->seq;
      my $positions = $enz->positions($tempseq);
      push @refiltered, $potenz if (scalar keys %{$positions} == 0);
    }
    if (defined $cpostpend)
    {
      my $tempseq = $potenz->seq . $cpostpend->chew3;
      my $positions = $enz->positions($tempseq);
      push @refiltered, $potenz if (scalar keys %{$positions} == 0);
    }
  }
  @filtered = @refiltered;
  $outenz = scalar @filtered ? $filtered[0] : $outenz;

  my $outenzseq = $ORE->replace_ambiguous_bases($outenz->seq);
  my $rand = $ORE->random_dna(-length => $outenz->inside_cut);
  $outenzseq = $outenzseq . $rand;
  my $finalseq = defined $cprepend ?  $outenzseq . $seq : $seq . $outenzseq;
  #Check if the randomly generated sequence contains the site twice !
  my $positions = $GD->positions(-sequence => $outenzseq, -query => $finalseq);
  while (scalar keys %{$positions} != 1)
  {
    $outenzseq = $ORE->replace_ambiguous_bases($outenz->seq);
    $rand = $ORE->random_dna(-length => $outenz->inside_cut);
    $outenzseq = $outenzseq . $rand;
    $finalseq = defined $cprepend ?  $outenzseq . $seq : $seq . $outenzseq;
    $positions = $GD->positions(-sequence => $finalseq, -query => $outenz);
  }

  return [$outenz->id, $outenzseq];
}

sub bands
{
  my ($enzarr, $sequence) = @_;
  my @arr = @{$enzarr};
  my %seens = ();
  my $seqlen = length $sequence;
  my %pos = (0 => 1, $seqlen => 1);
  foreach my $enz (@arr)
  {
    next if (exists $seens{$enz});
    $seens{$enz}++;
    my $epos = $enz->positions($sequence);
    foreach my $offset (keys %{$epos})
    {
      $pos{$offset}++;
    }
  }
  my @temp = sort {$a <=> $b} keys %pos;
  my @bands;
  my $x = 0;
  while ($x < (scalar @temp) - 1)
  {
    #push @bands, $temp[$x+1] - $temp[$x] + 1 - $x;
    push @bands, substr $sequence, $temp[$x], $temp[$x+1] - $temp[$x] + 1 - $x;
    $x++;
  }
  return \@bands;
}

sub find_obscuring_band
{
  my ($target, $bandarr, $range) = @_;
  my $flag = undef;
  my @bands = @{$bandarr};
  foreach my $band (@bands)
  {
    my $blen = length $band;
    my $dist = abs $blen - $target;
    if ($dist <= $range)
    {
      $flag = $band;
      last;
    }
  }
  return $flag;
}

__END__

=head1 NAME

  DNAssemble_120_Carve.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences.

    -c,  --cloningvector: A vector to serve as the cloning plasmid.

    -d,  --destinationvector : A vector to serve as the delivery plasmid.

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