#!/usr/bin/env perl

use Bio::DNAssemble;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_200_Check_' . $VERSION;
my $DNAS = '_200';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'             => \$p{INPUT},
      'output:s'            => \$p{OUTPUT},
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

my $poolsource = q{subconstructs/pool};
my $subsource = q{subconstruct};
my $allsource = q{subconstructs/pool/subconstruct};
my %vectors = ();

################################################################################
################################## CHECKING  ###################################
################################################################################

# Look through all the building blocks and propagate constructibility
#
my @bbs = $design->get_constructs(-kind => 'building_block');
my %magicbbs = ();
foreach my $bb (@bbs)
{
  my $method = $bb->method;

  my $bbname = $bb->id();
  my @sublist = $bb->xmlnode->findnodes($allsource);
  my @conlist = ();
  if (! scalar @sublist)
  {
    $method = 'magic';
  }
  foreach my $sub (@sublist)
  {
    my $subid = $sub->textContent;
    @conlist = $design->get_constructs(-kind => 'assembly_oligo', -id => $subid);
    if (! scalar @conlist)
    {
      $method = 'magic';
    }
  }
  my $seq = $bb->sequence();
  my $vecname = $bb->vector();
  if (! exists $vectors{$vecname})
  {
    $vectors{$vecname} = $DNA->load_vector(-name => $vecname);
  }
  my $vec = $vectors{$vecname};
  my $fpseq = substr($seq, 0, 150);
  my $tpseq = substr($seq, -150, 150);
  if ($fpseq !~ $vec->chew5)
  {
    print {$OUTH} $vec->chew5, "\n";
    print {$OUTH} "$bbname: vector sequence is missing 5 prime!\n";
    $design->weep(-construct => $bb, -note => 'missing 5 prime ' . $vecname . ' sequence' );
    $method = 'magic';
  }
  if ($tpseq !~ $vec->chew3)
  {
    print {$OUTH} $vec->chew3, "\n";
    print {$OUTH} "$bbname: vector sequence is missing 3 prime!\n";
    $design->weep(-construct => $bb, -note => 'missing 3 prime ' . $vecname . ' sequence' );
    $method = 'magic';
  }
  #$bb->method($method);
  if ($method eq 'magic')
  {
    print {$OUTH} "$bbname requires magic\n";
    $magicbbs{$bbname} = $method;
    next;
  }
  
}

# Look through all the intermediates and propagate constructibility
#
my @intermediates = $design->get_constructs(-kind => 'intermediate');
my %magicints = ();
foreach my $intermediate (@intermediates)
{
  my $method = $intermediate->method();
  my $intname = $intermediate->id();
  my @sublist = $intermediate->xmlnode->findnodes($allsource);
  my @bbs;
  if (! scalar @sublist)
  {
    $method = 'magic';
  }
  foreach my $sub (@sublist)
  {
    my $subid = $sub->textContent;
    if (exists $magicbbs{$subid})
    {
      $method = 'magic';
    }
    push @bbs, $design->get_constructs(-id => $subid, -kind => 'building_block');
  }
  my $seq = $intermediate->sequence();
  my $vecname = $intermediate->vector();
  if (! exists $vectors{$vecname})
  {
    $vectors{$vecname} = $DNA->load_vector(-name => $vecname);
  }
  my $vec = $vectors{$vecname};
  my $fpseq = substr($seq, 0, 150);
  my $tpseq = substr($seq, -150, 150);
  if ($fpseq !~ $vec->chew5)
  {
    print {$OUTH} $vec->chew5, "\n";
    print {$OUTH} "$intname: vector sequence is missing 5 prime!\n";
    $design->weep(-construct => $intermediate, -note => 'missing 5 prime ' . $vecname . ' sequence' );
    $method = 'magic';
  }
  if ($tpseq !~ $vec->chew3)
  {
    print {$OUTH} $vec->chew3, "\n";
    print {$OUTH} "$intname: vector sequence is missing 3 prime!\n";
    $design->weep(-construct => $intermediate, -note => 'missing 3 prime ' . $vecname . ' sequence' );
    $method = 'magic';
  }
  #$intermediate->method($method);
  if ($method eq 'magic')
  {
    $magicints{$intname} = $method;
    print {$OUTH} "$intname requires magic\n";
    $design->weep(-construct => $intermediate, -note => 'missing subconstruct');
    next;
  }
}


# Look through all the chunks and propagate constructibility
#
my @chunks = $design->get_constructs(-kind => 'chunk');
my $enzsource = q{reaction/enzyme};
foreach my $chunk (@chunks)
{
  my $method = $chunk->method();
  my $chname = $chunk->id();
  my @poollist = $chunk->xmlnode->findnodes($poolsource);
  if (! scalar @poollist)
  {
    $method = 'magic';
  }
  foreach my $pool (@poollist)
  {
    my @subslist = $pool->findnodes($subsource);
    my $subid = $subslist[0]->textContent;
    if (exists $magicints{$subid})
    {
      $method = 'magic';
    }
    my @conlist = $design->get_constructs(-id => $subid, -kind => 'intermediate');
    my $sub = $conlist[0];
    my $seq = $sub->sequence();
    my @enzslist = $pool->findnodes($enzsource);
    my %enzes;
    foreach my $enznode (@enzslist)
    {
      next if ($enznode->hasAttribute('type') && $enznode->getAttribute('type') eq 'vprob');
      $enzes{$enznode->getAttribute('name')}++;
    }
    my $flag = 0;
    $flag = 2 if (scalar keys %enzes == 1);
    foreach my $enz (keys %enzes)
    {
      my $pos = $DNA->GD->positions(-sequence => $seq, -query => $RES->{$enz});
      my $seecount = scalar keys %{$pos};
      if ($seecount != $enzes{$enz} && $flag && $seecount != $flag)
      {
        print {$OUTH} "$subid: $enz is not appearing the right number of times:";
        print {$OUTH} "\t$seecount times instead of $enzes{$enz}\n";
        $design->weep(-construct => $sub, -note => 'weird enzyme pattern for ' . $enz);
        $method = 'magic';
      }
    }
  }
  my $seq = $chunk->sequence();
  my $vecname = $chunk->vector();
  if (! exists $vectors{$vecname})
  {
    $vectors{$vecname} = $DNA->load_vector(-name => $vecname);
  }
  my $vec = $vectors{$vecname};
  my $fpseq = substr($seq, 0, 150);
  my $tpseq = substr($seq, -150, 150);
  if ($fpseq !~ $vec->chew5)
  {
    print {$OUTH} $vec->chew5, "\n";
    print {$OUTH} "$chname: vector sequence is missing 5 prime!\n";
    $design->weep(-construct => $chunk, -note => 'missing 5 prime ' . $vecname . ' sequence' );
    $method = 'magic';
  }
  if ($tpseq !~ $vec->chew3)
  {
    print {$OUTH} $vec->chew3, "\n";
    print {$OUTH} "$chname: vector sequence is missing 3 prime!\n";
    $design->weep(-construct => $chunk, -note => 'missing 3 prime ' . $vecname . ' sequence' );
    $method = 'magic';
  }
  #$chunk->method($method);
  if ($method eq 'magic')
  {
    print {$OUTH} "$chname requires magic\n";
    $design->weep(-construct => $chunk, -note => 'missing subconstruct');
  }
}

################################################################################
################################## REPORTING ###################################
################################################################################
$design->dump_xml($p{OUTPUT}, 1);

print {$OUTH} "\n\n";
print {$OUTH} "Wrote $p{OUTPUT}\n\n";
print {$OUTH} $DNA->attitude() . " brought to you by $DNAV\n\n";

print {$OUTH} "\n\n******* $DNAV FINISHED\n\n";
close $OUTH;
exit;


################################################################################
################################# SUBROUTINES ##################################
################################################################################


__END__


=head1 NAME

  DNAssemble_200_Check.pl

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