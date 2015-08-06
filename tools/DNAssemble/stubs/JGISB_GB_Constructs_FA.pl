#!/usr/bin/env perl

use Bio::GeneDesign;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_GB_Constructs_FA_$VERSION";

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT},
      'output=s'          => \$p{OUTPUT},
      'skip_adds'         => \$p{SKIP}
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();

#The input file must exist and be a format we care to read.
die "\n JGISB_ERROR: You must supply an input file.\n"
  if (! $p{INPUT});

my ($iter, $filename, $suffix) = $GD->import_seqs($p{INPUT});

#The output path must exist, and we'll need it to end with a slash
$p{OUTPUT} = $p{OUTPUT} || ".";
$p{OUTPUT} .= "/" if (substr($p{OUTPUT}, -1, 1) !~ /[\/]/);
die "\n GDERROR: $p{OUTPUT} does not exist.\n"
  if ($p{OUTPUT} && ! -e $p{OUTPUT});


################################################################################
################################# CONFIGURING ##################################
################################################################################
my $bpcount = 0;
my $count = 0;
my @finals = ();

################################################################################
################################## ARRAYING  ###################################
################################################################################
while ( my $obj = $iter->next_seq() )
{
  my $collec = $obj->annotation;
  my @keys = $collec->get_all_annotation_keys();
  my @comments = $collec->get_Annotations('comment');
  my $destvec = undef;
  my $origname = undef;
  foreach my $comment (@comments)
  {
    my $ctext = $comment->as_text;
    if ($ctext =~ m{destination_vector \= (.+)}ms)
    {
      $destvec = $GD->load_vector(-name => $1)
    }
    elsif ($ctext =~ m{original_name \= (.+)}ms)
    {
      $origname = $1;
    }
  }
  my $sequence = $obj->seq;
  $count++;
  $bpcount += length $sequence;
  if (! $p{SKIP} && $destvec)
  {
    my $pre = lc $destvec->chew5;
    $sequence = $pre . $sequence;
    $sequence .= lc $destvec->chew3;
  }
  my $id = $obj->id;
  $id .= q{ (} . $origname . q{)} if ($origname);
  my $construct = Bio::Seq->new(-id => $id, -seq => $sequence);
  push @finals, $construct;
}

my $outputfilename = $filename . "_constructs.fasta";
$GD->export_seqs(
        -filename   => $outputfilename,
        -path       => $p{OUTPUT},
        -format     => "fasta",
        -sequences  => \@finals
);

print "\nPrinted $count building blocks and $bpcount bases.\n\n";

print "brought to you by $JGISBV\n\n";

exit;

__END__

=head1 NAME

  JGISB_GB_Constructs_FA.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

  Converts a genbank file containing constructs to a fasta file

=head1 USAGE

=head1 ARGUMENTS

Required arguments:

  -i,   --input : A genbank file to use as an input

Optional arguments:

  -o,   --output: Where should files be written
  -s,   --skip_adds: Do not add suffixes or prefixes to construct sequences
  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012, JGI Syn Bio developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of the Joint Genome Institute, the Lawrence Berkeley National
Laboratory, the Department of Energy, and the JGI developers may not be used to
endorse or promote products derived from this software without specific prior
written permission.

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