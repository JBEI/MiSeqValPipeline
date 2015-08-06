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
my $DNAV = 'DNAssemble_100_CreateXML_' . $VERSION;
my $DNAS = '_100';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
  'input=s'     => \$p{INPUT},
  'name=s'      => \$p{NAME},
  'output:s'    => \$p{OUTPUT},
  'logfile:s'   => \$p{LOGPATH},
  'help'        => \$p{HELP}
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
if (! $p{INPUT} || ! -e $p{INPUT})
{
  print {$OUTH} "\nDNA_ERROR: No input file provided or extant.\n";
  print {$OUTH} "\n\n******* $DNAV BAILED\n\n";
  close $OUTH;
  croak;
}
my ($iterator, $filename, $filesuffix) = $DNA->import_seqs($p{INPUT});
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.xml};


################################################################################
################################# CONFIGURING ##################################
################################################################################
my $design = $DNA->create_design(-name => $p{NAME});

my @seqobjs;
while ( my $obj = $iterator->next_seq() )
{
  push @seqobjs, $obj;
}

@seqobjs = sort {length $a->seq <=> length $b->seq} @seqobjs;

################################################################################
################################## CLEANING  ###################################
################################################################################
my %seqhsh;
my $count = 0;
Readonly my $PAD => 3;
foreach my $obj (@seqobjs)
{
  my $name = $obj->id;
  
  #Verify sequence legitimacy
  my ($seq, $lefts) = $DNA->purify($obj->seq);
  my @remain = @{$lefts};
  my $reason = undef;
  if (scalar @{$lefts})
  {
    print {$OUTH} "$name had bad characters (@remain) in its sequence\n";
    $reason = 'forbidden characters (';
    $reason .= join q{, }, @remain;
    $reason .= q{)};
  }
  
  #Check for duplicates
  $seq = uc $seq;
  if (! exists $seqhsh{$seq})
  {
    $seqhsh{$seq} = $name;
  }
  else
  {
    print {$OUTH} "$name is identical to $seqhsh{$seq} and will ";
    print {$OUTH} "be skipped\n";
    next;
  }

  #Create construct
  my $construct = $design->create_construct({
    id       => $name,
    kind     => 'unknown',
    sequence => $seq,
  });
  if ($reason)
  {
    $design->warn(-construct => $construct, -note => $reason);
  }
  $count++;
  
  #Pull subfeatures
  my @seqfeats = $obj->get_SeqFeatures;
  $construct->add_features(\@seqfeats);
}

print {$OUTH} "$count constructs added to xml file.\n";

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

__END__

=head1 NAME

  DNAssemble_100_CreateXML.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION


=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing DNA sequences. Genbank, FASTA, etc

    -n,  --name : a name for the root element of the design

  Optional arguments:

    -o,  --output : a filepath for dumping output

    -l,  --logfile : A logfile to update

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