#!/usr/bin/env perl

use Bio::DNAssemble;
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use Readonly;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

our $VERSION = '1.00';
my $DNAV = 'DNAssemble_140_Cluster_' . $VERSION;
my $DNAS = '_140';

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'             => \$p{INPUT},
      'output:s'            => \$p{OUTPUT},
      'percid:f'            => \$p{PERCID},
      'logfile:s'           => \$p{LOGPATH},
      'help'                => \$p{HELP}
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
my $filename = $design->filename();
$filename = $1 if ($filename =~ m{(.+)\_}msix);
$p{OUTPUT} = $p{OUTPUT} || $filename . $DNAS . q{.xml};

Readonly my $PERCID =>  0.8;
$p{PERCID} = $p{PERCID} || $PERCID;
my %vectors = ();

my @CLEANUP;


################################################################################
################################# CLUSTERING  ##################################
################################################################################
# Cluster building blocks
#
my $bbpath = $filename . q{_BBS.fasta};
$design->constructs_to_fasta(
  -kind => 'building_block',
  -path => $bbpath,
  -sort => 'size'
);
my %changedbbs;
if (-e $bbpath)
{
  my $ucfilename = $filename .'_BBS.uc';
  my $cmd = "usearch -cluster_smallmem $bbpath -id $p{PERCID} -uc $ucfilename";
  my $output = $DNA->safeopen($cmd);
  my $bbclusters = ucparse($ucfilename);
  my @keys = keys %{$bbclusters};
  my $clusternumber = scalar @keys;
  my $maxclustercount = 0;
  foreach my $cent (@keys)
  {
    my $count = scalar @{$bbclusters->{$cent}};
    $maxclustercount = $count if ($count > $maxclustercount);
  }
  if ($maxclustercount >= 2)
  {
    my %barcodes = $DNA->make_barcodes($maxclustercount - 1);
    my @barlist = map {[$_, $barcodes{$_}]} keys %barcodes;
    my $clustercount = 0;
    foreach my $cent (@keys)
    {
      $clustercount++;
      my @list = sort @{$bbclusters->{$cent}};
      pop @list;
      my $matchcount = 0;
      foreach my $nodename (@list)
      {
        my @bbs = $design->get_constructs(-kind => 'building_block', -id => $nodename);
        my $bb = $bbs[0];

        my $name = $bb->id();
        my $seq = $bb->sequence();
        my $cvname = $bb->vector();
        if (! exists $vectors{$cvname})
        {
          my $clonvec = $DNA->load_vector(-name => $cvname);
          $vectors{$cvname} = $clonvec;
        }
        my $cvec = $vectors{$cvname};
        my $cvecseq = $cvec->chew5;
        my $start = length $cvecseq;

        my $nseq = $cvecseq . $barlist[$matchcount]->[1] . substr $seq, $start;
        $bb->sequence($nseq);
        $bb->barcode($barlist[$matchcount]->[0]);
        $bb->istart($bb->istart() + $start);
        $changedbbs{$nodename} = [$cvec, $barlist[$matchcount]];
        $matchcount++;
      }
    }
    print {$OUTH} "$clusternumber BB clusters (max size $maxclustercount)\n";
  }
  push @CLEANUP, $bbpath;
  push @CLEANUP, $ucfilename;
}

# Update intermediate sequences, if they themselves are not deliverables
# and the changed building block is in the first pool
#
foreach my $bbname (keys %changedbbs)
{
  my @poolnodes = $design->get_pools(-subid => $bbname);
  foreach my $poolnode (@poolnodes)
  {
    next if ($poolnode->getAttribute('number') != 1);
    my $subconstruct = $poolnode->parentNode;
    my $int = Bio::DNAssemble::Construct->new(-node => $subconstruct->parentNode);
    next if ($int->deliverable());
    my ($cvec, $barcode) = @{$changedbbs{$bbname}};
    my $intseq = $int->sequence();
    my $cvecseq = $cvec->chew5;
    my $start = length $cvecseq;
    my $nseq = $cvecseq . $barcode->[1] . substr $intseq, $start;
    $int->barcode($barcode->[0]);
    $int->sequence($nseq);
  }
}

# Cluster intermediates
#
my $intpath = $filename . q{_INTS.fasta};
$design->constructs_to_fasta(
  -kind => 'intermediate',
  -path => $intpath,
  -sort => 'size'
);
if (-e $intpath)
{
  my $ucfilename = $filename . '_INTS.uc';
  my $cmd = "usearch -cluster_smallmem $intpath -id $p{PERCID} -uc $ucfilename";
  my $output = $DNA->safeopen($cmd);
  my $intclusters = ucparse($ucfilename);
  my @keys = keys %{$intclusters};
  my $clusternumber = scalar @keys;
  my $maxclustercount = 0;
  foreach my $cent (@keys)
  {
    my $count = scalar @{$intclusters->{$cent}};
    $maxclustercount = $count if ($count > $maxclustercount);
  }
  if ($maxclustercount >= 2)
  {
    # Intermediates can't be barcoded, (unless we back calc to bbs, maybe)
    # We have to sort into screening groups.
    #
    my @screengroups = ();
    my $clustercount = 0;
    foreach my $cent (@keys)
    {
      $clustercount++;
      my @list = sort @{$intclusters->{$cent}};
      my $thiscount = scalar @list;
      my $matchcount = 0;
      foreach my $nodename (@list)
      {
        my @ints = $design->get_constructs(-kind => 'intermediate', id => $nodename);
        my $int = $ints[0];
        $int->screening_group(q{i.} . $matchcount);
        $screengroups[$matchcount]++;
        $matchcount++;
      }
    }
  }
  print {$OUTH} "$clusternumber INT clusters (max size $maxclustercount)\n";
  push @CLEANUP, $intpath;
  push @CLEANUP, $ucfilename;
}

# Cluster chunks
#
my $chpath = $filename . q{_CHUNKS.fasta};
$design->constructs_to_fasta(
  -kind => 'chunk',
  -path => $chpath,
  -sort => 'size'
);
if (-e $chpath)
{
  my $ucfilename = $filename . '_CHUNKS.uc';
  my $cmd = "usearch -cluster_smallmem $chpath -id $p{PERCID} -uc $ucfilename";
  my $output = $DNA->safeopen($cmd);
  my $chclusters = ucparse($ucfilename);
  my @keys = keys %{$chclusters};
  my $clusternumber = scalar @keys;
  my $maxclustercount = 0;
  foreach my $cent (@keys)
  {
    my $count = scalar @{$chclusters->{$cent}};
    $maxclustercount = $count if ($count > $maxclustercount);
  }
  if ($maxclustercount >= 2)
  {
    # Chunks can't be barcoded, (unless we back calc to bbs, maybe)
    # We have to sort into screening groups.
    #
    my @screengroups = ();
    my $clustercount = 0;
    foreach my $cent (@keys)
    {
      $clustercount++;
      my @list = sort @{$chclusters->{$cent}};
      my $thiscount = scalar @list;
      my $matchcount = 0;
      foreach my $nodename (@list)
      {
        my @chs = $design->get_constructs(-kind => 'chunk', id => $nodename);
        my $ch = $chs[0];
        $ch->screening_group(q{c.} . $matchcount);
        $screengroups[$matchcount]++;
        $matchcount++;
      }
    }
  }
  print {$OUTH} "$clusternumber CHUNK clusters (max size $maxclustercount)\n";
  push @CLEANUP, $chpath;
  push @CLEANUP, $ucfilename;
}

################################################################################
################################## REPORTING ###################################
################################################################################
$DNA->cleanup(\@CLEANUP);

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
sub ucparse
{
  my ($path) = @_;
  open my $UCFILE, '<', $path;
  my $ref = do { local $INPUT_RECORD_SEPARATOR = <$UCFILE> };
  close $UCFILE;
  my %clusters;
  my @lines = split qq{\n}, $ref;
  foreach my $line (@lines)
  {
    my @cols = split q{\t}, $line;
    if ($cols[0] eq q{S})
    {
      my $name = $cols[8];
      $clusters{$name} = [$name];
    }
    elsif ($cols[0] eq q{H})
    {
      my $name = $cols[9];
      my $hit = $cols[8];
      push @{$clusters{$name}}, $hit;
    }
  }
  return \%clusters;
}


__END__


=head1 NAME

  DNAssemble_140_Cluster.pl

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