#!/usr/bin/env perl

use Bio::GeneDesign;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '5.00';
my $GDV = "GD_Design_Oligos_JGI_$VERSION";
my $GDS = "_OL";

local $| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'       => \$p{INPUT},
      'output=s'      => \$p{OUTPUT},
      'olilen:i'      => \$p{OLILEN},
      'laptemp:i'     => \$p{LAPTEMP},
      'poolsize:i'    => \$p{POOLSIZE},
      'maxpoolnum:i'  => \$p{MAXPOOLNUM},
      'maxolilen:i'   => \$p{MAXOLILEN},
      'minolilen:i'   => \$p{MINOLILEN},
      'minlaplen:i'   => \$p{MINLAPLEN},
      'alternate'     => \$p{ALTERNATE},
      'tolerance:f'   => \$p{TMTOL},
      'verbose'       => \$p{VERBOSE},
      'format=s'      => \$p{FORMAT},
      'help'          => \$p{HELP}
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();

#The input file must exist and be a format we care to read.
die "\n GDERROR: You must supply an input file.\n"
  if (! $p{INPUT});
my ($iterator, $filename, $suffix) = $GD->import_seqs($p{INPUT});

$p{FORMAT} = $p{FORMAT} || $suffix || "genbank";

#The output path must exist, and we'll need it to end with a slash
$p{OUTPUT} = $p{OUTPUT} || ".";
$p{OUTPUT} .= "/" if (substr($p{OUTPUT}, -1, 1) !~ /[\/]/);
die "\n GDERROR: $p{OUTPUT} does not exist.\n"
  if ($p{OUTPUT} && ! -e $p{OUTPUT});

################################################################################
################################# CONFIGURING ##################################
################################################################################
my @seqstowrite;
my @seqstoignore;
my @seqstofail;

$p{TMTOL} = $p{TMTOL}  ||  .5;
my $alg = $p{ALTERNATE}  ? 'JHU' : 'JGI';

my $totlen = 0;
my $totnum = 0;

################################################################################
################################## CHOPPING  ###################################
################################################################################
while ( my $obj = $iterator->next_seq() )
{
  my $chunkname = $obj->id;
  my $chseq = uc $obj->seq;
  my $chlen = length($chseq);
  
  print "Working on $chunkname...\n";
  
  my @seqfeats = $obj->get_SeqFeatures;
  
  #Does this sequence have building blocks?
  my @bbs = grep {$_->primary_tag eq "building_block"} @seqfeats;
  unless (scalar(@bbs))
  {
    print "\t$chunkname has no annotated building blocks... skipping!\n";
    push @seqstoignore, $obj;
    push @seqstowrite, $obj;
    next;
  }
  
  #Does this sequence already have assembly oligos?
  my @cos = grep {$_->primary_tag eq "assembly_oligo"} @seqfeats;
  if (scalar(@cos))
  {
    print "\t$chunkname already has assembly oligos... skipping!\n";
    push @seqstoignore, $obj;
    push @seqstowrite, $obj;
    next;
  }
  
  #Chop!
  foreach my $bb (@bbs)
  {
    my $bbname = join q{}, $bb->get_tag_values('name');
    my $oligos = $GD->chop_oligos(
        -building_block   => $bb,
        -algorithm        => $alg,
        -oligo_len_min    => $p{MINOLILEN},
        -oligo_len_max    => $p{MAXOLILEN},
        -oligo_len        => $p{OLILEN},
        -overlap_tm       => $p{LAPTEMP},
        -overlap_len_min  => $p{MINLAPLEN},
        -tm_tolerance     => $p{TMTOL},
        -pool_size        => $p{POOLSIZE},
        -pool_num_max     => $p{MAXPOOLNUM},
        -verbose          => $p{VERBOSE}
    );
    my $num = scalar(@$oligos);
    unless ($num > 0)
    {
      print "\t$bbname cannot be made into oligos... skipping\n\n";
      push @seqstofail, $bbname;
      next;
    }
    foreach my $oligo (@$oligos)
    {
      $obj->add_SeqFeature($oligo);
      $totlen += length(join(q{}, $oligo->get_tag_values("sequence")));
      $totnum ++;
    }
  }
  push @seqstowrite, $obj;
}

#Maybe report

print "\n\n";
my $outputfilename = $filename  . $GDS . "." . $p{FORMAT};
my $reportpath = $p{OUTPUT} . $filename . "_OLreport.txt";
open (my $REP, '>', $reportpath);
if (scalar @seqstoignore)
{
  print $REP $_->id . " : IGNORED FOR OLIGOS\n" foreach @seqstoignore;
}
if (scalar @seqstofail)
{
  print $REP $_ . " : FAILED OLIGOS\n" foreach @seqstofail;
}
close $REP;
if (scalar @seqstowrite)
{
  $GD->export_seqs(
          -filename   => $outputfilename,
          -path       => $p{OUTPUT},
          -format     => $p{FORMAT},
          -sequences  => \@seqstowrite);
}
print "\n\n";

print "\n";
print "Wrote " . $p{OUTPUT} . "$outputfilename\n";
print "Wrote $reportpath\n";
print "\n";
print $GD->attitude() . " brought to you by $GDV\n\n";

exit;

__END__

=head1 NAME

  GD_Design_Oligos.pl

=head1 VERSION

  Version 5.00

=head1 DESCRIPTION

    The Design_Oligos script will break each nucleotide sequence it is
    given into a set of overlapping assembly oligonucleotides. It uses EMBOSS
    palindrome to check for potential secondary structures.

    Output will be named by the identification of the part, and will be tagged
    with the _OL suffix. Default values are assumed for every parameter not
    provided.

    You must provide either a target oligo length or a minimum maximum range.
    You can provide both.

=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a file containing building block sequences.


  Optional arguments:
    --output : a path in which to deposit oligo sequences.
    --minlaplen : minimum length of overlap allowed
    --laptemp : Tm in temperature degrees C of oligo overlaps (default 70)
    --olilen : the length in bp of assembly oligos (default 180)
    --minolilen : minimum length of assembly oligos permitted (default 45)
    --maxolilen : maximum length of assembly oligos permitted (default 200)
    --tolerance : amount of +/- variation allowed in Tm (default 2.5)
    --poolsize : oligos will be pooled; GD will make bridging primers between
                  pools of the specified size
    --maxpoolnum : the maximum number of pools to create.
    --alternate : design oligos such that they alternate on opposite strands;
                the default is not facing, that is, all oligos except the very
                last are on the same strand
    --format : default genbank
    --verbose : show deliberations happening
    -h,  --help: display this message

=cut