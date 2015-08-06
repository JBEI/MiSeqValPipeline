#!/usr/bin/env perl

use Bio::GeneDesign;
use Bio::SeqFeature::Generic;
use Bio::Seq;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

use strict;
use warnings;

my $VERSION = '1.00';
local $| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'input=s'         => \$p{INPUT},
      'output=s'        => \$p{OUTPUT},
      'format=s'        => \$p{FORMAT}
);
pod2usage() if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################

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

$p{FORMAT} = $p{FORMAT} || "genbank";
$p{EXT} = $p{FORMAT} eq "fasta" ? ".fasta" : ".gb";

################################################################################
################################### CARVING ####################################
################################################################################
print "gosh, beginning\n";

while ( my $obj = $iter->next_seq() )
{
  
  my $collec = $obj->annotation();
  my @keys = $collec->get_all_annotation_keys();
  my @comments = $collec->get_Annotations('comment');
  my $destvecname = undef;
  my $origname = undef;
  foreach my $comment (@comments)
  {
    my $ctext = $comment->as_text();
    if ($ctext =~ m{destination_vector \= (.+)}ms)
    {
      $destvecname = $1;
    }
    elsif ($ctext =~ m{original_name \= (.+)}ms)
    {
      $origname = $1;
    }
  }
  my $chunkname = $origname || $obj->id();
  if (! defined $destvecname)
  {
    print "Skipping $chunkname; no destination vector\n";
    next;
  }
  my $vector = $GD->load_vector(-name => $destvecname);
  if (! defined $vector->chew5 || ! defined $vector->chew3)
  {
    die "\n ERROR: $destvecname is not annotated for chewback ligation\n";
  }
  
  my $vobj = $vector->seqobj;
  my $orient = $vector->chew5loc > $vector->chew3loc  ? -1  : 1;
  my $fcoord = $orient == -1
      ? $vector->chew3loc + length $vector->chew3
      : $vector->chew5loc + length $vector->chew5;
  my $vseq = $vobj->seq;
  
  
  my $backbone = $vobj->clone();
  my $chseq = $obj->seq;
  my $fpseq = substr($vseq, 0, $fcoord);
  my $tpseq = substr($vseq, $fcoord - 1);
  $backbone->seq($fpseq . $chseq . $tpseq);
  $backbone->id($vector->name . "_" . $chunkname);
  my $Ffeat = Bio::SeqFeature::Generic->new
  (
    -primary => "insert",
    -start => $fcoord +1,
    -end => $fcoord + length($chseq),
    -tag => {label => $chunkname}
  );
  $backbone->add_SeqFeature($Ffeat) || die "no go on the feat add boss\n";
  
  
  my $outputfilename = $backbone->id . $p{EXT};
  my $outputpath = $p{OUTPUT} . "/" . $outputfilename;
  open (my $OUTFH, ">" . $outputpath );
  my $FOUT = Bio::SeqIO->new(-fh=>$OUTFH, -format=>$p{FORMAT}, -width => 80);
  if ($p{FORMAT} eq "genbank")
  {

  }
  $FOUT->write_seq($backbone);
  close $OUTFH;
}

exit;

__END__

=head1 NAME

  JGISB_ForceFeedVector.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

  My head hurts.

=head1 USAGE

  Take every insert in JBEI.FASTA and, in the current directory, create a new
  vector insert genbank file
    ./JGISB_ForceFeedVector.pl -i JBEI.FASTA -f genbank -o .

=head1 ARGUMENTS

  Required arguments:

    -i,  --input : a FASTA or Genbank file containing insert DNA sequences that
            are annotated to a vector.

  Optional arguments:

    -o, --output : a path in which to deposit building block sequences.
    -f, --format : what format should output be in?

=cut
