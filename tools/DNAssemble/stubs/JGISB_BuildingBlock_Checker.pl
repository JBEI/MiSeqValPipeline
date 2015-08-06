#!/usr/bin/env perl

use Bio::GeneDesign;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_BuildingBlock_Checker_$VERSION";

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT}
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $GD = Bio::GeneDesign->new();
my $RES = $GD->set_restriction_enzymes(-enzyme_set => 'standard_and_IIB');

#The input file must exist and be a format we care to read.
die "\n JGISB_ERROR: You must supply an input file.\n" if (! $p{INPUT});
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
my %bbreport = ();
my %chreport = ();

################################################################################
################################## ARRAYING  ###################################
################################################################################
while ( my $obj = $iter->next_seq() )
{
  my $chname = $obj->id;
  my @seqfeats = $obj->get_SeqFeatures;
    
  my $collec = $obj->annotation;
  my @keys = $collec->get_all_annotation_keys();
  my @comments = $collec->get_Annotations('comment');
  my $destvec = undef;
  foreach my $comment (@comments)
  {
    my $ctext = $comment->as_text;
    if ($ctext =~ m{destination_vector \= (.+)}ms)
    {
      my $destvecname = $1;
      $destvec = $GD->load_vector(-name => $destvecname);
    }
  }
  #Does this sequence have building blocks?
  my @bbs = grep {$_->primary_tag eq "building_block"} @seqfeats;
  my $bbcount = scalar @bbs;
  unless ($bbcount > 0)
  {
    print ($obj->id . " has no annotated building blocks... skipping!\n");
    next;
  }

  my @ints = grep {$_->primary_tag eq "intermediate"} @seqfeats;
  my %realints = ();
  foreach my $int (@ints)
  {
    my $name = join q{}, $int->get_tag_values("name");
    $realints{$name} = $int;
  }
  my %ints;
  
  #
  # Check the sanity of each individual building block
  #
  my %subsets;
  my $subsetflag = 0;
  my $count = 0;
  my $clonvec = undef;
  foreach my $bb (@bbs)
  {
    $count++;
    my $bbname = join("", $bb->get_tag_values("name"));
    if ($bb->has_tag('intermediate'))
    {
      my $subset = join q{}, $bb->get_tag_values('intermediate');
      $subsetflag++;
      push @{$ints{$subset}}, $bb;
    }
    elsif ($subsetflag != 0)
    {
      my $problem = "Building block may be missing a intermediate attribute!";
      push @{$bbreport{$bbname}}, $problem;
    }
    my $bbseq = $bb->seq->seq;
    if ($bb->has_tag("prefix"))
    {
      $bbseq = join("", $bb->get_tag_values("prefix")) . $bbseq;
    }
    if ($bb->has_tag("suffix"))
    {
      $bbseq = $bbseq . join("", $bb->get_tag_values("suffix"));
    }
    $bbseq =~ s/\s//g;
    my $bblen = length($bbseq);
    $bpcount += $bblen;

    my $prefixoffset = 0;
    my $suffixinset = 0;
    
    if (defined $destvec)
    {
      if ($count == 1 && $bbseq !~ $destvec->chew5)
      {
        my $problem = "First building block does not have destvec 5' sequence!";
        push @{$bbreport{$bbname}}, $problem;
      }
      elsif ($count == $bbcount && $bbseq !~ $destvec->chew3)
      {
        my $problem = "Last building block does not have destvec 3' sequence!";
        push @{$bbreport{$bbname}}, $problem;
      }
      elsif ($count != 1 && $count != $bbcount)
      {
        if ($bbseq =~ $destvec->chew5 || $bbseq =~ $destvec->chew3)
        {
          my $problem = "Middle building block contains destvec sequence!";
          push @{$bbreport{$bbname}}, $problem;
        }
      }
    }
    
    if ($bb->has_tag('vector'))
    {
      my $cloningvector = join q{}, $bb->get_tag_values('vector');
      if (! defined $clonvec)
      {
        $clonvec = $GD->load_vector(-name => $cloningvector);
      }
      elsif ($cloningvector ne $clonvec->name)
      {
        my $problem = "Cloning vector is not the same as first bb in chunk!";
        push @{$bbreport{$bbname}}, $problem;
      }
      $prefixoffset += length($clonvec->chew5);
      $suffixinset += length($clonvec->chew3);
      if ($bbseq !~ $clonvec->chew5)
      {
        my $problem = "cloning vector 5' region is not present!";
        push @{$bbreport{$bbname}}, $problem;
      }
      if ($bbseq !~ $clonvec->chew3)
      {
        my $problem = "cloning vector 3' region is not present!";
        push @{$bbreport{$bbname}}, $problem;
      }
    }
  }
  
  my $subsetcount = scalar keys %ints;
  my $subsetseen = 0;
  #
  # Check the sanity of subsets if applicable
  #
  if ($subsetcount && $bbcount <= 5)
  {
    my $problem = "Subset design with fewer than six building blocks!";
    push @{$chreport{$chname}}, $problem;
  }
  if ($subsetcount != scalar keys %ints)
  {
    my $problem = "Subset missing";
    push @{$chreport{$chname}}, $problem;
  }
  foreach my $subset (sort keys %realints)
  {
    my $int = $realints{$subset};
    $subsetseen++;
    my @set = @{$ints{$subset}};
    my %enzes = ();
    my ($fenz, $tenz) = (undef, undef);
    if (scalar @set > 5)
    {
      my $problem = "$subset has more than five non vector members!";
      push @{$chreport{$chname}}, $problem;
    }
    my $subsetseq = $int->seq->seq;
    
    my ($prefix, $suffix) = (q{}, q{});
    if ($int->has_tag("prefix"))
    {
      $prefix = join q{}, $int->get_tag_values("prefix");
      $prefix =~ s{\s}{}g;
    }
    if ($int->has_tag("suffix"))
    {
      $suffix = join q{}, $int->get_tag_values("suffix");
      $suffix =~ s{\s}{}g;
    }
    $subsetseq = $prefix . $subsetseq . $suffix;

    if (! $int->has_tag('enz_f'))
    {
      my $problem = "No five prime enzyme on $subset!";
      push @{$chreport{$chname}}, $problem;
    }
    if (! $int->has_tag('enz_t'))
    {
      my $problem = "No three prime enzyme on $subset!";
      push @{$chreport{$chname}}, $problem;
    }
    $fenz = join q{}, $int->get_tag_values('enz_f');
    $tenz = join q{}, $int->get_tag_values('enz_t');
    $enzes{$fenz}++;
    $enzes{$tenz}++;
    my $fpos = $GD->positions(-sequence => $subsetseq, -query => $RES->{$fenz});
    my @foffsets = keys %{$fpos};
    if (scalar @foffsets != $enzes{$fenz})
    {
      my $problem = "$subset: $fenz is not appearing the right number of times!";
      push @{$chreport{$chname}}, $problem;
    }
    my $tpos = $GD->positions(-sequence => $subsetseq, -query => $RES->{$tenz});
    my @toffsets = keys %{$tpos};
    if (scalar @toffsets != $enzes{$tenz})
    {
      my $problem = "$subset: $tenz is not appearing the right number of times!";
      push @{$chreport{$chname}}, $problem;
    }
    
    my $fpseq = substr($subsetseq, 0, 150);
    my $tpseq = substr($subsetseq, -150, 150);
    if ($clonvec)
    {
      if ($fpseq !~ $clonvec->chew5)
      {
        my $problem = "$subset: cloning vector sequence is missing 5 prime!";
        push @{$chreport{$chname}}, $problem;
      }
      
      if ($tpseq !~ $clonvec->chew3)
      {
        my $problem = "$subset: cloning vector sequence is missing 3 prime!";
        push @{$chreport{$chname}}, $problem;
      }
    }
    if ($destvec)
    {
      if ($subsetseen == 1 && $fpseq !~ $destvec->chew5)
      {
        my $problem = "$subset: dest vector sequence is missing 5 prime!";
        push @{$chreport{$chname}}, $problem;
      }
      elsif ($subsetseen == $subsetcount && $tpseq !~ $destvec->chew3)
      {
        my $problem = "$subset: dest vector sequence is missing 3 prime!";
        push @{$chreport{$chname}}, $problem;
      }
    }
  }
}

foreach my $bbname (sort keys %bbreport)
{
  my $report = $bbname . qq{\n\t};
  $report .= join qq{\n\t}, @{$bbreport{$bbname}};
  print "$report\n";
}

foreach my $chname (sort keys %chreport)
{
  my $report = $chname . qq{\n\t};
  $report .= join qq{\n\t}, @{$chreport{$chname}};
  print "$report\n";
}

if (! scalar keys %bbreport && ! scalar keys %chreport)
{
  print "\n\nNo problems found.\n";
}

print "brought to you by $JGISBV\n\n";

exit;

__END__

=head1 NAME

  JGISB_BuildingBlock_Checker.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

  Checks the sanity of Building Blocks

=head1 USAGE

=head1 ARGUMENTS

Required arguments:

  -i,   --input : A genbank file to use as an input

Optional arguments:

  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013, JGI Syn Bio developers
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