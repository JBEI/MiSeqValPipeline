#!/usr/bin/env perl

use Bio::GeneDesign;
use Spreadsheet::WriteExcel;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_Create_Tracker_$VERSION";

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT},
      'output=s'          => \$p{OUTPUT},
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();
my $RES = $GD->set_restriction_enzymes(-enzyme_set => 'standard_and_IIB');

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
my $xlfile = $filename . "_TRACKING.xls";
my $xlpath = $p{OUTPUT} . $xlfile;
my $track = Spreadsheet::WriteExcel->new($xlpath);
$track->compatibility_mode();
my $bbcount = 0;
my %intermediates = ();
my @finalbbs = ();
my %objkeys = ();
my %objvecs = ();
my %enzlist = ();
my %enzunits = ();
my %vectors = ();
my %vstats = ();
my %vseqs = ();


my $header = $track->add_format();
$header->set_bold();
$header->set_bg_color('yellow');
   
################################################################################
################################## GATHERING ###################################
################################################################################
while ( my $obj = $iter->next_seq() )
{
  my $id = $obj->id;
  my $collec = $obj->annotation;
  my @comments = $collec->get_Annotations('comment');
  my $destvecname = undef;
  my $origname = undef;
  foreach my $comment (@comments)
  {
    my $ctext = $comment->as_text;
    if ($ctext =~ m{destination_vector \= (.+)}ms)
    {
      $destvecname = $1;
    }
    elsif ($ctext =~ m{original_name \= (.+)}ms)
    {
      $origname = $1;
    }
  }
  $objkeys{$id} = $origname;
  $objvecs{$id} = $destvecname;
  
  my @seqfeats = $obj->get_SeqFeatures;
  
  
  #Does this sequence have assembly oligos?
  my @ols = grep {$_->primary_tag eq 'assembly_oligo'} @seqfeats;
  my %ollu;
  unless (scalar(@ols))
  {
    print ("Fragment " . $id . " has no annotated assembly oligos... skipping!\n");
    next;
  }
  foreach my $ol (@ols)
  {
    my $olname = join q{}, $ol->get_tag_values("name");
    my @splitz = split m{\.}, $olname;
    my $bbname = $splitz[0] . q{.} . $splitz[1];
    if (exists $ollu{$bbname})
    {
      push @{$ollu{$bbname}}, $olname;
    }
    else
    {
      $ollu{$bbname} = [$olname];
    }
  }
  
  #Does this sequence have building blocks?
  my @bbs = grep {$_->primary_tag eq "building_block"} @seqfeats;
  unless (scalar(@bbs))
  {
    print ("Part " . $id . " has no annotated building blocks... skipping!\n");
    next;
  }
  foreach my $bb (@bbs)
  {
    $bbcount++;
    my $bbname = join("", $bb->get_tag_values("name"));
    my $bbseq = $bb->seq->seq;
    if ($bb->has_tag("BBprefix") && ! $p{SKIP})
    {
      $bbseq = join("", $bb->get_tag_values("BBprefix")) . $bbseq;
    }
    if ($bb->has_tag("BBsuffix") && ! $p{SKIP})
    {
      $bbseq = $bbseq . join("", $bb->get_tag_values("BBsuffix"));
    }
    $bbseq =~ s/\s//g;
    my $bbgc = join(q{}, $bb->get_tag_values('gcp'));
    my $bblen = length($bbseq);
    my $bbvec = $bb->has_tag("vector")  ? join("", $bb->get_tag_values("vector")) : q{};
    $vectors{$bbvec} = 1 if ($bbvec ne q{});
    my $intname = undef;
    my $skipflag = 0;
    if (! exists $ollu{$bbname})
    {
      print "$bbname has no annotated assembly oligos... skipping!\n";
      $skipflag++;
    }
    #push @finalbbs, [$bbcount, $bbname, $bbseq, $bblen, $bbgc, $bbvec, $intname, $skipflag];
    push @finallbbs, $bb;
  }
}

#Does this sequence have intermediates?
my @ints = grep {$_->primary_tag eq "intermediate"} @seqfeats;
unless (scalar(@ints))
{
  print ("Part " . $id . " has no intermediates!\n");
}
foreach my $int (@ints)
{
  if ($int->has_tag('enz_t') || $int->has_tag('enz_f'))
  {
    my $enz_f = join q{}, $int->get_tag_values('enz_f');
    my $enz_t = join q{}, $int->get_tag_values('enz_t');
    $enzlist{$enz_f} += 1;
    $enzlist{$enz_t} += 1;
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
}


foreach my $vector (keys %vectors)
{
  my $vecobj = $GD->load_vector(-name => $vector);
  $vectors{$vector} = $vecobj;
  my $newstart = $vecobj->chew3loc();
  my $tplen = length $vecobj->chew3;
  my $fplen = length $vecobj->chew5;
  my $oldseq = $vecobj->seq();
  my $newseq = substr $oldseq, $newstart - 2 + $tplen;
  $newseq .= substr $oldseq, 0, $newstart - 1 - $fplen;
  $vseqs{$vector} = $newseq;
  $vstats{$vector} = $GD->restriction_status(-sequence => $vseqs{$vector});
}

################################################################################
################################## REPORTING ###################################
################################################################################
my $currrow = 1;
my $currsheet = $track->add_worksheet("SynBio_set_track");
$currsheet->write(0, 0, "Sort", $header);
$currsheet->write(0, 1, "Gene Name", $header);
$currsheet->write(0, 2, "Seq", $header);
$currsheet->write(0, 3, "Size", $header);
$currsheet->write(0, 4, "GC", $header);
$currsheet->write(0, 5, "Cloning Vector", $header);
$currsheet->write(0, 6, "Assembly", $header);
$currsheet->write(0, 7, "Comments", $header);
$currsheet->write(0, 8, "correct_clone", $header);
$currsheet->write(0, 9, "mini_prep", $header);
$currsheet->write(0, 10, "Glycerol", $header);
$currsheet->write(0, 11, "Delivery_date", $header);
foreach my $bb (@finalbbs)
{
  next if ($bb->[7]);
  for (my $x = 0; $x <= 5; $x++)
  {
    $currsheet->write($currrow, $x, $bb->[$x]);
  }
  $currrow++;
}
$currrow = 1;

if (scalar keys %enzlist)
{
  $currsheet = $track->add_worksheet("Parts_to_Intermediates");
  $currsheet->write(0, 0, "Sort", $header);
  $currsheet->write(0, 1, "Part", $header);
  $currsheet->write(0, 2, 'Intermediate', $header);
  foreach my $bb (@finalbbs)
  {
    $currsheet->write($currrow, 0, $bb->[0]);
    $currsheet->write($currrow, 1, $bb->[1]);
    $currsheet->write($currrow, 2, $bb->[6]);
    $currrow++;
  }

  $currrow = 1;

  $currsheet = $track->add_worksheet("Intermediates");
  $currsheet->write(0, 0, "Intermediate", $header);
  $currsheet->write(0, 1, "First reaction", $header);
  $currsheet->write(0, 2, "Second reaction", $header);
  $currsheet->write(0, 3, "Seqlen", $header);
  $currsheet->write(0, 4, "Seq and veclen", $header);
  $currsheet->write(0, 5, "Seq", $header);
  foreach my $intname (sort keys %intermediates)
  {
    $currsheet->write($currrow, 0, $intname);
    my $seq = $intermediates{$intname}->{seq};
    my $vec = $intermediates{$intname}->{cvec};
    my $clonvec = $vectors{$vec};
    my $seqlen = length $seq;
    my $wholeseq = $seq . $vseqs{$vec};
    my $veclen = length $wholeseq;
  
    $currsheet->write($currrow, 3, $seqlen);
    $currsheet->write($currrow, 4, $veclen);
    $currsheet->write($currrow, 5, $seq);
    my $enzf = $intermediates{$intname}->{enzf};
    my $enzt = $intermediates{$intname}->{enzt};
    my $renzf = $RES->{$enzf};
    my $renzt = $RES->{$enzt};
    my $renzttemp = $renzt->temp();
    my $renzftemp = $renzf->temp();

    ##Check to see if cutting will obscure vector band
    my $bands = bands([$renzf, $renzt], $vseqs{$vec});
    my $vprob = find_obscuring_band($seqlen, $bands, 500);
  
    ##Determine units etc, optimal v suboptimal conditions
    my $bestbuffer = $renzf->acceptable_buffer($renzt, 50);
    my $firstrxn = q{};
    my $secondrxn = q{};
    if (defined $bestbuffer && $renzttemp == $renzftemp)
    {
      my $units = $renzf->units($bestbuffer, $wholeseq);
      $enzunits{$enzf} += $units;
      $firstrxn .= $enzf . ' x ' . $units . q{ };
      if ($enzf ne $enzt)
      {
        my $rtunits = $renzt->units($bestbuffer, $wholeseq);
        $enzunits{$enzt} += $rtunits;
        $firstrxn .= ' + ' . $enzt . ' x ' . $rtunits . q{ };
      }

      if (defined $vprob)
      {
        my $istats = $GD->restriction_status(-sequence => $seq);
        my $pstats = $GD->restriction_status(-sequence => $vprob);
        my @cands;
        my @ttiers;
        my @otiers;
        foreach my $enz (keys %{$RES})
        {
          next if ($istats->{$enz} > 0);
          next if ($pstats->{$enz} == 0);
          push @cands, $enz;
        }
        my @sort = sort {$RES->{$a}->score <=> $RES->{$b}->score} @cands;
        my $want = undef;
        foreach my $enz (@sort)
        {
          my $bands = bands([$RES->{$enz}], $vprob);
          my $sameflag = find_obscuring_band($seqlen, $bands, 500);
          next if defined $sameflag;
          my $etemp = $RES->{$enz}->temp();
          my $buff = $RES->{$enz}->buffers->{$bestbuffer};
          if ($etemp ne $renzttemp || ! $buff || $buff < 50)
          {
            push @otiers, $enz;
            next;
          }
          if (exists $enzlist{$enz})
          {
            $want = $enz;
            last;
          }
          push @ttiers, $enz;
        }
        $want = $want || $ttiers[0];
        if (defined $want)
        {
          my $units = $RES->{$want}->units($bestbuffer, $wholeseq);
          $enzunits{$want} += $units;
          $firstrxn .= ' + ' . $want . ' x ' . $units . q{ };
          $enzlist{$want}++;
        }
        elsif (scalar @otiers)
        {
          $want = $otiers[0];
          my $renz = $RES->{$want};
          my $usebuffer = $renz->acceptable_buffer($renz);
          my $units = $renz->units($usebuffer, $wholeseq);
          $enzunits{$want} += $units;
          $secondrxn .= $want . ' x ' . $units . q{ in };
          $secondrxn .= $usebuffer . ' @' . $renz->temp();
        }
        else
        {
          print "Intermediate $intname will have a vector band problem\n";
          $secondrxn .= ' Vector band will obscure';
        }
      }
      $firstrxn .= q{in } . $bestbuffer . ' @' . $renzftemp;
      $currsheet->write($currrow, 1, $firstrxn);
      $currsheet->write($currrow, 2, $secondrxn);
    }
    else
    {
      my $fbuffers = $renzf->common_buffers($renzf);
      my $fusebuff = $fbuffers->[0];
      my $rfunits = $renzf->units($fusebuff, $wholeseq);
      $enzunits{$enzf} += $rfunits;
      $firstrxn .= $enzf . ' x ' . $rfunits;
      my $firstconds = q{ in } . $fusebuff . q{ @ } . $renzftemp . q{};

      my $tbuffers = $renzt->common_buffers($renzt);
      my $tusebuff = $tbuffers->[0];
      my $rtunits = $renzt->units($tusebuff, $wholeseq);
      $enzunits{$enzt} += $rtunits;
      $secondrxn .= $enzt . ' x ' . $rtunits;
      my $secondconds = q{ in } . $tusebuff . q{ @ } . $renzttemp . q{};

      if (defined $vprob)
      {
        my $istats = $GD->restriction_status(-sequence => $seq);
        my $pstats = $GD->restriction_status(-sequence => $vprob);
        my @cands;
        my @ftiers;
        my @stiers;
        foreach my $enz (keys %{$RES})
        {
          next if ($istats->{$enz} > 0);
          next if ($pstats->{$enz} == 0);
          push @cands, $enz;
        }
        my @sort = sort {$RES->{$a}->score <=> $RES->{$b}->score} @cands;
        my $want = undef;
        my $rxnnum = 0;
        foreach my $enz (@sort)
        {
          my $renz = $RES->{$enz};
          my $bands = bands([$renz], $vprob);
          my $sameflag = find_obscuring_band($seqlen, $bands, 500);
          next if defined $sameflag;
          my $etemp = $renz->temp();
          my $ebuffs = $renz->buffers;
          next if ($etemp ne $renzftemp && $etemp ne $renzttemp);
          my $swit = $etemp eq $renzftemp ? 1 : 0;
          my $buff = $swit ? $ebuffs->{$fusebuff} : $ebuffs->{$tusebuff};
          next if (! $buff || $buff < 50);
          if (exists $enzlist{$enz})
          {
            if ($swit)
            {
              $want = $enz;
              $rxnnum = 1;
            }
            else
            {
              $want = $enz;
              $rxnnum = 2;
            }
            last;
          }
          elsif ($swit)
          {
            push @ftiers, $enz;
          }
          else
          {
            push @stiers, $enz;
          }
        }
        if ($rxnnum == 0)
        {
          if (defined $ftiers[0] && defined $stiers[0])
          {
            my $fcost = $RES->{$ftiers[0]}->score();
            my $scost = $RES->{$stiers[0]}->score();
            if ($fcost < $scost)
            {
              $rxnnum = 1;
              $want = $ftiers[0];
            }
            else
            {
              $rxnnum = 2;
              $want = $stiers[0];
            }
          }
          else
          {
            my $swit = scalar @ftiers;
            $want = $swit  ? $ftiers[0]  : $stiers[0];
            $rxnnum = $swit ? 1 : 2;
          }
        }
        if ($rxnnum == 1)
        {
          my $renz = $RES->{$want};
          my $units = $renz->units($fusebuff, $wholeseq);
          $enzunits{$want} += $units;
          $firstrxn .= ' + ' . $want . ' x ' . $units . q{ };
          $enzlist{$want}++;
        }
        elsif ($rxnnum == 2)
        {
          my $renz = $RES->{$want};
          my $units = $renz->units($tusebuff, $wholeseq);
          $enzunits{$want} += $units;
          $secondrxn .= ' + ' . $want . ' x ' . $units . q{ };
          $enzlist{$want}++;
        }
        else
        {
          print "Intermediate $intname will have a vector band problem\n";
        }
      }
      $currsheet->write($currrow, 1, $firstrxn . $firstconds);
      $currsheet->write($currrow, 2, $secondrxn . $secondconds);
    }
    $currrow++;
  }

  $currrow = 1;
  $currsheet = $track->add_worksheet("Enzymes");
  $currsheet->write(0, 0, "Enzyme", $header);
  $currsheet->write(0, 1, "Times Used", $header);
  $currsheet->write(0, 2, "Units Minimum", $header);
  $currsheet->write(0, 3, "Cost", $header);
  $currsheet->write(0, 4, "Recognizes", $header);
  foreach my $enz (sort {$enzlist{$b} <=> $enzlist{$a}} keys %enzlist)
  {
    $currsheet->write($currrow, 0, $enz);
    $currsheet->write($currrow, 1, $enzlist{$enz});
    $currsheet->write($currrow, 2, $enzunits{$enz});
    my $cost = sprintf("%.2f", $RES->{$enz}->score() * $enzunits{$enz});
    $currsheet->write($currrow, 3, $cost);
    $currsheet->write($currrow, 4, $RES->{$enz}->cutseq());
    $currrow++;
  }
}


$currrow = 1;

$currsheet = $track->add_worksheet("Usernames");
$currsheet->write(0, 0, "Sort", $header);
$currsheet->write(0, 1, "User's Name", $header);
$currsheet->write(0, 2, "Destination Vector", $header);
foreach my $id (sort keys %objkeys)
{
  $currsheet->write($currrow, 0, $id);
  $currsheet->write($currrow, 1, $objkeys{$id}) if (defined $objkeys{$id});
  $currsheet->write($currrow, 2, $objvecs{$id}) if (defined $objvecs{$id});
  $currrow++;
}

$currrow = 1;
$currsheet = $track->add_worksheet("PCA Primer cocktail");

$currsheet = $track->add_worksheet("Pacbio_references");
$currsheet->write(0, 0, "Sort", $header);
$currsheet->write(0, 1, "Gene Name", $header);
$currsheet->write(0, 2, "Copy_Paste_txt", $header);
foreach my $bb (@finalbbs)
{
  $currsheet->write($currrow, 0, $bb->[0]);
  $currsheet->write($currrow, 1, $bb->[1]);
  my $ref = ">" . $bb->[1] . "\n" . $bb->[2];
  $currsheet->write($currrow, 2, $ref);
  $currrow++;
}
$currrow = 1;

$currsheet = $track->add_worksheet("PacBio Data Raw");
$currsheet->write(0, 0, "Sort", $header);
$currsheet->write(0, 1, "Gene Name", $header);
foreach my $bb (@finalbbs)
{
  $currsheet->write($currrow, 0, $bb->[0]);
  $currsheet->write($currrow, 1, $bb->[1]);
  $currrow++;
}
$currrow = 1;

$currsheet = $track->add_worksheet("PacBio_results");
$currsheet->write(0, 0, "Sort", $header);
$currsheet->write(0, 1, "Gene Name", $header);
foreach my $bb (@finalbbs)
{
  $currsheet->write($currrow, 0, $bb->[0]);
  $currsheet->write($currrow, 1, $bb->[1]);
  $currrow++;
}
$currrow = 1;

$track->close();
print "\n\n";

print "Wrote $xlpath\n";

print "brought to you by $JGISBV\n\n";

exit;




__END__

=head1 NAME

  JGISB_Create_Tracker.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

  Converts a genbank file containing building blocks to an excel tracking sheet

=head1 USAGE

=head1 ARGUMENTS

Required arguments:

  -i,   --input : A genbank file to use as an input

Optional arguments:

  -o,   --output: Where should files be written
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