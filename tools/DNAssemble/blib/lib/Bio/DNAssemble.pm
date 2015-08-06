#
# DNAssemble engine
#

=head1 NAME

Bio::DNAssemble

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <smrichardson@lbl.gov>

=cut

package Bio::DNAssemble;
use base qw(Bio::Root::Root);

use Bio::DNAssemble::ConfigData;
use Bio::DNAssemble::Design;
use Bio::DNAssemble::Oligo qw(:GD);
use Bio::DNAssemble::Vector;
use Bio::GeneDesign;
use File::Basename;
use Bio::SeqIO;
use IPC::Open2;
use English qw( -no_match_vars );
use Carp;
use autodie qw(open close unlink);

use strict;
use warnings;

my $VERSION = 1.00;

=head1 SOURCES

=cut



=head1 CONSTRUCTORS

=head2 new

Returns an initialized Bio::DNAssemble object.

This function reads the ConfigData written at installation, imports the
relevant sublibraries, and sets the relevant paths.

    my $DNA = Bio::DNAssemble->new();

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  bless $self, $class;

  $self->{script_path} = Bio::DNAssemble::ConfigData->config('script_path');
  $self->{conf} = Bio::DNAssemble::ConfigData->config('conf_path');
  $self->{conf} .= q{/} unless substr($self->{conf}, -1, 1) eq q{/};

  $self->{tmp_path} = Bio::DNAssemble::ConfigData->config('tmp_path');
  $self->{tmp_path} .= q{/} unless substr($self->{tmp_path}, -1, 1) eq q{/};

  $self->{EMBOSS} = Bio::DNAssemble::ConfigData->config('EMBOSS_support');
  if ($self->{EMBOSS})
  {
    require Bio::DNAssemble::Palindrome;
    import Bio::DNAssemble::Palindrome qw(:GD);
  }

  $self->{BLAST} = Bio::DNAssemble::ConfigData->config('BLAST_support');
  if ($self->BLAST)
  {
    require Bio::DNAssemble::Blast;
    import Bio::DNAssemble::Blast qw(:GD);
  }

  $self->{vmatch} = Bio::DNAssemble::ConfigData->config('vmatch_support');
  if ($self->{vmatch})
  {
    require Bio::DNAssemble::Vmatch;
    import Bio::DNAssemble::Vmatch qw(:GD);
  }

  $self->{version} = $VERSION;
  $self->{GD} = Bio::GeneDesign->new();
  return $self;
}

=head1 ACCESSORS

=cut

=head2 EMBOSS

returns a value if EMBOSS_support was vetted and approved during installation.

=cut

sub EMBOSS
{
  my ($self) = @_;
  return $self->{'EMBOSS'};
}

=head2 BLAST

returns a value if BLAST_support was vetted and approved during installation.

=cut

sub BLAST
{
  my ($self) = @_;
  return $self->{'BLAST'};
}

=head2 vmatch

returns a value if vmatch_support was vetted and approved during installation.

=cut

sub vmatch
{
  my ($self) = @_;
  return $self->{'vmatch'};
}

=head2 GD

return the GD object

=cut

sub GD
{
  my ($self) = @_;
  return $self->{'GD'};
}

=head1 FUNCTIONS

=cut

=head2 safeopen

run a command and wait for it to finish. return the output handle as a string.

=cut

sub safeopen
{
  my ($self, $cmd) = @_;
  my ($inh, $outh) = (undef, undef);
  my $pid = open2($outh, $inh, $cmd) || croak "oops on $cmd: $OS_ERROR";
  waitpid $pid, 0;
  my $parse = <$outh>;
  return $parse;
}

=head2 filter_uniqueness

NO TEST

=cut

sub filter_uniqueness
{
  my ($self, @args) = @_;

  my ($seqs) = $self->_rearrange([qw(sequences)], @args);

  $self->throw("no argument provided to filter_uniqueness")
	  unless $seqs;

  my $seqarr = _filter_uniqueness( $seqs );
  return $seqarr;
}

=head2 filter_palindromes

=cut

sub filter_palindromes
{
  my ($self, @args) = @_;

  $self->throw("EMBOSS palindrome is not available")
    unless $self->EMBOSS;

  my ($seqobjs, $minpallen, $maxpallen, $gaplimit, $mismatches)
    = $self->_rearrange([qw(
        sequences
        minpallen
        maxpallen
        gaplimit
        mismatches)], @args
  );

  $self->throw("no nucleotide sequences provided")
	  unless $seqobjs;

  $self->throw("sequences argument is not an array reference")
    unless ref($seqobjs) eq "ARRAY";

  $minpallen = $minpallen || 8;
  $maxpallen = $maxpallen || 20;
  my $maxlen = undef;
  foreach my $seqobj (@$seqobjs)
  {
    $self->throw(ref($seqobj) . " is not a Bio::Seq object")
      unless $seqobj->isa("Bio::Seq");

    $self->throw("$seqobj is not a nucleotide sequence")
      unless $seqobj->alphabet eq "dna";
    my $len = $seqobj->length();
    $maxlen = $len if (! $maxlen || $len > $maxlen);
  }
#  $maxpallen = sprintf "%.0f", .5 * $maxlen if ($maxpallen > .5 * $maxlen);

  my $writedir = $self->{tmp_path} || "./";
  $gaplimit = $gaplimit  || 10;
  $mismatches = $mismatches || 1;

  my $seqarr = _filter_palindromes( $seqobjs,
                                  $minpallen,
                                  $maxpallen,
                                  $gaplimit,
                                  $mismatches,
                                  $writedir);
  return $seqarr;
}

=head2 filter_vmatch

=cut

sub filter_vmatch
{
  my ($self, @args) = @_;

  $self->throw("vmatch filtering is not available")
    unless $self->vmatch;

  my ($parent, $seqobjs, $mismatches, $revcom)
    = $self->_rearrange([qw(
        parent
        sequences
        mismatches
        revcom)], @args
  );

  $self->throw("no nucleotide sequences provided")
	  unless $seqobjs;

  $self->throw("sequences argument is not an array reference")
    unless ref($seqobjs) eq "ARRAY";

  foreach my $seqobj (@$seqobjs)
  {
    $self->throw("sequence argument is not a Bio::Seq object")
      unless $seqobj->isa("Bio::Seq");

    $self->throw("$seqobj is not a nucleotide sequence")
      unless $seqobj->alphabet eq "dna";
  }

  $self->throw("no parent sequence provided")
	  unless $parent;

  $self->throw("parent argument is not a Bio::Seq object")
    unless $parent->isa("Bio::Seq");

  $self->throw("$parent is not a nucleotide sequence")
    unless $parent->alphabet eq "dna";

  my $writedir = $self->{tmp_path} || "./";
  $mismatches = $mismatches || 7;
  $revcom = $revcom || 1;

  my $seqarr = _filter_vmatch(  $parent,
                                $seqobjs,
                                $mismatches,
                                $revcom,
                                $writedir);
  return $seqarr;
}

=head2 search_vmatch

=cut

sub search_vmatch
{
  my ($self, @args) = @_;

  $self->throw("vmatch searching is not available")
    unless $self->vmatch;

  my ($parent, $seqobj, $mismatches, $revcom)
    = $self->_rearrange([qw(
        parent
        sequence
        mismatches
        revcom)], @args
  );

  $self->throw("no nucleotide sequence provided")
	  unless $seqobj;

  $self->throw("sequence argument is not a Bio::Seq object")
    unless $seqobj->isa("Bio::Seq");

  $self->throw("$seqobj is not a nucleotide sequence")
    unless $seqobj->alphabet eq "dna";

  $self->throw("no parent sequence provided")
	  unless $parent;

  $self->throw("parent argument is not a Bio::Seq object")
    unless $parent->isa("Bio::Seq");

  $self->throw("$parent is not a nucleotide sequence")
    unless $parent->alphabet eq "dna";

  my $writedir = $self->{tmp_path} || "./";
  $mismatches = $mismatches || 1;
  $revcom = $revcom || 1;

  my $hits = _search_vmatch(  $parent,
                                $seqobj,
                                $mismatches,
                                $revcom,
                                $writedir);
  return $hits;
}

=head2 filter_blast

=cut

sub filter_blast
{
  my ($self, @args) = @_;

  $self->throw("BLAST+ filtering is not available")
    unless $self->BLAST;

  my ($parent, $seqobjs, $percent, $revcom)
    = $self->_rearrange([qw(
        parent
        sequences
        percent_identity
        revcom)], @args
  );

  $self->throw("no nucleotide sequences provided")
	  unless $seqobjs;

  $self->throw("sequences argument is not an array reference")
    unless ref($seqobjs) eq "ARRAY";

  foreach my $seqobj (@$seqobjs)
  {
    $self->throw(ref($seqobj) . " is not a Bio::Seq object")
      unless $seqobj->isa("Bio::Seq");

    $self->throw("$seqobj is not a nucleotide sequence")
      unless $seqobj->alphabet eq "dna";
  }

  $self->throw("no parent sequence provided")
	  unless $parent;

  $self->throw(ref($parent) . " is not a Bio::Seq object")
    unless $parent->isa("Bio::Seq");

  $self->throw("$parent is not a nucleotide sequence")
    unless $parent->alphabet eq "dna";

  my $writedir = $self->{tmp_path} || "./";
  $percent = $percent || 75;
  $revcom = $revcom || 1;

  my $seqarr = _filter_blast( $parent,
                              $seqobjs,
                              $percent,
                              $writedir);
  return $seqarr;
}

=head2 purify

=cut

sub purify
{
  my ($self, @args) = @_;

  my ($sequence) = $self->_rearrange([qw(sequence)], @args);

  my $origseq = $sequence;
  $sequence =~ s{A}{}gmsix;
  $sequence =~ s{T}{}gmsix;
  $sequence =~ s{C}{}gmsix;
  $sequence =~ s{G}{}gmsix;

  my %disallowed = ();
  my $seqlen = length $sequence;
  foreach my $x (0 .. $seqlen - 1)
  {
    $disallowed{substr $sequence, $x, 1}++;
  }
  my @badchars = sort keys %disallowed;
  foreach my $bad (@badchars)
  {
    $origseq =~ s{$bad}{}gmsi;
  }
  return ($origseq, \@badchars);
}

=head2 make_barcodes

=cut

sub make_barcodes
{
  my ($self, $num) = @_;
  my %barcodes = ();
  while (scalar keys %barcodes < $num)
  {
    my $lastchar = undef;
    my $seq = q{};
    my $key = q{};
    foreach my $y (1 .. 6)
    {
      my $char = $self->{GD}->random_dna(-length => 1, -gc_percentage => 50);
      while ($lastchar && $char eq $lastchar)
      {
        $char = $self->{GD}->random_dna(-length => 1, -gc_percentage => 50);
      }
      $seq .= $char x 3;
      $key .= $char;
      $lastchar = $char;
    }
    $barcodes{$key} = $seq;
  }
  return %barcodes;
}

=head2 choose_overlaps

=cut

sub choose_overlaps
{
  my ($self, @args) = @_;
  my ($seqobj, $laplist, $laplenmin, $spanlen, $spanlenmin, $spanlenmax)
    = $self->_rearrange([qw(
      construct_sequence
      overlap_list
      overlap_len_min
      span_len
      span_len_min
      span_len_max)], @args
  );
  my $return = {error => undef, warning => undef, overlaps => []};
  
  $self->throw("No construct sequence object defined for choose_overlaps")
    if (! defined $seqobj);
  
  #Filter overlaps for palindromes and mispriming
  #If vmatch or blast are not available, filter for uniqueness
  my $lapcount = scalar @{$laplist};
  if ($self->EMBOSS)
  {
    $laplist = $self->filter_palindromes(-sequences => $laplist);
  }
  if ($self->vmatch)
  {
    $laplist = $self->filter_vmatch(-sequences => $laplist, -parent => $seqobj);
  }
  elsif ($self->BLAST)
  {
    $laplist = $self->filter_blast(-sequences => $laplist, -parent => $seqobj);
  }
  else
  {
    $laplist = $self->filter_uniqueness(-sequences => $laplist);
  }
  $lapcount = scalar @{$laplist};

  #Flatten overlap annotations for speed, figure out the average length of the
  #overlaps that passed filtering
  my %LAPS;
  my $seqlen = length $seqobj->seq;
  my $mask = 'o' x $seqlen;
  foreach my $lap (@{$laplist})
  {
    my $tm     = join q{}, map {$_->value} $lap->get_Annotations('Tm');
    my $lefann = join q{}, map {$_->value} $lap->get_Annotations('left');
    my $rigann = join q{}, map {$_->value} $lap->get_Annotations('right');
    $LAPS{$lap->seq} = [$tm, $lap->seq, $lefann, $rigann];
    my $laplen = $rigann - $lefann + 1;
    substr $mask, $lefann - 1, $laplen, 'i' x $laplen;
  }

  #Here we should assess spread and decide if this is possible
  my $gapsize = $spanlenmax - ($laplenmin * 2);
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
    $return->{error} = 'No oligos: ' . join q{, }, @errors;
    return $return;
  }

  my @overlaps = values %LAPS;
  my $avg = 0;
  $avg += length $_->[1] foreach @overlaps;
  $avg = sprintf '%.0f', $avg / (scalar @overlaps);

  #Decide what size spanners to go for and how many.
  #Adjust the target size to avoid length outliers.
  my $tarlen = $spanlen || sprintf '%.0f', ($spanlenmin + $spanlenmax) / 2;
  my $tarnum = sprintf '%.0f', $seqlen / ($tarlen - $avg);
  my $tarspanlen = sprintf '%.0f', $seqlen / $tarnum;
  my $diff = $seqlen - (($tarnum * $tarspanlen) - ($avg * ($tarnum - 1)));
  #print "\ttarget: $tarnum spanners size $tarspanlen bp rem $diff\n";
  if (abs($diff) >= $tarnum)
  {
    my $rem = sprintf '%.0f', $diff / $tarnum;
    $tarspanlen = $tarspanlen + $rem;
    $diff = $diff - ($tarnum * $rem);
  }
  #print "\t final: $tarnum spanners size $tarspanlen bp rem $diff\n";

  if ($lapcount == 0 && $tarnum > 1)
  {
    $return->{error} = 'No suitable overlaps found';
    return $return;
  }
  
  #Pick overlaps properly spaced so as to obey length parameters.
  #Reroll if length parameters are violated by a choice.
  my $rbound = $seqlen - $tarspanlen;
  my ($lpos, $rpos) = (1, 1);
  my %fails;
  my $redoflag = 0;
  my $counter = 1;
  my $amideadyet = 1;
  my @chosenlaps;
  while ($lpos < $rbound && $amideadyet < 500)
  {
    $amideadyet++;
    my $rtarget = ($lpos + $tarspanlen);
    my @laps = grep {! exists $fails{$_->[1]}} @overlaps;
    @laps = grep {$_->[3] <= $lpos + $spanlenmax} @laps;
    @laps = grep {$_->[3] >= $lpos + $spanlenmin} @laps;
    @laps = sort {abs($a->[3] - $rtarget) <=> abs($b->[3] - $rtarget) || length $a->[1] <=> length $b->[1] } @laps;
    if (! scalar @laps)
    {
      #discarding last choice
      my $disgrace = pop @chosenlaps;
      $fails{$disgrace->[1]}++ if (defined $disgrace);
      #print "\t\tadjusting from $tarspanlen ";
      $tarspanlen = $counter % 2 == 0 ? $tarspanlen + $counter : $tarspanlen - $counter;
      #print "to $tarspanlen ";
      $counter++;
      if (scalar @chosenlaps)
      {
        $rpos = $chosenlaps[-1]->[3];
        $lpos = $chosenlaps[-1]->[2];
      }
      else
      {
        ($rpos, $lpos) = (1, 1);
        %fails = ();
      }
    }
    else
    {
      my $lap = $laps[0];
      #print "\tChoosing overlap $lap->[1] ($lap->[2]..$lap->[3])\n";
      push @chosenlaps, $lap;
      $rpos = $chosenlaps[-1]->[3];
      $lpos = $chosenlaps[-1]->[2];
      $redoflag = 0;
      if ($lpos >= $rbound)
      {
        if ((($seqlen - $lpos) + 1 < $spanlenmin) || (($seqlen - $lpos) + 1 > $spanlenmax))
        {
          #print "\t\t\tdiscarding last choice\n";
          my $disgrace = pop @chosenlaps;
          $fails{$disgrace->[1]}++;
          if (scalar @chosenlaps)
          {
            $rpos = $chosenlaps[-1]->[3];
            $lpos = $chosenlaps[-1]->[2];
          }
        }
      }
    }
    if ($amideadyet >= 500)
    {
      #print "DNAWARNING: Repeated search over 500 times!\n";
      $return->{error} = 'No oligos: no path, stuck';
      return $return;
    }
  }
  $return->{overlaps} = \@chosenlaps;
  return $return;
}

=head2 import_seqs

NO TEST

=cut

sub import_seqs
{
  my ($self, $path) = @_;
  $self->throw("$path does not exist.") if (! -e $path);
  my $iterator = Bio::SeqIO->new(-file => $path)
    || $self->throw("Can't parse $path!");

  my ($filename, $dirs, $suffix) = fileparse($path, qr/\.[^.]*/x);
  $suffix = substr($suffix, 1) if (substr($suffix, 0, 1) eq ".");
  $suffix = 'fasta' if ($suffix eq 'fa');

  return ($iterator, $filename, $suffix);
}

=head2 load_design_from_file

=cut

sub load_design_from_file
{
  my ($self, $path) = @_;
  $self->throw("No path provided.") if (! defined $path);
  $self->throw("$path does not exist.") if (! -e $path);
  my $design = Bio::DNAssemble::Design->new(-path => $path);
  return $design;
}

=head2 create_design

=cut

sub create_design
{
  my ($self, @args) = @_;
  my ($name) = $self->_rearrange([qw(name)], @args);
  my $design = Bio::DNAssemble::Design->new(-name => $name);
  return $design;
}

=head2 list_vectors

NO TEST

=cut

sub list_vectors
{
  my ($self) = @_;
  my $vpath = $self->{conf} . 'vectors/';
  my $vecs = {};
  opendir (my $VECDIR, $vpath) || croak "can't opendir $vpath";
  foreach my $file (readdir($VECDIR))
  {
    my $name = basename($file);
    $name =~ s{\.[^.]+$}{}x;
    if ($file =~ /\.genbank\Z/x)
    {
      $vecs->{$name} = $vpath . $file;
    }
  }
  closedir($VECDIR);
  return $vecs;
}

=head2 load_vector

NO TEST

=cut

sub load_vector
{
  my ($self, @args) = @_;

  my ($vname, $vpath) = $self->_rearrange([qw(name path)], @args);

  $vname = $vname || "custom_vector";

  if ($vname && ! $vpath)
  {
    $vpath = $self->{conf} . "vectors/" . $vname . ".genbank";
    $self->throw("Can't find $vname in configuration $vpath")
      unless (-e $vpath);
  }
  my $vector = Bio::DNAssemble::Vector->new(-name => $vname, -path => $vpath);

  return $vector;
}

=head2 _stripdown

=cut

sub _stripdown
{
  my ($self, $seqarg, $type, $enz_allowed) = @_;

  $enz_allowed = $enz_allowed || 0;
  my @seqs = ref($seqarg) eq "ARRAY" ? @$seqarg  : ($seqarg);
  my @list;
  foreach my $seq (@seqs)
  {
    my $str = $seq;
    my $ref = ref($seq);
    if ($ref)
    {
      my $bit = $self->_checkref($seq, $enz_allowed);
      $self->throw("object $ref is not a compatible object $bit") if ($bit < 1);
      $str = ref($seq->seq)  ? $seq->seq->seq  : $seq->seq;
    }
    push @list, uc $str;
  }
  return \@list if ($type eq "ARRAY");
  return $list[0];
}

=head2 _checkref

=cut

sub _checkref
{
  my ($self, $pobj, $enz_allowed) = @_;
  my $ref = ref $pobj;
  return -1 if (! $ref);
  $enz_allowed = $enz_allowed || 0;
  my ($bioseq, $bioseqfeat) = (0, 0);
  $bioseq = $pobj->isa("Bio::Seq");
  $bioseqfeat = $pobj->isa("Bio::SeqFeatureI");
  if ($enz_allowed)
  {
    $enz_allowed = $pobj->isa("Bio::DNAssemble::RestrictionEnzyme");
  }
  return $bioseq + $bioseqfeat + $enz_allowed;
}

=head1 PLEASANTRIES

=head2 pad

    my $name = 5;
    my $nice = $DNA->pad($name, 3);
    $nice == "005" || die;

    $name = "oligo";
    $nice = $DNA->pad($name, 7, "_");
    $nice == "__oligo" || die;

Pads an integer with leading zeroes (by default) or any provided set of
characters. This is useful both to make reports pretty and to standardize the
length of designations.

=cut

sub pad
{
  my ($self, $num, $thickness, $chars) = @_;
  my $t = $num;
  $chars = $chars || "0";
  $t = $chars . $t while (length($t) < $thickness);
  return $t;
}

=head2 attitude

    my $adverb = $DNA->attitude();

Ask DNAssemble how it handled your request.

=cut

sub attitude
{
  my @adverbs = qw(Elegantly Energetically Enthusiastically Excitedly Daintily
    Deliberately Diligently Dreamily Courageously Cooly Cleverly Cheerfully
    Carefully Calmly Briskly Blindly Bashfully Absentmindedly Awkwardly
    Faithfully Ferociously Fervently Fiercely Fondly Gently Gleefully Gratefully
    Gracefully Happily Helpfully Heroically Honestly Joyfully Jubilantly
    Jovially Keenly Kindly Knowingly Kookily Loftily Lovingly Loyally
    Majestically Mechanically Merrily Mostly Neatly Nicely Obediently Officially
    Optimistically Patiently Perfectly Playfully Positively Powerfully
    Punctually Properly Promptly Quaintly Quickly Quirkily Rapidly Readily
    Reassuringly Righteously Sedately Seriously Sharply Shyly Silently Smoothly
    Solemnly Speedily Strictly Successfully Suddenly Sweetly Swiftly Tenderly
    Thankfully Throroughly Thoughtfully Triumphantly Ultimately Unabashedly
    Utterly Upliftingly Urgently Usefully Valiantly Victoriously Vivaciously
    Warmly Wholly Wisely Wonderfully Yawningly Zealously Zestfully
  );
  my $num = scalar @adverbs;
  my $index = (sprintf "%.0f", rand($num)) % $num;
  return $adverbs[$index];
}

=head2 endslash

=cut

sub endslash
{
  my ($self, $path) = @_;
  if (-e $path && ! -d $path)
  {
    return $path;
  }
  if (substr($path, -1, 1) ne q{/})
  {
    $path .= q{/};
  }
  return $path;
}

=head2 cleanup

=cut

sub cleanup
{
  my ($self, $patharr) = @_;
  my @paths = grep {-e $_} @{$patharr};
  foreach my $path (@paths)
  {
    unlink $path;
  }
  return 1;
}

1;

__END__

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013, DNAssemble developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of Johns Hopkins, the Joint Genome Institute, the Lawrence Berkeley
National Laboratory, the Department of Energy, and the DNAssemble developers may
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