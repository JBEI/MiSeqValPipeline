#
# DNAssemble module for design objects
#

=head1 NAME

Bio::DNAssemble::Design

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

DNAssemble object that represents a design

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::DNAssemble::Design;

use File::Basename;
use XML::LibXML;
use Bio::DNAssemble::Construct;
use Bio::DNAssemble::Error;
use Bio::DNAssemble::Issue;
use Bio::DNAssemble::Warning;
use Bio::SeqIO;
use Bio::Seq;
use Carp;
use autodie qw(open close);

use strict;
use warnings;

use base qw(Bio::Root::Root);

our $VERSION = 1.00;

my @cks = qw(unknown building_block chunk assembly_oligo stitching_oligo intermediate);
my %CONKINDS = map {$_ => 1} @cks;

=head1 CONSTRUCTOR METHODS

=head2 new

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($name, $path) = $self->_rearrange([qw(NAME PATH)], @args);

  if (defined $path)
  {
    $self->throw("$path does not exist.") if (! -e $path);

    my ($filename, $filedirs, $filesuffix) = fileparse($path, qr/\.[^.]*/x);
    my $parser = XML::LibXML->new();
    $self->{doc} = $parser->load_xml(location => $path)
      || $self->throw("Can't parse $path!");
    $self->{root} = $self->{doc}->documentElement();
    my @cons = $self->{root}->getElementsByTagName('constructs');
    $self->{cons} = $cons[0];
    my @wars = $self->{root}->getElementsByTagName('warnings');
    $self->{warnings} = $wars[0];
    my @isss = $self->{root}->getElementsByTagName('issues');
    $self->{issues} = $isss[0];
    my @errs = $self->{root}->getElementsByTagName('errors');
    $self->{errors} = $errs[0];
    $self->{filepath} = $path;
    $self->{filename} = $filename;
    $self->{filesuffix} = $filesuffix;
  }
  else
  {
    $self->{doc} = XML::LibXML::Document->new('1.0', 'utf-8');
    $self->{root} = $self->{doc}->createElement('design');
    $self->{doc}->setDocumentElement($self->{root});
    $self->{cons} = $self->{doc}->createElement('constructs');
    $self->{root}->appendChild($self->{cons});
    $self->{root}->setAttribute('id', $name) if (defined $name);
    $self->{warnings} = $self->{doc}->createElement('warnings');
    $self->{root}->appendChild($self->{warnings});
    $self->{issues} = $self->{doc}->createElement('issues');
    $self->{root}->appendChild($self->{issues});
    $self->{errors} = $self->{doc}->createElement('errors');
    $self->{root}->appendChild($self->{errors});
  }

  return $self;
}

=head1 FUNCTIONAL METHODS

=head2 dump_xml

=cut

sub dump_xml
{
  my ($self, $path, $swit) = @_;;
  $swit = $swit || 0;
  open my $XML, '>', $path;
  print {$XML} $self->{doc}->toString($swit);
  close $XML;
  $self->{filepath} = $path;
  return 1;
}

=head1 CONSTRUCT METHODS

=head2 create_construct

=cut

sub create_construct
{
  my ($self, $atts) = @_;
  my $construct = Bio::DNAssemble::Construct->new(
    -doc => $self->{doc},
    -attributes => $atts
  );
  $self->{cons}->appendChild($construct->xmlnode());
  return $construct;
}

=head2 add_construct

=cut

sub add_construct
{
  my ($self, $construct) = @_;
  my $node = $construct->xmlnode();
  $self->{doc}->importNode($node);
  $self->{cons}->appendChild($node);
  return 1;
}

=head2 delete_construct

=cut

sub delete_construct
{
  my ($self, $construct) = @_;
  $self->{cons}->removeChild($construct->xmlnode());
  return;
}

=head2 get_constructs

=cut

sub get_constructs
{
  my ($self, @args) = @_;

  my ($kind, $id) = $self->_rearrange([qw(KIND ID)], @args);
  
  $self->throw("$kind is not a kind of construct\n")
    if ($kind && ! exists $CONKINDS{$kind});
    
  my $rulesource = '//construct';
  if ($kind || $id)
  {
    $rulesource .= q{[};
  }
  if ($id)
  {
    $rulesource .= q{@id = "} . $id . q{"};
  }
  if ($kind && $id)
  {
    $rulesource .= q{ and };
  }
  if ($kind)
  {
    $rulesource .= q{@kind = "} . $kind . q{"};
  }
  if ($kind || $id)
  {
    $rulesource .= q{]};
  }
  my @nodes = $self->{root}->findnodes($rulesource);
  my @constructs;
  foreach my $node (@nodes)
  {
    push @constructs, Bio::DNAssemble::Construct->new(-node => $node);
  }
  return @constructs;
}

=head2 constructs_to_fasta

=cut

sub constructs_to_fasta
{
  my ($self, @args) = @_;
  
  my ($kind, $path, $sort) = $self->_rearrange([qw(KIND PATH SORT)], @args);
    
  my @constructs = $self->get_constructs(-kind => $kind);
  my @objects = ();
  foreach my $construct (@constructs)
  {
    my $name = $construct->id();
    my $seq  = $construct->sequence();
    push @objects, Bio::Seq->new(-id => $name, -seq => $seq);
  }
  if ($sort && $sort eq 'size')
  {
    @objects = sort {length $b->seq <=> length $a->seq} @objects;
  }
  elsif ($sort && $sort eq 'name')
  {
    @objects = sort {$a->id cmp $b->id} @objects;
  }
  if ($path)
  {
    open my $OUTFH, '>', $path;
    my $FOUT = Bio::SeqIO->new(-fh => $OUTFH, -format => 'fasta');
    $FOUT->write_seq($_) foreach (@objects);
    close $OUTFH;
    return 1;
  }
  else
  {
    my $str = q{};
    foreach my $obj (@objects)
    {
      $str .= '>' . $obj->id . "\n" . $obj->seq . "\n";
    }
    return $str;
  }
}

=head2 compress_constructs

=cut

sub compress_constructs
{
  my ($self) = @_;

  my @compressibles = qw(building_block assembly_oligo stitching_oligo intermediate);
  my %switches;
  foreach my $kind (@compressibles)
  {
    my @constructs = $self->get_constructs(-kind => $kind);
    my %seqhsh;
    foreach my $construct (@constructs)
    {
      my $name = $construct->id();
      my $seq = $construct->sequence();
      if (! exists $seqhsh{$seq})
      {
        $seqhsh{$seq} = $name;
      }
      else
      {
        $switches{$name} = $seqhsh{$seq};
        $self->{cons}->removeChild($construct->xmlnode());
      }
    }
  }
  foreach my $oldname (keys %switches)
  {
    my $rulesource = '//subconstruct[. = "' . $oldname . q{"]};
    my @subnodes = $self->{cons}->findnodes($rulesource);
    foreach my $subconstruct (@subnodes)
    {
      my @whats = $subconstruct->childNodes;
      foreach my $what (@whats)
      {
        my $piecename = $what->data;
        if (exists $switches{$piecename})
        {
          $what->setData($switches{$piecename});
        }
      }
    }
  }
  return keys %switches;
}

=head1 POOL METHODS

=head2 get_pools

=cut

sub get_pools
{
  my ($self, @args) = @_;

  my ($sid, $method) = $self->_rearrange([qw(SUBID METHOD)], @args);
      
  my $rulesource = '//pool';
  if ($sid || $method)
  {
    $rulesource .= q{[};
  }
  if ($sid)
  {
    $rulesource .= q{subconstruct = "} . $sid . q{"};
  }
  if ($method && $sid)
  {
    $rulesource .= q{ and };
  }
  if ($method)
  {
    $rulesource .= q{@method = "} . $method . q{"};
  }
  if ($method || $sid)
  {
    $rulesource .= q{]};
  }
  my @nodes = $self->{root}->findnodes($rulesource);
  return @nodes;
}

=head1 WARNING METHODS

=head2 create_warning

=cut

sub create_warning
{
  my ($self, $atts) = @_;
  my $warning = Bio::DNAssemble::Warning->new(
    -doc => $self->{doc},
    -attributes => $atts,
  );
  $self->{warnings}->appendChild($warning->xmlnode());
  return $warning;
}

=head2 add_warning

=cut

sub add_warning
{
  my ($self, $warning) = @_;
  my $node = $warning->xmlnode();
  $self->{doc}->importNode($node);
  $self->{warnings}->appendChild($node);
  return 1;
}

=head2 warn

=cut

sub warn
{
  my ($self, @args) = @_;

  my ($construct, $note) = $self->_rearrange([qw(construct note)], @args);
  
  my $cid = $construct->id();
  my @warnings = $self->get_warnings(-sid => $cid);
  my $warning = $warnings[0];
  if (! defined $warning)
  {
    $warning = $self->create_warning({sid => $cid, reason => $note});
  }
  else
  {
    $warning->reason([$note]);
  }
  return;
}

=head2 get_warnings

=cut

sub get_warnings
{
  my ($self, @args) = @_;
  my ($sid) = $self->_rearrange([qw(SID)], @args);
  my $rulesource = '//warning';
  if ($sid)
  {
    $rulesource .= q{[@sid = "} . $sid . q{"]};
  }
  my @nodes = $self->{root}->findnodes($rulesource);
  my @warnings;
  foreach my $warning (@nodes)
  {
    push @warnings, Bio::DNAssemble::Error->new(-node => $warning);
  }
  return @warnings;
}

=head2 warnings_as_hash

=cut

sub warnings_as_hash
{
  my ($self) = @_;
  my @warnings = $self->{root}->findnodes(q{//warning});
  my %results = ();
  foreach my $warning (@warnings)
  {
    my $sid = $warning->getAttribute('sid');
    my $text = q{};
    my @sublist = $warning->findnodes(q{reason});
    my @reasons = map {$_->textContent} @sublist;
    $text .= join q{; }, @reasons;
    $results{$sid} = $text;
  }
  return \%results;
}

=head1 ERROR METHODS

=head2 create_error

=cut

sub create_error
{
  my ($self, $atts) = @_;
  my $error = Bio::DNAssemble::Error->new(
    -doc => $self->{doc},
    -attributes => $atts,
  );
  $self->{errors}->appendChild($error->xmlnode());
  return $error;
}

=head2 add_error

=cut

sub add_error
{
  my ($self, $error) = @_;
  my $node = $error->xmlnode();
  $self->{doc}->importNode($node);
  $self->{errors}->appendChild($node);
  return 1;
}

=head2 weep

=cut

sub weep
{
  my ($self, @args) = @_;

  my ($construct, $note) = $self->_rearrange([qw(construct note)], @args);
  
  my $cid = $construct->id();
  my @errors = $self->get_errors(-sid => $cid);
  my $error = $errors[0];
  if (! defined $error)
  {
    $error = $self->create_error({sid => $cid, reason => $note});
  }
  else
  {
    $error->reason([$note]);
  }
  return;
}

=head2 get_errors

=cut

sub get_errors
{
  my ($self, @args) = @_;
  my ($sid) = $self->_rearrange([qw(SID)], @args);
  my $rulesource = '//error';
  if ($sid)
  {
    $rulesource .= q{[@sid = "} . $sid . q{"]};
  }
  my @nodes = $self->{root}->findnodes($rulesource);
  my @errors;
  foreach my $error (@nodes)
  {
    push @errors, Bio::DNAssemble::Error->new(-node => $error);
  }
  return @errors;
}

=head2 errors_as_hash

=cut

sub errors_as_hash
{
  my ($self) = @_;
  my @errors = $self->{root}->findnodes(q{//error});
  my %results = ();
  foreach my $error (@errors)
  {
    my $sid = $error->getAttribute('sid');
    my $text = q{};
    my @sublist = $error->findnodes(q{reason});
    my @reasons = map {$_->textContent} @sublist;
    $text .= join q{; }, @reasons;
    $results{$sid} = $text;
  }
  return \%results;
}


=head1 ISSUE METHODS

=head2 create_issue

=cut

sub create_issue
{
  my ($self, $atts) = @_;
  my $issue = Bio::DNAssemble::Issue->new(
    -doc => $self->{doc},
    -attributes => $atts,
  );
  $self->{issues}->appendChild($issue->xmlnode());
  return $issue;
}

=head2 add_issue

=cut

sub add_issue
{
  my ($self, $issue) = @_;
  my $node = $issue->xmlnode();
  $self->{doc}->importNode($node);
  $self->{issues}->appendChild($node);
  return 1;
}

=head2 get_issues

=cut

sub get_issues
{
  my ($self, @args) = @_;
  my ($sid) = $self->_rearrange([qw(SID)], @args);
  my $rulesource = '//issue';
  if ($sid)
  {
    $rulesource .= q{[@sid = "} . $sid . q{"]};
  }
  my @nodes = $self->{root}->findnodes($rulesource);
  my @issues;
  foreach my $issue (@nodes)
  {
    push @issues, Bio::DNAssemble::Issue->new(-node => $issue);
  }
  return @issues;
}

=head2 remove_issue

=cut

sub remove_issue
{
  my ($self, $issue) = @_;
  my $node = $issue->xmlnode();
  $self->{issues}->removeChild($node);
  return 1;
}

=head2 critique

=cut

sub critique
{
  my ($self, @args) = @_;

  my ($construct, $attributes) = $self->_rearrange([qw(construct attributes)], @args);
  my $cid = $construct->id();
  $attributes->{sid} = $cid;
  my $error = $self->create_issue($attributes);
  return;
}

=head1 ACCESSOR METHODS

Methods for setting and accessing design attributes

=head2 name

=cut

sub name
{
  my ($self) = @_;
  return $self->{root}->getAttribute('id') || q{};
}

=head2 filename

=cut

sub filename
{
  my ($self) = @_;
  return $self->{filename};
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
