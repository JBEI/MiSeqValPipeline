#
# DNAssemble module for objects based on XML nodes
#

=head1 NAME

Bio::DNAssemble::Error

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

DNAssemble object that represents a synthesized piece of DNA

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::DNAssemble::XMLNode;

use Bio::GeneDesign::Basic qw(:GD);
use Carp;

use strict;
use warnings;

use base qw(Bio::Root::Root);

our $VERSION = 1.00;


=head1 CONSTRUCTOR METHODS

=head2 new

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
}

=head1 ATTRIBUTE METHODS

Methods for setting and accessing construct attributes

=head2 get_attribute

=cut

sub get_attribute
{
  my ($self, $key) = @_;
  return undef if (! defined $key);
  if (exists $self->{attrs}->{$key})
  {
    return $self->{node}->getAttribute($key);
  }
  elsif (exists $self->{texts}->{$key})
  {
    return $self->_fetch_from_textnode($key);
  }
  elsif (exists $self->{textm}->{$key})
  {
    my @values = $self->_list_textnodes($key);
    return join q{; }, @values;
  }
  else
  {
    $self->throw('Do not know what to do with attribute of type ' . $key);
  }
  return;
}

=head2 set_attributes

=cut

sub set_attributes
{
  my ($self, $atthsh) = @_;
  return if (! defined $atthsh);
  my %sethsh = %{$atthsh};
  foreach my $key (keys %sethsh)
  {
    if ($key eq 'sequence')
    {
      $self->_set_sequence($sethsh{$key});
    }
    elsif (exists $self->{attrs}->{$key})
    {
      $self->{node}->setAttribute($key, $sethsh{$key});
    }
    elsif (exists $self->{texts}->{$key})
    {
      $self->_set_to_textnode($key, $sethsh{$key});
    }
    elsif (exists $self->{textm}->{$key})
    {
      my $ref = ref $sethsh{$key};
      if ($ref eq 'ARRAY')
      {
        my @values = @{$sethsh{$key}};
        $self->_add_textnode($key, $_) foreach @values;
      }
      else
      {
        $self->_add_textnode($key, $sethsh{$key});
      }
    }
    else
    {
      $self->throw('Do not know what to do with attribute of type ' . $key);
    }
  }
  return;
}

=head2 _fetch_from_textnode

=cut

sub _fetch_from_textnode
{
  my ($self, $tagname) = @_;
  my @list = $self->{node}->getChildrenByTagName($tagname);
  if (! scalar @list)
  {
    return q{};
  }
  return $list[0]->textContent;
}

=head2 _list_textnodes

=cut

sub _list_textnodes
{
  my ($self, $tagname) = @_;
  my @list = $self->{node}->getChildrenByTagName($tagname);
  my @values = ();
  foreach my $node (@list)
  {
    push @values, $node->textContent;
  }
  return @values;
}

=head2 _set_to_textnode

=cut

sub _set_to_textnode
{
  my ($self, $tagname, $value) = @_;
  my @list = $self->{node}->getChildrenByTagName($tagname);
  if (! scalar @list)
  {
    $self->{node}->appendTextChild($tagname, $value);
  }
  else
  {
    my @subs = $list[0]->childNodes();
    $subs[0]->setData($value);
  }
  return;
}

=head2 _add_textnode

=cut

sub _add_textnode
{
  my ($self, $tagname, $value) = @_;
  $self->{node}->appendTextChild($tagname, $value);
  return;
}

=head2 _set_sequence

=cut

sub _set_sequence
{
  my ($self, $sequence) = @_;
  $sequence = uc $sequence;
  $self->_set_to_textnode('sequence', $sequence);
  my $scount = _count($sequence);
  $self->_set_to_textnode('GCp', $scount->{'GCp'});
  $self->_set_to_textnode('seqlen', $scount->{'len'});
  return;
}

=head1 ACCESSOR METHODS

=head2 xmlnode

=cut

sub xmlnode
{
  my ($self) = @_;
  return $self->{node};
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
