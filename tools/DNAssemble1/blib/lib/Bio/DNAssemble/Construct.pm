#
# DNAssemble module for vector objects
#

=head1 NAME

Bio::DNAssemble::Construct

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

DNAssemble object that represents a synthesized piece of DNA

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::DNAssemble::Construct;

use Carp;
use Bio::DNAssemble::Pool;
use Bio::DNAssemble::Feature;

use strict;
use warnings;

use base qw(Bio::DNAssemble::XMLNode);

our $VERSION = 1.00;

my %ATTRS = map {$_ => 1} qw(id kind method deliverable orient start end);
my %TEXTS = map {$_ => 1} qw(seqlen sequence GCp barcode uid vector
                             screening_group istart ilen
);
my %TEXTM = map {$_ => 1} qw();

=head1 CONSTRUCTOR METHODS

=head2 new

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{attrs} = \%ATTRS;
  $self->{texts} = \%TEXTS;
  $self->{textm} = \%TEXTM;

  my ($doc, $atts, $node) = $self->_rearrange([qw(DOC ATTRIBUTES NODE)], @args);
  
  if ($doc)
  {
    $self->{doc} = $doc;
    $self->{node} = $self->{doc}->createElement('construct');
    if ($atts)
    {
      $self->set_attributes($atts);
    }
  }
  elsif ($node)
  {
    $self->{doc} = $node->ownerDocument();
    $self->{node} = $node;
  }
  return $self;
}

=head1 ASSEMBLY METHODS

=head2 add_pool

=cut

sub add_pool
{
  my ($self, $atts) = @_;
  my $rulesource = './subconstructs';
  my @nodes = $self->{node}->findnodes($rulesource);
  if (! defined $nodes[0])
  {
    $self->{subconstructs} = $self->{doc}->createElement('subconstructs');
    $self->{node}->appendChild($self->{subconstructs});
  }
  else
  {
    $self->{subconstructs} = $nodes[0];
  }
  my $pool = Bio::DNAssemble::Pool->new(
    -construct => $self,
    -attributes => $atts
  );
  $self->{subconstructs}->appendChild($pool->xmlnode);
  return $pool;
}

=head2 get_pools

=cut

sub get_pools
{
  my ($self, @args) = @_;

  my ($method) = $self->_rearrange([qw(METHOD)], @args);
  
  my $rulesource = './subconstructs/pool';
  if ($method)
  {
    $rulesource .= q{[@method = "} . $method . q{"]};
  }
  my @nodes = $self->{node}->findnodes($rulesource);
  my @pools;
  foreach my $pool (@nodes)
  {
    push @pools, Bio::DNAssemble::Pool->new(-node => $pool);
  }
  return @pools;
}

=head2 add_subconstructs

=cut

sub add_subconstructs
{
  my ($self, $arr) = @_;

  my $rulesource = './subconstructs';
  my @nodes = $self->{node}->findnodes($rulesource);
  if (! defined $nodes[0])
  {
    $self->{subconstructs} = $self->{doc}->createElement('subconstructs');
    $self->{node}->appendChild($self->{subconstructs});
  }
  else
  {
    $self->{subconstructs} = $nodes[0];
  }
  my @constructs = @{$arr};
  foreach my $construct (@constructs)
  {
    $self->{subconstructs}->appendTextChild('subconstruct', $construct->id());
  }
  return;
}

=head2 get_subconstructs

NOTE: NOT EFFICIENT

=cut

sub get_subconstruct_ids
{
  my ($self, @args) = @_;

  my ($kind) = $self->_rearrange([qw(KIND)], @args);
  my $prulesource = './subconstructs/pool/subconstruct';
  my @subidnodes = $self->{node}->findnodes($prulesource);
  my @features;
  my %idhsh;
  foreach my $textnode (@subidnodes)
  {
    $idhsh{$textnode->textContent}++;
  }
  return keys %idhsh;
  #my $srulesource = '//constructs';
  #if ($kind)
  #{
  #  $srulesource .= q{[@kind = "} . $kind . q{"]};
  #}
  #my @subnodes = $self->{doc}->findnodes($srulesource);
  #my @subconstructs;
  #foreach my $subnode (@subnodes)
  #{
  #  my $subconstruct = Bio::DNAssemble::Construct->new(-node => $subnode);
  #  next if (! exists $idhsh{$subconstruct->id});
  #  print "Fetching ", $subconstruct->id, "\n";
  #  push @subconstructs, $subconstruct;
  #}
  #return @subconstructs;
}

=head2 add_features

=cut

sub add_features
{
  my ($self, $arr) = @_;
  my $rulesource = './features';
  my @nodes = $self->{node}->findnodes($rulesource);
  if (! defined $nodes[0])
  {
    $self->{features} = $self->{doc}->createElement('features');
    $self->{node}->appendChild($self->{features});
  }
  else
  {
    $self->{features} = $nodes[0];
  }
  my @seqfeatures = @{$arr};
  foreach my $seqfeature (@seqfeatures)
  {
    my $feature = Bio::DNAssemble::Feature->new(
      -construct => $self,
      -seqfeature => $seqfeature,
    );
    $self->{features}->appendChild($feature->xmlnode);
  }
  return;
}

=head2 get_features

=cut

sub get_features
{
  my ($self, @args) = @_;

  my ($ptag) = $self->_rearrange([qw(PRIMARY_TAG)], @args);
  
  my $rulesource = './features/feature';
  if ($ptag)
  {
    $rulesource .= q{[@primary_tag = "} . $ptag . q{"]};
  }
  my @nodes = $self->{node}->findnodes($rulesource);
  my @features;
  foreach my $feature (@nodes)
  {
    push @features, Bio::DNAssemble::Feature->new(-node => $feature);
  }
  return @features;
}

=head1 ATTRIBUTE METHODS

Methods for setting and accessing construct attributes

=head2 sequence

=cut

sub sequence
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    return $self->_set_sequence($value);
  }
  return $self->get_attribute('sequence');
}

=head2 deliverable

=cut

sub deliverable
{
  my ($self, $value) = @_;
  my $key = 'deliverable';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 kind

=cut

sub kind
{
  my ($self, $value) = @_;
  my $key = 'kind';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 method

=cut

sub method
{
  my ($self, $value) = @_;
  my $key = 'method';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 barcode

=cut

sub barcode
{
  my ($self, $value) = @_;
  my $key = 'barcode';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 uid

=cut

sub uid
{
  my ($self, $value) = @_;
  my $key = 'uid';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 vector

=cut

sub vector
{
  my ($self, $value) = @_;
  my $key = 'vector';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 seqlen

=cut

sub seqlen
{
  my ($self, $value) = @_;
  my $key = 'seqlen';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 GCp

=cut

sub GCp
{
  my ($self, $value) = @_;
  my $key = 'GCp';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 istart

=cut

sub istart
{
  my ($self, $value) = @_;
  my $key = 'istart';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 start

=cut

sub start
{
  my ($self, $value) = @_;
  my $key = 'start';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 end

=cut

sub end
{
  my ($self, $value) = @_;
  my $key = 'end';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 ilen

=cut

sub ilen
{
  my ($self, $value) = @_;
  my $key = 'ilen';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 screening_group

=cut

sub screening_group
{
  my ($self, $value) = @_;
  my $key = 'screening_group';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 id

=cut

sub id
{
  my ($self, $value) = @_;
  my $key = 'id';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
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
