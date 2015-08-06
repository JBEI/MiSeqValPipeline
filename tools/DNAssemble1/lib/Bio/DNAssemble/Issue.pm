#
# DNAssemble module for error objects
#

=head1 NAME

Bio::DNAssemble::Issue

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

DNAssemble object that represents a synthesized piece of DNA

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::DNAssemble::Issue;

use Carp;

use strict;
use warnings;

use base qw(Bio::DNAssemble::XMLNode);

our $VERSION = 1.00;

my %ATTRS = map {$_ => 1} qw(sid number);
my %TEXTS = map {$_ => 1} qw(start end kind repeat);
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
    $self->{node} = $doc->createElement('issue');
    if ($atts)
    {
      $self->set_attributes($atts);
    }
  }
  elsif ($node)
  {
    $self->{node} = $node;
  }
  return $self;
}

=head1 ATTRIBUTE METHODS

Methods for setting and accessing error attributes

=head2 sid

=cut

sub sid
{
  my ($self, $value) = @_;
  my $key = 'sid';
  if (defined $value)
  {
    return $self->set_attributes({$key => $value});
  }
  return $self->get_attribute($key);
}

=head2 number

=cut

sub number
{
  my ($self, $value) = @_;
  my $key = 'number';
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


=head2 repeat

=cut

sub repeat
{
  my ($self, $value) = @_;
  my $key = 'repeat';
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

1;

__END__

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, DNAssemble developers
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
