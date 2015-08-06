#
# DNAssemble module for vector objects
#

=head1 NAME

Bio::DNAssemble::Vector

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

DNAssemble object that represents a vector

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::DNAssemble::Vector;

use Bio::SeqIO;
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

  my ($name, $path) = $self->_rearrange([qw(NAME PATH)], @args);

  my $iterator = Bio::SeqIO->new(-file => $path)
    || $self->throw("Can't parse $path!");

  my @VECS;
  while ( my $obj = $iterator->next_seq() )
  {
    push @VECS, $obj;
  }
  $self->throw("Too many vectors in file at $path") if (scalar(@VECS) > 1);
  $self->throw("No vectors in file at $path") if (scalar(@VECS) == 0);

  $self->{name} = $name;
  $self->{seqobj} = $VECS[0];
  $self->{seq} = $VECS[0]->seq;
  $self->{len} = length $self->{seq};

  my @feats = $self->{seqobj}->get_SeqFeatures;
  $self->{features} = \@feats;


  my @es = grep {$_->primary_tag eq "protein_bind"} @feats;
  my @exes = map {join(q{}, $_->get_tag_values("label"))} @es;

  my @ps = grep {$_->primary_tag eq "primer_bind"} @feats;
  my %dsubs = map {join(q{}, $_->get_tag_values("label")) => $_} @ps;
  if (exists $dsubs{"CLH5"})
  {
    $self->{chew5} = uc $dsubs{"CLH5"}->seq->seq;
    $self->{chew5loc} = $dsubs{"CLH5"}->start;
  }
  else
  {
    carp ("$name has no CLH5 sequence annotated");
    $self->{chew5} = q{};
    $self->{chew5loc} = 0;
  }
  if (exists $dsubs{"CLH3"})
  {
    $self->{chew3} = uc $dsubs{"CLH3"}->seq->seq;
    $self->{chew3loc} = $dsubs{"CLH3"}->start;
  }
  else
  {
    carp ("$name has no CLH3 sequence annotated");
    $self->{chew3} = q{};
    $self->{chew3loc} = 0;
  }

  $self->{enzyme_list} = \@exes;

  return $self;
}

=head1 FUNCTIONAL METHODS

=head2 clone

=cut

sub clone
{
   my ($self) = @_;
   my $copy;
   foreach my $key (keys %$self)
   {
     if (ref $self->{$key} eq "ARRAY")
     {
       $copy->{$key} = [@{$self->{$key}}];
     }
     elsif (ref $self->{$key} eq "HASH")
     {
       $copy->{$key} = {%{$self->{$key}}};
     }
     else
     {
      $copy->{$key} = $self->{$key};
     }
   }
   bless $copy, ref $self;
   return $copy;
}

=head1 ACCESSOR METHODS

Methods for setting and accessing vector attributes

=head2 name

=cut

sub name
{
  my ($self) = @_;
  return $self->{'name'};
}

=head2 id

=cut

sub id
{
  my ($self) = @_;
  return $self->{'name'};
}

=head2 seq

=cut

sub seq
{
  my ($self) = @_;
  return $self->{'seq'};
}

=head2 len

=cut

sub len
{
  my ($self) = @_;
  return $self->{'len'};
}

=head2 seqobj

=cut

sub seqobj
{
  my ($self) = @_;
  return $self->{'seqobj'};
}

=head2 chew5

=cut

sub chew5
{
  my ($self) = @_;
  return $self->{'chew5'};
}

=head2 chew5loc

=cut

sub chew5loc
{
  my ($self) = @_;
  return $self->{chew5loc};
}

=head2 chew3

=cut

sub chew3
{
  my ($self) = @_;
  return $self->{'chew3'};
}

=head2 chew3loc

=cut

sub chew3loc
{
  my ($self) = @_;
  return $self->{chew3loc};
}

=head2 enzyme_list

=cut

sub enzyme_list
{
  my ($self) = @_;
  return @{$self->{'enzyme_list'}};
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
