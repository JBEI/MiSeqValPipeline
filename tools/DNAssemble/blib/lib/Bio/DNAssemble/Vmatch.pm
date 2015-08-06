#
# DNAssemble using Vmatch
#

=head1 NAME

DNAssemble::Vmatch

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

Vmatch interface

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::DNAssemble::Vmatch;
require Exporter;

use Digest::MD5 qw(md5_hex);
use Carp;
use Bio::Grep;

use strict;
use warnings;

our $VERSION = 1.00;

use base qw(Exporter);
our @EXPORT_OK = qw(
  _filter_vmatch
  _search_vmatch
  $VERSION
);
our %EXPORT_TAGS =  (GD => \@EXPORT_OK);

=head1 Functions

=head2 _filter_vmatch()

=cut

sub _filter_vmatch
{
  my ($parent, $seqarrs, $mismatches, $revcom, $writedir) = @_;

  my $rand = md5_hex(md5_hex(time().{}.rand().$$));
  my $band = md5_hex(md5_hex(time().{}.rand().$$));

  #Write the parent to a FASTA file
  my $dbpath = $writedir . "$rand.FASTA";
  open (my $POUT, '>', $dbpath) || croak($!);
  print $POUT ">" . $parent->id . "\n" . $parent->seq . "\n";
  close $POUT;

  #Write the queries to a FASTA file
  my $datapath = $writedir . "$band.FASTA";
  open (my $DOUT, '>', $datapath) || croak($!);
  print $DOUT ">" . $_->id . "\n" . $_->seq . "\n" foreach (@$seqarrs);
  close $DOUT;

  #Make the vmatch database
  my $sbe = Bio::Grep->new('Vmatch');
  $sbe->settings->datapath($writedir);
  $sbe->generate_database({ file => $dbpath, description => $parent->id});

  #Run vmatch on everything in the array
  $sbe->search
  ({
    database            => "$rand.FASTA",
    query_file          => $datapath,
    complete            => 1,
    mismatches          => $mismatches,
    direct_and_rev_com  => 1,
    online              => 1}
  );
  my %hits = map {$_->id => 0} @$seqarrs;
  while ( my $res = $sbe->next_res )
  {
    $hits{$res->query->id}++;
  }
  #Clean up and return
  system "rm $writedir$rand.FASTA.*" || croak ("Can't clean up? $!");
  system "rm $dbpath" || croak ("Can't clean up? $!");
  system "rm $datapath" || croak ("Can't clean up? $!");
  my @cleans = grep {$hits{$_->id} == 1} @$seqarrs;
  return \@cleans;
}

=head2 _search_vmatch()

=cut

sub _search_vmatch
{
  my ($parent, $seq, $mismatches, $revcom, $writedir) = @_;

  my $rand = md5_hex(md5_hex(time().{}.rand().$$));
  my $band = md5_hex(md5_hex(time().{}.rand().$$));

  #Write the parent to a FASTA file
  my $dbpath = $writedir . "$rand.FASTA";
  open (my $POUT, '>', $dbpath) || croak($!);
  print $POUT q{>} . $parent->id . "\n" . $parent->seq . "\n";
  close $POUT;

  #Write the queries to a FASTA file
  my $datapath = $writedir . "$band.FASTA";
  open (my $DOUT, '>', $datapath) || croak($!);
  print $DOUT q{>} . $seq->id . "\n" . $seq->seq . "\n";
  close $DOUT;

  #Make the vmatch database
  my $sbe = Bio::Grep->new('Vmatch');
  $sbe->settings->datapath($writedir);
  $sbe->generate_database
  (
    {
      file => $dbpath,
      description => $parent->id
    }
  );

  #Run vmatch on everything in the array
  $sbe->search
  (
    {
      database            => "$rand.FASTA",
      query_file          => $datapath,
      complete            => 1,
      mismatches          => $mismatches,
      direct_and_rev_com  => 1,
      online              => 1
    }
  );
  my @hits;
  while ( my $res = $sbe->next_res )
  {
    push @hits, $res->sequence();
  }
  #Clean up and return
  system "rm $writedir$rand.FASTA.*" || croak ("Can't clean up? $!");
  system "rm $dbpath" || croak ("Can't clean up? $!");
  system "rm $datapath" || croak ("Can't clean up? $!");
  return \@hits;
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

