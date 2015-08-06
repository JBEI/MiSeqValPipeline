#
# DNAssemble using EMBOSS palindrome
#

=head1 NAME

DNAssemble::Palindrome

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

EMBOSS Palindrome interface

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::DNAssemble::Palindrome;
require Exporter;

use Digest::MD5 qw(md5_hex);
use Bio::Factory::EMBOSS;
use Bio::Tools::EMBOSS::Palindrome;

use strict;
use warnings;

our $VERSION = 1.00;

use base qw(Exporter);
our @EXPORT_OK = qw(
  _filter_palindromes
  $VERSION
);
our %EXPORT_TAGS =  (GD => \@EXPORT_OK);

=head1 Functions

=head2 _filter_palindromes()

=cut

sub _filter_palindromes
{
  my ($seqarr, $minlen, $maxlen, $gaplimit, $mismatches, $writedir) = @_;

  #Write out everything in the seqarr to a FASTA file
  my $rand = md5_hex(md5_hex(time().{}.rand().$$));
  my $path = $writedir . "$rand.FASTA";
  open (my $FOUT, '>',  $path) || croak ("can't write tmp pal file: $!");
  print $FOUT ">" . $_->id . "\n" . $_->seq . "\n" foreach (@$seqarr);
  close $FOUT;
  #Run the palindrome program and write out to a .pal file
  my $factory = Bio::Factory::EMBOSS->new();
  my $app = $factory->program('palindrome');
  my $out = $writedir . "$rand.pal";
  $app->run({
    -sequence      => $path,
    -minpallen     => $minlen,
    -maxpallen     => $maxlen,
    -gaplimit      => $gaplimit,
    -outfile       => $out,
    -nummismatches => $mismatches
  });

  #Parse the .pal file and keep only sequences without secondary structure
  my $parser = Bio::Tools::EMBOSS::Palindrome->new(-file => $out);
  my %hpids;
  while( my $seq = $parser->next_seq )
  {
    my @hpins = $seq->get_SeqFeatures;
    $hpids{$seq->id}++ if (scalar @hpins);
    #foreach my $feat ( @hpins )
    #{
    #print $seq->id;
    #my ($a, $b, $c, $d) = ($feat->start + $seq->id, $feat->end + $seq->id,
    #  $feat->hstart + $seq->id, $feat->hend + $seq->id);
    #print " has a repeat from $a..$b and $c..$d\n\n";
    #}
  }

  #Clean up and return
  system "rm $out" || croak ("Can't clean up? $!");
  system "rm $path" || croak ("Can't clean up? $!");
  my @cleans = grep {! exists $hpids{$_->id}} @$seqarr;
  return \@cleans;
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

