#
# DNAssemble using Blast
#

=head1 NAME

DNAssemble::Blast

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

Vmatch interface

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::DNAssemble::Blast;
require Exporter;

use Bio::Tools::Run::StandAloneBlastPlus;

use strict;
use warnings;

our $VERSION = 1.00;

use base qw(Exporter);
our @EXPORT_OK = qw(
  _filter_blast
  $VERSION
);
our %EXPORT_TAGS =  (GD => \@EXPORT_OK);

=head1 Functions

=head2 _filter_blast()

=cut

sub _filter_blast
{
  my ($parent, $seqarrs, $percent, $writedir) = @_;

  #Make the BLAST factory
  my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
                          -program => 'blastn',
                          -db_dir  => $writedir,
                          -db_data => $parent);

  #BLAST everything in the array
  my %hits = map {$_->id => 0} @$seqarrs;
  my %objs = map {$_->id => $_} @$seqarrs;
  $factory->run(
     -method => "blastn",
      -query => $seqarrs,
      -method_args => [ -word_size => 4,
                        -perc_identity => $percent,
                        -gapextend => 2,
                        -gapopen => 0,
                        -penalty => -1]);


  while (my $result = $factory->next_result)
  {
    my $name = $result->query_name();
    my $qlen = length($objs{$name}->seq);
    while( my $hit = $result->next_hit())
    {
      while( my $hsp = $hit->next_hsp())
      {
        my $hlen = length($hsp->hit_string);
        $hits{$name}++ if ($hlen >= 0.7*$qlen && $hsp->evalue < 0.0001);
      }
    }
  }

  #Clean up and return
  $factory->cleanup;
  my @cleans = grep {$hits{$_->id} == 1} @$seqarrs;
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

