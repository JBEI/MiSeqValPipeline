#!/usr/bin/env perl

use Bio::GeneDesign;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_GB_EnzymeList_$VERSION";

$| = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'              => \$p{HELP},
      'input=s'           => \$p{INPUT},
      'output=s'          => \$p{OUTPUT}
);

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $GD = Bio::GeneDesign->new();
$GD->set_restriction_enzymes(-enzyme_set => "standard");
my $RES = $GD->enzyme_set;

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
my $count = 1;
my %enzes = ();

################################################################################
################################## ARRAYING  ###################################
################################################################################
while ( my $obj = $iter->next_seq() )
{
  my @seqfeats = $obj->get_SeqFeatures;
  foreach my $feature (@seqfeats)
  {
    if ($feature->has_tag("enz_t"))
    {
      my $enzname = join("", $feature->get_tag_values("enz_t"));
      $enzname =~ s/\s//g;
      $enzes{$enzname}++;
      $count++;
    }
    if ($feature->has_tag("enz_f"))
    {
      my $enzname = join("", $feature->get_tag_values("enz_f"));
      $enzname =~ s/\s//g;
      $enzes{$enzname}++;
      $count++;
    }
  }
}

my $total = scalar(keys %enzes);

my $enzfile = $p{OUTPUT} . $filename . "_enzymes.txt";
open (my $ENZ, ">" . $enzfile ) || die "can't write to $enzfile, $!";
print $ENZ "# Enzymes needed for $filename\n";
print $ENZ "$_\t$enzes{$_}\t" . $RES->{$_}->cutseq . "\n" foreach sort keys %enzes;
close $ENZ;


print "\nNeed $count instances of $total enzymes\n\n";
print "\nPrinted an enzyme file at $enzfile\n\n";

print "brought to you by $JGISBV\n\n";

exit;

__END__

=head1 NAME

  JGISB_GB_EnzymeList_FA.pl

=head1 VERSION

  Version 1.00

=head1 DESCRIPTION

List the enzymes needed in a genbank file

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