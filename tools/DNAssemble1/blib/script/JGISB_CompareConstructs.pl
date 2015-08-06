#!/usr/bin/env perl

#
# SGE commands
#
#$ -l h_rt=14400
#$ -l ram.c=5G
#$ -pe pe_slots 1
#$ -cwd
#$ -V
#$ -P gentech-rnd.p

use Getopt::Long;
use Pod::Usage;
use IPC::Open2;
use CGI qw(:standard);
use English qw( -no_match_vars );
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
use Env qw(SCRATCH);

use Bio::GeneDesign;
use Bio::SearchIO;
use EnvironmentModules;

use strict;
use warnings;

our $VERSION = '1.00';

module('load blast+');
module('load perl');
module('load python');
module('load openmpi');
module('load taskfarmermq/2.2');

local $OUTPUT_AUTOFLUSH = 1;

#
# Get command line arguments
#
my %p = ();
GetOptions (
      'help'            => \$p{HELP},
      'input=s'         => \$p{INPUT},
      'slots=i'         => \$p{SLOTS},
      'threshold=i'     => \$p{THRESH},
);
if ($p{HELP})
{
  pod2usage(-verbose=>99);
}


################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $GD = Bio::GeneDesign->new();

#
# The input file must exist and be a format we care to read.
#
if (! $p{INPUT})
{
  die "\n JGISB_ERROR: You must supply an input file.\n";
}
my ($iter, $filename, $suffix) = $GD->import_seqs($p{INPUT});

#
# Slots requested defaults to 16. I'm not going to do a lot of policing here.
#
$p{SLOTS} = $p{SLOTS} || 16;

#
# Warning for similarity threshold at 80 by default.
#
$p{THRESH} = $p{THRESH} || 80;
if ($p{THRESH} <= 0 || $p{THRESH} >= 100)
{
  die "\n JGISB_ERROR: The threshold for similariy must be between 0 and 100.\n";
}

#
# The temporary file path processing will be written to
#
my $tmpdir = $SCRATCH;
$tmpdir = $GD->endslash($tmpdir);
my $procdir = $tmpdir . $filename . q{_ps/};
mkdir $procdir unless (-e $procdir);
my $procprefix = $procdir . $filename;
print "Will process in $procdir\n";


################################################################################
################################# CONFIGURING ##################################
################################################################################
#
# Collect all the objects, anonymize their names, split their sequences to file
#
my @registry = ();
my $count = 0;
my @objs = ();
my %seqobjs = ();
my %sequences = ();
my %seqpaths = ();
my %seqlens = ();
my %compares;
my %compareslinks = ();
my %key;
while ( my $obj = $iter->next_seq() )
{
  $count++;
  my $oldid = $obj->id;
  my $newid = 'p' . $GD->pad($count, 3, '0');
  $obj->id($newid);
  $sequences{$newid} = $obj->seq;
  $seqlens{$newid} = length $obj->seq;
  $seqobjs{$newid} = $obj;
  $key{$newid} = $oldid;

  #
  # Write a single object to the processing directory
  #
  my $ofilename = $newid . '.fasta';
  my $filepath = $procdir . $ofilename;
  $seqpaths{$newid} = $filepath;
  $GD->export_seqs(
    -filename => $ofilename,
    -path => $procdir,
    -sequences => [$obj],
    -format => 'fasta'
  );
  push @objs, $obj;
  push @registry, $filepath;
}
my @names = sort keys %sequences;

################################################################################
################################## COMPARING ###################################
################################################################################
#
# Create the blast commands; only compare each sequence once
#
my $tasklist;
my @paths;
for my $x (0 ..  $count - 1)
{
  my $ida = $names[$x];
  my $seqapath = $seqpaths{$ida};
  for my $y ($x + 1 .. $count - 1)
  {
    my $idb = $names[$y];
    next if ($idb eq $ida);
    my $seqbpath = $seqpaths{$idb};
    
    my $jid = $ida . "_" . $idb;
    my $oname = $procdir . "psrez_$jid.xml";
    my $cmd = "blastn -query $seqapath -subject $seqbpath -out $oname -outfmt 5";
    push @paths, $oname;
    $tasklist .= $cmd . q{:} . $oname . ":0\n";
  }
}


#
# Create the queue file that lists the tasks
#
my $queue = $filename . '_q';
my $qfile = $procprefix . '_tasklist.txt';
_filedump($qfile, $tasklist);
push @registry, $qfile;
push @registry, $qfile . '.done';


#
# Create worker request file that creates N workers for the queue.
#
my $wfile = $procprefix . '_worker.sh';
my $wcmd = "\#!/bin/bash -l\nmpirun -n 8 tfmq-worker -b 5 -t 60 -q $queue\n";
#my $wcmd = "\#!/bin/bash -l\ntfmq-worker -b 5 -t 60 -q $queue\n";
_filedump($wfile, $wcmd);
push @registry, $wfile;

#
# Submit the request for the workers and capture the sge job id.
#
my $sgeflags = "-l ram.c=5.25G -pe pe_fill $p{SLOTS} -cwd -V";
#my $sgeflags = "-l ram.c=5.25G -pe pe_slots 16 -cwd -V";
my $subregex = qr{Your job (\d+) \(\"}ms;

my $jname = $filename . "_blast";
my $wsub = "qsub $sgeflags -N $jname $wfile";
my $wparse = _safeopen($wsub);


#
# Start the client, capture the system job id.
#
my $ccmd = "tfmq-client -i $qfile -q $queue";
my $cid = open2(my $cout, my $cin, $ccmd);


################################################################################
################################### ASSEMBLY ###################################
################################################################################
#
# Wait for the workers to terminate and then compile the results
#
print "$count constructs gathered for comparison\n";
print $wparse;
print "Waiting for blast...\n";
waitpid $cid, 0;


#
# Gather all the blast output files and parse
#
my $zip = Archive::Zip->new();
my $zipname = $filename . '_archive.zip';
opendir (my $PDIR, $procdir) || die ("can't opendir $procdir, $OS_ERROR");
my @rezzes = grep { $_ =~ m{psrez_}} readdir $PDIR;
my $reportstring = q{};

foreach my $res (@rezzes)
{
  my $respath = $procdir . q{/} . $res;
  push @registry, $respath;
  my @precomp = split m{\.}, $res;
  my @comps = split m{\_}, $precomp[0];
  my ($firstid, $secondid) = ($comps[1], $comps[2]);
  my $firstlen = $seqlens{$firstid};
  my $secondlen = $seqlens{$secondid};
  if ($firstlen < $secondlen)
  {
    ($firstlen, $secondlen) = ($secondlen, $firstlen);
    ($firstid, $secondid) = ($secondid, $firstid);
  }
  my ($first, $second) = ($key{$firstid}, $key{$secondid});
  my $hitlen = 0;
  my $rez = Bio::SearchIO->new(-file => $respath, -format => 'blastxml');
  while (my $result = $rez->next_result())
  {
    while (my $hit = $result->next_hit() )
    {
      if (scalar($hit->hsps()))
      {
        my $thislen = $hit->length();
        my $mismatch = $thislen - $hit->matches('id');
        $hitlen += $thislen - $mismatch;
      }
    }
  }
  my $perc = sprintf "%.0f", 100 * ($hitlen / $secondlen);
  $perc = 100 if ($perc > 100);
  $compares{$firstid}->{$secondid} = $perc;
  $compares{$secondid}->{$firstid} = $perc;
  if ($perc >= $p{THRESH})
  {
    $reportstring .= "$key{$firstid} too similar to $key{$secondid} ($perc%)\n";
    my $graph = $GD->make_dotplot(
      -first      => $seqobjs{$firstid},
      -second     => $seqobjs{$secondid},
      -window     => 10,
      -stringency => 10
    );
    my $pngname = $precomp[0] . '.png';
    my $path = $procdir . $pngname;
    $compareslinks{$firstid}->{$secondid} = $pngname;
    $compareslinks{$secondid}->{$firstid} = $pngname;
    push @registry, $path;
    open   (my $IMG, '>', $path) || die ("JGISB_ERROR: Cannot write dot plot: $OS_ERROR");
    binmode $IMG;
    print   $IMG $graph;
    close   $IMG;
    $zip->addFile($path, $pngname);
    $zip->addFile($respath, $precomp[0] . '.xml');
  }
}

################################################################################
################################## REPORTING ###################################
################################################################################
my %COLORS = (
  100 => "9900000",
   90 => "#a31919",
   80 => "#ad3232",
   70 => "#b74c4c",
   60 => "#c16666",
   50 => "#cc7f7f",
   40 => "#d69999",
   30 => "#e0b2b2",
   20 => "#eacccc",
   10 => "#f4e5e5",
    0 => "#ffffff",
);
my @ranges = keys %COLORS;

my $style = <<"END";
<!--
td,th
{
  text-align: center;
}
-->
END

my $html = start_html(-style => {-code => $style});
$html .= "<table border=\"1\" align=\"center\">\n";
$html .= "<tr>\n\t<td>color key (% similarity)</td>\n";
foreach my $status (sort {$a <=> $b} keys %COLORS)
{
  next if $status == 100;
  my $bgcolor = $COLORS{$status};
  my $range = $status . q{-};
  $range .= $status + 10;
  $html .= "<td bgcolor=\"$bgcolor\">$range</td>\n";
}
$html .= "</tr>\n";
$html .= "</table>\n";
$html .= "<table border=\"1\" align=\"center\">\n";
$html .= "<tr>\n\t<td>Seq ID</td>\n";
foreach my $name (@names)
{
  $html .= "\t<th>$name</th>\n";
}
$html .= "</tr>\n";
foreach my $namer (@names)
{
  $html .= "<tr>\n";
  $html .= "\t<th>$namer</th>\n";
  foreach my $namec (@names)
  {
    if ($namec eq $namer)
    {
      my $color = $COLORS{100};
      $html .= "\t<td bgcolor=\"$color\">100</td>\n";
      last;
    }
    my $val = $compares{$namer}->{$namec};
    my $bin = bin(\@ranges, $val);
    my $color = $COLORS{$bin};
    $html .= "\t<td bgcolor=\"$color\">";
    if (exists $compareslinks{$namer}->{$namec})
    {
      $html .= "<a href=\"" . $compareslinks{$namer}->{$namec} . "\">";
      $html .= "$val</a></td>\n";
    }
    else
    {
      $html .= "$val</td>\n";
    }
  }
  $html .= "</tr>";
}
$html .= "</table>";
$html .= "<table border=\"1\" align=\"center\">\n";
$html .= "<tr>\n\t<th>Code</th><th>Name</th></tr>\n";
foreach my $name (@names)
{
  my $code = $key{$name};
  $html .= "<tr><td>$name</td><td>$code</td></tr>\n";
}
$html .= "</table>\n";

$html .= end_html();
my $htmlpath = $procdir . $filename . '.html';
open (my $HTML, '>', $htmlpath) || die "Can't open $htmlpath: $OS_ERROR\n";
print $HTML $html;
close $HTML;
$zip->addFile($htmlpath, $filename . '.html');
push @registry, $htmlpath;

my $reportpath = $procdir . $filename . '.txt';
open (my $REPORT, '>', $reportpath) || die "Can't open $reportpath: $OS_ERROR\n";
print $REPORT $reportstring;
close $REPORT;
$zip->addFile($reportpath, $filename . '.txt');
push @registry, $reportpath;

$zip->writeToFileNamed($zipname) == AZ_OK || die "JGISB_ERROR: Zip write error $OS_ERROR";

foreach my $path (@registry)
{
  unlink $path;
}
rmdir $procdir;

exit;

sub bin
{
  my ($range, $value) = @_;
  my @ranges = sort {$a <=> $b} @{$range};
  my $rangesize = scalar @ranges;
  for my $x (1 .. $rangesize - 1)
  {
    if ($value >= $ranges[$x-1] && $value <= $ranges[$x])
    {
      return $ranges[$x-1];
    }
  }
  return $ranges[-1];
}

sub _filedump
{
  my ($path, $data) = @_;
  open (my $FH, '>', $path) || die ("can't write to $path, $OS_ERROR");
  print {$FH} $data;
  close $FH;
  return;
}

sub _safeopen
{
  my ($cmd) = @_;
  my ($inh, $outh) = (undef, undef);
  my $pid = open2($outh, $inh, $cmd) || die ("oops on $cmd: $OS_ERROR");
  my $parse = <$outh>;
  return $parse;
}

__END__

=head1 NAME

JGISB_PrescreenConstructs.pl

=head1 VERSION

1.00

=head1 DESCRIPTION

Inspects a set of sequences and provides you with a graphical overview of their
pairwise likenesses. Dotplots are made for those over a certain similarity
threshold.

=head1 USAGE

First load the modules
module load blast+
module load python
module load perl
module load taskfarmermq/2.1

Submit a set of genbank files
qsub -N myjobname JGISB_PrescreenConstructs.pl -i mychunks.fasta

Submit a set of FASTA files to be Batch 31 and email me as the job progresses
qsub -N myjobname -M myemail@lbl.gov -m abe JGISB_PrescreenConstructs.pl
  -i mychunks.FASTA

=head1 REQUIRED ARGUMENTS

  --input : A sequence file to use as an input. This file may contain one
          or more construct sequences that are destined for building block and
          oligo design.

=head1 OPTIONS

  --slots : How many slots should worker tasks request on genepool? default 8.
  --help : Display this message

=head1 DEPENDENCIES

If you are running this on genepool, you must load the module perl first. If you
didn't do that, I am surprised that you can read this.

This application requires taskfarmermq v2.2, which requires the python module.

=head1 BUGS AND LIMITATIONS

NO BUGS! HAHAHAHAHAHA!

=head1 INCOMPATIBILITIES

As long as you use it EXACTLY like I showed you, and on genepool, there are no
incompatibilities.

=head1 AUTHOR

Sarah Richardson
SMRichardson@lbl.gov

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2014, JGI Syn Bio developers
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
