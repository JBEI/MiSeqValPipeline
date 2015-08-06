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
use English qw( -no_match_vars );
use Carp;
use Env qw(SCRATCH);

use Bio::Annotation::Comment;
use Bio::GeneDesign;
use EnvironmentModules;

use strict;
use warnings;

our $VERSION = '1.00';

module('load mysql');
module('load EMBOSS');
module('load vmatch');
module('load blast+');
module('load perl');
module('load python');
module('load genedesign');
module('load openmpi');
module('load taskfarmermq/2.1');

local $OUTPUT_AUTOFLUSH = 1;

#
# Get command line arguments
#
my %p = ();
GetOptions (
      'help'            => \$p{HELP},
      'batch=i'         => \$p{BATCH},
      'input=s'         => \$p{INPUT},
      'slots=i'         => \$p{SLOTS},
      'output=s'        => \$p{OUTPUT},

      'bblength:i'      => \$p{TARBBLEN},
      'bbmaxlength:i'   => \$p{MAXBBLEN},
      'bbminlength:i'   => \$p{MINBBLEN},
      'bblap:i'         => \$p{TARBBLAP},
      'cloningvector=s' => \$p{CLONVEC},
      'destvector=s'    => \$p{DESTVEC},
      'stitch:i'        => \$p{STITCH},

      'ollength:i'      => \$p{TAROLLEN},
      'olmaxlength:i'   => \$p{MAXOLLEN},
      'olminlength:i'   => \$p{MINOLLEN},
      'ollaptemp:i'     => \$p{LAPTEMP},
      'poolsize:i'      => \$p{POOLSIZE},
      'maxpoolnum:i'    => \$p{MAXPOOLNUM},
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
  croak "\n JGISB_ERROR: You must supply an input file.\n";
}
my ($iter, $filename, $suffix) = $GD->import_seqs($p{INPUT});


#
# The output path must exist, and we'll need it to end with a slash
#
$p{OUTPUT} = $p{OUTPUT} || q{.};
$p{OUTPUT} = $GD->endslash($p{OUTPUT});
if (! -e $p{OUTPUT})
{
  croak "\n JGISB_ERROR: $p{OUTPUT} does not exist.\n";
}


#
# The batch name defaults to 001; any input is padded to three digits
#
$p{BATCH} = $p{BATCH} || 1;
my $batch = 'Batch' . $GD->pad($p{BATCH}, 3, '0');

#
# Slots requested defaults to 8. I'm not going to do a lot of policing here.
#
$p{SLOTS} = $p{SLOTS} || 16;

#
# The temporary file path processing will be written to
#
my $tmpdir = $SCRATCH || $p{OUTPUT};
$tmpdir = $GD->endslash($tmpdir);
my $procdir = $tmpdir . $batch . q{/};
mkdir $procdir unless (-e $procdir);
my $procprefix = $procdir . $filename;
print "Will process in $procdir and write in $p{OUTPUT}\n";


################################################################################
################################# CONFIGURING ##################################
################################################################################
#
# Figure out what the BB and OL commands look like
#
my $ofstr = "--output $procdir --format genbank";

my %bopts = (
  MAXBBLEN => '--maxlength',    TARBBLEN => '--length',   TARBBLAP => '--lap',
  MINBBLEN => '--minlength',    STITCH => '--stitch',
  CLONVEC => '--cloningvector', DESTVEC => '--destvector'
);
my %oopts = (
  TAROLLEN => '--olilen', MAXOLLEN => '--maxolilen', MINOLLEN => '--minolilen',
  LAPTEMP => '--laptemp', POOLSIZE => '--poolsize', MAXPOOLNUM => '--maxpoolnum'
);

my $bcmdstd = $ofstr;
foreach my $key (sort grep {defined $p{$_}} keys %bopts)
{
  $bcmdstd .= q{ } . $bopts{$key} . q{ } . $p{$key};
}

my $ocmdstd = $ofstr;
foreach my $key (sort grep {defined $p{$_}} keys %oopts)
{
  $ocmdstd .= q{ } . $oopts{$key} . q{ } . $p{$key};
}


#
# Collect all the objects, anonymize their names by the batch name
#
my $count = 0;
my @objs = ();
while ( my $obj = $iter->next_seq() )
{
  $count++;
  my $oldid = $obj->id;
  my $newid = $batch . '_p' . $GD->pad($count, 3, '0');
  $obj->id($newid);
    
  my $collec = $obj->annotation;
  my $comment = Bio::Annotation::Comment->new();
  $comment->text("original_name = $oldid");
  $collec->add_Annotation('comment', $comment);
  $obj->annotation($collec);
  
  push @objs, $obj;
}


################################################################################
################################## PROCESSING ##################################
################################################################################
#
# Split the original input and write them out to process longest first
#
my $btasklist;
my $otasklist;
my %paths = ();
foreach my $obj (sort {$b->length <=> $a->length} @objs)
{
  my $id = $obj->id;
  my $prefix = $procdir . $id;

  #
  # Write a single object to the processing directory
  #
  my $ofilename = $id . '.genbank';
  my $filepath = $procdir . $ofilename;
  $GD->export_seqs(
    -filename => $ofilename,
    -path => $procdir,
    -sequences => [$obj]
  );
  
  #
  # Customize the BB and OL commands and determine what the output files are
  #
  my $bcmd = "GD_Design_Building_Blocks.pl $bcmdstd --input $filepath";
  $paths{$id}->{BBOUT} = $prefix . '_BB.genbank';
  $paths{$id}->{BBREP} = $prefix . '_BBreport.txt';

  my $ocmd = "GD_Design_Oligos.pl $ocmdstd --input " . $paths{$id}->{BBOUT};
  $paths{$id}->{OLOUT} = $prefix . '_BB_OL.genbank';
  $paths{$id}->{OLREP} = $prefix . '_BB_OLreport.txt';


  #
  # Determine the tasks required and add them to the list
  #
  my $btask = $bcmd . q{:};
  $btask .= $paths{$id}->{BBOUT} . q{,} . $paths{$id}->{BBREP};
  $btasklist .= $btask . ":0\n";

  my $otask = $ocmd . q{:};
  $otask .= $paths{$id}->{OLOUT} . q{,} . $paths{$id}->{OLREP};
  $otasklist .= $otask . ":0\n";
}


#
# Create the queue files that list the tasks for each queue
#
my $bqueue = $filename . '_q1';
my $qbfile = $procprefix . '_tasklist1.txt';
_filedump($qbfile, $btasklist);

my $oqueue = $filename . '_q2';
my $qofile = $procprefix . '_tasklist2.txt';
_filedump($qofile, $otasklist);


#
# Create worker request files that create N workers for the two queues.
#
my $shopen = "\#!/bin/bash -l\nmpirun -n $p{SLOTS} ";

my $wbfile = $procprefix . '_BBworker.sh';
my $wbcmd = $shopen . "tfmq-worker -q $bqueue\n";
_filedump($wbfile, $wbcmd);

my $wofile = $procprefix . '_OLworker.sh';
my $wocmd = $shopen . "tfmq-worker -q $oqueue\n";
_filedump($wofile, $wocmd);


#
# Submit the request for the BB workers and capture the sge job id.
# Submit the request for the OL workers to wait for the BB workers to finish.
#
my $sgeflags = "-l ram.c=5.25G -pe pe_fill $p{SLOTS} -cwd -V";
my $subregex = qr{Your job (\d+) \(\"}ms;

my $bjid = undef;
my $bjname = $filename . "_BB";
my $wbsub = "qsub $sgeflags -N $bjname $wbfile";
my $wbparse = _safeopen($wbsub);
$bjid = $1 if ($wbparse =~ $subregex);

my $ojid = undef;
my $ojname = $filename . "_OL";
my $wosub = "qsub -hold_jid $bjid $sgeflags -N $ojname $wofile";
my $woparse = _safeopen($wosub);
$ojid = $1 if ($woparse =~ $subregex);


#
# Start the two clients, one for the BB queue and one for the OL queue
# Capture the job ids.
#
my $cbcmd = "tfmq-client -i $qbfile -q $bqueue";
my $cbid = open2(my $cbout, my $cbin, $cbcmd);

my $cocmd = "tfmq-client -i $qofile -q $oqueue";
my $coid = open2(my $coout, my $coin, $cocmd);


################################################################################
################################### ASSEMBLY ###################################
################################################################################
#
# Wait for the OL workers to terminate and then compile the results
#
print "$count constructs gathered\n";
print $wbparse;
print $woparse;
print "Waiting for oligos...\n";
waitpid $coid, 0;

my @cats = ();
my @reps = ();
foreach my $obj (@objs)
{
  my $id = $obj->id;
  my $BREP = $paths{$id}->{BBREP};
  push @reps, $BREP if (-e $BREP);
  
  my $OREP = $paths{$id}->{OLREP};
  push @reps, $OREP if (-e $OREP);
  
  my $OOUT = $paths{$id}->{OLOUT};
  if (-e $OOUT)
  {
    push @cats, $OOUT;
  }
  else
  {
    my $BOUT = $paths{$id}->{BBOUT};
    push @cats, $BOUT if (-e $BOUT);
  }
}

my $outfile = $p{OUTPUT} . $batch . '_BB_OL.genbank';
print "Writing results to $outfile\n";
system "cat @cats > $outfile";

my $repfile = $p{OUTPUT} . $batch . '_report.txt';
print "Writing reports to $repfile\n";
system "cat @reps > $repfile";


my $JGISBV = "JGISB_DesignConstructs_$VERSION";
print "\n\n" . $GD->attitude() . " brought to you by $JGISBV\n\n";


exit;


sub _filedump
{
  my ($path, $data) = @_;
  open (my $FH, '>', $path) || croak ("can't write to $path, $OS_ERROR");
  print {$FH} $data;
  close $FH;
  return;
}

sub _safeopen
{
  my ($cmd) = @_;
  my ($inh, $outh) = (undef, undef);
  my $pid = open2($outh, $inh, $cmd) || croak ("oops on $cmd: $OS_ERROR");
  my $parse = <$outh>;
  return $parse;
}

__END__

=head1 NAME

JGISB_DesignConstructs.pl

=head1 VERSION

1.00

=head1 DESCRIPTION

Carves a set of chunks into building blocks and then chops those building blocks
into oligos. But all embarassingly parallel like.

=head1 USAGE

First load the modules
module load mysql
module load EMBOSS
module load vmatch
module load blast+
module load python
module load perl
module load taskfarmermq/2.1
module load genedesign

Submit a set of genbank files
qsub -N myjobname /path/to/JGISB_DesignConstructs.pl -i mychunks.genbank -o .

Submit a set of FASTA files to be Batch 31 and email me as the job progresses
qsub -N myjobname -M myemail@lbl.gov -m abe /path/to/JGISB_DesignConstructs.pl
  -i mychunks.FASTA -o . -b 31

=head1 REQUIRED ARGUMENTS

  --input : A sequence file to use as an input. This file may contain one
          or more construct sequences that are destined for building block and
          oligo design.

=head1 OPTIONS

  --output: Where should files be written
  --batch: What batch number should these constructs be anonymized to?
  --slots : How many slots should worker tasks request on genepool? default 8.
  --help : Display this message
  
Options for Building Blocks:

  --BBlength: the length in bp of building blocks (default is 1000)
  --BBmaxlength: the maximum length a building block is allowed to be.
  --BBminlength: the minimum length a building block is allowed to be.
  --BBlap: the target overlap between building blocks. (default is 40)
  --cloningvector: A vector in config/vectors that is to serve as the
      cloning vector
  --destvector: A vector in config/vectors that is to serve as the
      destination vector
  --stitch: The target melting temperature of assembly primers for
      the assembly of building blocks.

Options for Oligos:

  --OLlaptemp : Tm in temperature degrees C of oligo overlaps (default 70)
  --OLlength : the length in bp of assembly oligos (default 180)
  --OLminlength : minimum length of assembly oligos permitted (default 45)
  --OLmaxlength : maximum length of assembly oligos permitted (default 200)
  --poolsize : oligos will be pooled; GD will make bridging primers between
                pools of the specified size
  --maxpoolnum : the maximum number of pools to create.

=head1 DEPENDENCIES

If you are running this on genepool, you must load the modules mysql, EMBOSS,
vmatch, blast+, perl, and genedesign first. If you didn't do that, I am
surprised that you can read this.

This application requires taskfarmermq v2.1, which requires the python module.

=head1 BUGS AND LIMITATIONS

NO BUGS! HAHAHAHAHAHA!

=head1 INCOMPATIBILITIES

As long as you use it EXACTLY like I showed you, and on genepool, there are no
incompatibilities.

=head1 AUTHOR

Sarah Richardson
SMRichardson@lbl.gov

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2013, JGI Syn Bio developers
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
