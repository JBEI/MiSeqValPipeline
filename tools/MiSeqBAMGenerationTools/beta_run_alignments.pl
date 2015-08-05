#!/usr/bin/env perl
use warnings; 
use strict;
use Modern::Perl;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw/uniq/;
use Cwd;
use YAML qw/LoadFile/;
use GT::Config;
use GT::Clutils;

umask 2;

my ( $help, $align_dir, $sge_project, $ref_file, $aln_args, $sge_args, $samXe_args, 
	$no_db_for_rg, $bwa_algorithm, $picard_path, $samtools_path, $bwa_path, $verbose, $do_analysis_task_stuff, $reseq_bin_dir, $DEBUG );

$DEBUG         = 0;
$align_dir     = 'bwa_dir';
$bwa_algorithm = 'mem';

$do_analysis_task_stuff = 0;

my $lib_file = '';

#Location of fastq_slice, guess_qual_format, and fastq_single_or_paired command-line commands
$reseq_bin_dir = '';

#Path to Picard JAR files
$picard_path = "";
$samtools_path = "";
$bwa_path = "";

GetOptions(
	'basedir=s'                => \$align_dir,
	'ref=s'                    => \$ref_file,
	'config_file|libfile=s'    => \$lib_file,
	'Project=s'                => \$sge_project,
	'verbose'                  => \$verbose,
	'aln-args=s'               => \$aln_args,
	'sampe-args=s'             => \$samXe_args,
	'sge-args=s'               => \$sge_args,
	'algo-bwa=s'               => \$bwa_algorithm,
	'picard_path=s'				=> \$picard_path,
	'samtools_path=s'			=> \$samtools_path,
	'bwa_path=s'				=> \$bwa_path,
	'at'                       => \$do_analysis_task_stuff,
	'weird_fq_name'            => \$no_db_for_rg,
	'reseqbindir=s'            	=> \$reseq_bin_dir,
	'DEBUG'                    => \$DEBUG,
	'help'                     => \$help,
);
#print "modules = $modules\n"; die;
sub usage {
	my $bin = basename $0;
	print STDERR <<EOH;
$bin [options] -c libs.info
e.g. $bin -algo backtrack -modules -c libs.inf  \# backtrack (aln + sampe/samse) algorithm using the 'modules' version of bwa
     $bin -c libs.info                          \# mem algorithm ( best unless read length < 70  ), avoiding modules version of bwa
          *tyically run w/o options on directory prepared with setup_dirs
 -algo-bwa       STRING   mem (default), backtrack or bwasw
 -aln-args       STRING   args for bwa aln ( e.g. -aln-args \'-M 4 -O 14\' )
 -sampe-args     STRING   args for bwa sampe or samse ( e.g. -sampe-args \'-M 4 -O 14\' )
 -sge-args       STRING   args for sge ( e.g. -sge-ags \"-l high.c\" ***USE THEM QUOTES***
 -weird_fq_name           (-weird ok also) Don't parse fastq filename for RG tag information, don't lookup in database either
 -config         libinfo  If reads not in run-id.lane.* (  run-id.lane.basecall-id.barcode ) 
                          this file needed to map read names to library names and samples
                          *also used if reads aren't in venonat
 -base           DIR      name of subdirectory to run in
 -noat                    SKIP analysis task updating
 -P              STRING   sge accounting project name for qsub
 -ref            FILE.fa  if no dir/base/config.yml, reference file to align against
                          must be already be indexed by bwa
 -v                       verbose, print more and ask run_bwa to print more
EOH
}

if ( $help || ( ! @ARGV && ! ( defined $lib_file && -e $lib_file ) )
   || !defined $align_dir || ( defined $ref_file && !-e $ref_file ) ) {    # ref from config usually now
	usage();
	exit;
}


my $no_mod_option = '-nomodules'; 

my $set_ref_by_config;

if ($ref_file) {
	# TODO: add updating configs, ** also, check all configs before launching anything **
	$set_ref_by_config = 0;
}
else {
	$set_ref_by_config = 1;
}

sub get_ref_from_config {
	my $config = GT::Config->new();
	$config->load();

	$ref_file = $config->{ref_fasta};
	if ( !defined $ref_file || !-e $ref_file ) {
		warn "ERROR: no ref_fasta entry in config.yml";
		usage();
		exit;
	}
	return $ref_file;
}

my $project_option = defined $sge_project ? "-P $sge_project " : '';
my @cmds = ();
my %cmds = ();

my $start_dir = getcwd;

#
#  get RG information from a text file instead of databases
#
# lib file format is:
# id, lib, full_path_fq, gen, sp, strain, isolate
#   names are ignored at the moment though
# switch to getting from config.yml
#  - updated by split_fq
#

my %libs     = ();
my @dir_list = @ARGV;

if ( $lib_file && -e $lib_file ) {
	my %dir_list = ();
	my $order    = 1;
	open( my $lfh, $lib_file ) or die $!;
	while (<$lfh>) {
		next if ( /^#/ || /^\s*$/ );
		my ( $pmo_info, $lib, $fp_fq ) = split;

		#TODO: should be shared function with slice_fastq
		my $fq_name = truncate_fastq_filename($fp_fq);

		if ( $DEBUG ) { print STDERR "libs{ $fq_name } = $lib\n"; }

		$libs{$fq_name} = $lib;
		$dir_list{$pmo_info} = $order++;
	}
	push @dir_list, sort { $dir_list{ $a } <=> $dir_list{ $b } } keys %dir_list;
}

sub truncate_fastq_filename {
	my ($name) = @_;
	$name = basename $name;
	$name =~ s/[.]fastq.gz$//;
	$name =~ s/[.]fq.gz$//;
	return ($name) if ($name);
}

####
# sliced look like
# pe-2096.7.1751.2@176000001.fq.gz
# pe-2096.7.1751.ABCDEF.2@176000001.fq.gz
####


for my $dir ( grep { -d } uniq @dir_list ) {
	my $pinfo = basename $dir;
	if ( !-e "$dir/$align_dir" ) {
		warn "skipping $dir -- no $align_dir\n";
	}

	my @bwas;
	chdir "$dir/$align_dir" or die "can't cd to $dir/$align_dir $!";

	if ($set_ref_by_config) {
		$ref_file = get_ref_from_config;
	}

	# load config for read info ( 
	# *** combine with above?

	my $config = LoadFile( 'config.yml' );
	for my $cfg_fq ( keys %{ $config->{fqs } } ) {
		my $trunc_fq_name = truncate_fastq_filename( $cfg_fq );
		if ( ! exists $libs{ $trunc_fq_name } ) {
			print STDERR "uh adding $trunc_fq_name for $cfg_fq = $config->{fqs}{$cfg_fq}->{lib}:\n";
			$libs{ $trunc_fq_name } = $config->{fqs}{$cfg_fq}->{library};
		}
	}
	
	if ( !-e "bam_dir" ) {
		mkdir 'bam_dir' or die "ERROR: failed to create $dir/bam_dir $!";
	}

	my ( @peRGs, @peR1s, @peR2s );  # ReadGroups, Read1s, Read2s

	my %run_lanes = ();

	#
	# get information for paired reads
	#
	for my $read1 ( glob("fastq_dir/pe-*.1@[0-9]*.fq.gz") ) {
		my $read2 = $read1;
		# change .1@ to .2@  - a bit excessive on the exact matching, sorry if it makes it look messy.

		$read2 =~ s/[.]1\@(\d+)[.]fq[.]gz$/.2\@$1.fq.gz/;

		if ( !-e $read1 or !-e $read2 ) { die "ERROR: can't find  $read1 or maybe $read2 in $dir" }

		my ( $id, $lib );
		if ( $no_db_for_rg ) { 
			( $id ) = split /\@/, basename $read1;
			$id =~ s/^pe-//;
			$id =~ s/[.]1$//;
			if ( ! exists $libs{ $id } ) {
				die "no lib for $id ( from $read1 ) in $dir";
			}
			$lib = $libs{ $id };
		}
		else {

			my ($read_basename) = split /\@/, basename $read1;    # for filename matching from libs.info file
			$read_basename =~ s/^pe-//;                           # remove pairing indicator
			$read_basename =~ s/[.]1$//;                          # remove read1/2 indicator

			if( $DEBUG ) { print STDERR  "check for libs{ $read_basename }\n"; }

			if (exists $libs{$read_basename}) {
				$lib = $libs{$read_basename};
			}

			die "no lib for $read1($read_basename) in $dir" unless ( defined $lib );
			
			#B/c $libs{$read_basename} must exist, otherwise program would have died already
			$id = $read_basename;
		}

		my $rg = "\'\@RG\\tID:$id\\tSM:$pinfo\\tPL:illumina\\tLB:$lib\'";

		push @peRGs, $rg;
		push @peR1s, $read1;
		push @peR2s, $read2;

	}

	#
	# get information for single ended reads
	#
	my ( @seRGs, @seR1s );

	for my $read1 ( glob("fastq_dir/se-*") ) {
		if ( !-e $read1 ) { die "fail on $read1 in $dir" }

		my ( $id, $lib );
		if ( $no_db_for_rg ) { 
			( $id ) = split /\@/, basename $read1;
			$id =~ s/^se-//;
			if ( ! exists $libs{ $id } ) {
				die "no lib for $id ( from $read1 ) in $dir";
			}
			$lib = $libs{ $id };
		}
		else {

			my ($read_basename) = split /\@/, basename $read1;    # for filename matching from libs.info file
			$read_basename =~ s/^se-//;                           # remove pairing indicator
			$read_basename =~ s/[.]1$//;                          # remove read1/2 indicator

			if( $DEBUG ) { print STDERR  "check for libs{ $read_basename }\n"; }

			if (exists $libs{$read_basename}) {
				$lib = $libs{$read_basename};
			}
			
			die "no lib for $read1($read_basename) in $dir" unless ( defined $lib );
			
			#B/c $libs{$read_basename} must exist, otherwise program would have died already
			$id = $read_basename;

		}
		my $rg = " \'\@RG\\tID:$id\\tSM:$pinfo\\tPL:illumina\\tLB:$lib\'";

		push @seRGs, $rg;
		push @seR1s, $read1;

	}

	my $mixed        = @seRGs && @peRGs  ? 1 : 0;
	my $sge_args_opt = defined $sge_args ? "-sge-args \'$sge_args\'" : '';
	my @hold_jids    = ();

	my $run_bwa_path_to_tools = "-picard_path ".$picard_path." -samtools_path ".$samtools_path." -bwa_path ".$bwa_path;

	if (@seRGs) {
		my $out = $mixed ? "$pinfo.single-end" : $pinfo;
		my $cmd = "perl $reseq_bin_dir/run_bwa.pl ".$run_bwa_path_to_tools." -algo $bwa_algorithm $sge_args_opt -ref $ref_file -read1 " . join( " ", @seR1s ) . " -rg " . join( " ", @seRGs ) . " ";
		$cmd .= " $project_option -out $out -verbose -thread 4 -mem 8";
		print STDERR "$pinfo $cmd\n";
		my $rval = `$cmd`;
		push @hold_jids, split /\s+/, $rval;
	}
	if (@peRGs) {
		my $out = $mixed ? "$pinfo.paired-end" : $pinfo;
		my $cmd = "perl $reseq_bin_dir/run_bwa.pl ".$run_bwa_path_to_tools." -algo $bwa_algorithm $sge_args_opt -ref $ref_file -read1 " . join( " ", @peR1s ) . " -read2 " . join( " ", @peR2s ) . " -rg " . join( " ", @peRGs ) . " ";
		$cmd .= " $project_option -out $out -verbose -thread 4 -mem 8";
		print STDERR "$pinfo $cmd\n";
		my $rval = `$cmd`;
		push @hold_jids, split /\s+/, $rval;
	}
	if ($mixed) {    # mixture :(
		my @cmds = ();

		my $jids = join( ",", @hold_jids );
		my $cmd = "java -jar $picard_path/picard.jar MergeSamFiles I=$pinfo.single-end.bam I=$pinfo.paired-end.bam O=$pinfo.bam ";
		$cmd .= "AS=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=4000000 CREATE_INDEX=true ";
		$cmd .= "SO=coordinate CO=FixMateInformation CO=MarkDuplicates";

		push @cmds, $cmd;
		push @cmds, "mv $pinfo.bai $pinfo.bam.bai";
		push @cmds, "rm $pinfo.single-end.bam $pinfo.single-end.bam.bai";
		push @cmds, "rm $pinfo.paired-end.bam $pinfo.paired-end.bam.bai";
		push @cmds, "$samtools_path/samtools flagstat $pinfo.bam > $pinfo.bam.flagstat";

		write_array_cmds( \@cmds, 0, "run_remerger_${$}" );

		$cmd = "./run_remerger_${$}.0";
		print STDERR "$pinfo $cmd\n";
		my ($rjid) = split /\s+/, `$cmd`;
		print "$rjid\n";
	}
	else {
		print join(",", @hold_jids), "\n";
	}
	print "\n";
	chdir $start_dir or die "can't cd back to $start_dir $!";
	sleep 5;
}