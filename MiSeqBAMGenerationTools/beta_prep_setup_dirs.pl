#!/usr/bin/env perl
use Modern::Perl;
use Getopt::Long;
use File::Basename;
use Cwd qw/realpath/;
use Data::Alias;
use GT::Config;

use YAML qw/LoadFile DumpFile/;

umask 2;

my ( $bwa_dir, $ref, $rnaseq_style, $config_file, $help );

$bwa_dir = 'bwa_dir';

GetOptions(
	'bwa_dir|base=s'   => \$bwa_dir,
	'ref_fasta=s' => \$ref,
	'rna'         => \$rnaseq_style,
	'help'        => \$help,
	'config=s'    => \$config_file, # for consistency
);

my $info_file = defined $config_file ? $config_file : shift;

if ( $help || !defined $ref || !-f $ref || !defined $info_file || !-e $info_file ) {
	my $bin = basename $0;
	print STDERR <<EOH;
Usage: $bin -ref /path/to/ref.fa config.info [ -base bwa_15x ]
  -rna_fasta        switch to rna style naming 
                 a) rewrites config file if given config file
                    saves original to config.bak.#
  -config    FILE   config file ( aka libr.info )
               created by beta_prep_lib_info
format is tab delimited
proposalId_projectId  LIBR  fastq.name.fq.gz  genus  species  strain
EOH
	exit;
}

$ref = realpath $ref unless ( $ref =~ /^\// );

open( my $ifh, $info_file ) or die $!;
my @info = <$ifh>;
close( $ifh ) || die $!;

if ( $rnaseq_style ) {
	# get uniq name for backup
	my $suffix = '';
	for my $i ( 1 .. 99 ) {
		if ( ! -e "$info_file.bak.$i" ) { 
			$suffix = "bak.$i";
			last;
		}
	}
    $suffix ||= "bak.99";

    # should write temp first then rename 
	rename( $info_file, "$info_file.$suffix" ) or die $!;

	my @rna_info = ();
	open( my $rna_fh, '>', $info_file ) or die $!;
	for my $line ( @info ) {
		my ( $id, $lib, @stuff ) = split /\t/, $line;
		$id .= '_' . $lib;
		$line = join("\t", $id, $lib, @stuff );
		print $rna_fh $line;
	}
	close( $rna_fh ) or die $!;
}

my $RETRIEVING_FROM_TAPE = 0;

my $SDM_ON_TAPE = 2;
my $SDM_FILE_READY = 0;
my $SDM_ERROR      = 1;

for (@info) {
	next if ( /^#/ );
	next unless ( /\S/ ); # skip blanks
	chomp;
	my ( $pmo, $lib, $fastq, $organism ) = split /\s+/, $_, 4;
	
	my ( $fq_path, $fp_fq );

	#Get path to FASTQ.GZ files
	$fp_fq = $fastq;
	$fastq = basename $fastq;

	print "$pmo\t$lib\t$fp_fq\t$organism\n";

	if ( !-e $pmo ) {
		makedir($pmo);
	}
	if ( !-e "$pmo/$bwa_dir" ) {
		makedir("$pmo/$bwa_dir");
	}
	if ( !-e "$pmo/fastq_dir" ) {
		makedir("$pmo/fastq_dir");
	}

	if ( ! -l "$pmo/fastq_dir/$fastq" ) {
		symlink( $fp_fq, "$pmo/fastq_dir/$fastq" ) or die "can't symlink( $fp_fq, $pmo/fastq_dir/$fastq ) $!";
	}

	# create or update config.yml

	my $config_file = "$pmo/$bwa_dir/config.yml";
	my $config;

	if ( -e $config_file ) {
		$config = LoadFile( $config_file );
	}

	$config->{ref_fasta} = $ref;
	if ( -e "$ref.fai" ) {
		$config->{ref_faidx} = "$ref.fai";
	}

	my $rg_id = $fastq;
	$rg_id =~ s/.fastq.gz$//;
	$rg_id =~ s/.fq.gz$//;

	alias my %config_fqs = %{ $config->{fqs} };

	$config_fqs{ $fastq }->{src_path}       = $fp_fq;
	$config_fqs{ $fastq }->{library}        = $lib;
	$config_fqs{ $fastq }->{organism}       = $organism;
	$config_fqs{ $fastq }->{rg_sample_name} = $pmo;
	$config_fqs{ $fastq }->{rg_id}          = $rg_id;

	DumpFile( $config_file, $config )
}

sub makedir {
	my ($dir) = @_;
	mkdir $dir or die "can't mkdir $dir $!";
	chmod 02775, $dir or die "can't chmod 02775 $dir $!";
}