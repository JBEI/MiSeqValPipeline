#!/usr/bin/env perl
use Modern::Perl;
use Cwd;
use GT::Clutils;
use GT::Utils;
use File::Basename;
use Getopt::Long;
use GT::Reseq::Utils;

umask 2;

#Location of fastq_slice, guess_qual_format, and fastq_single_or_paired command-line commands

my ( $base_dir, $fastq_dir, $fqsource_dir, $sge_project_opt, $config_file, $mainlib_dir, $reseq_bin_dir);

$base_dir     = 'bwa_dir';
$fastq_dir    = 'fastq_dir';    # in base_dir
$fqsource_dir = 'fastq_dir';    # in dir

my $help;

GetOptions(
	'basedir=s'         => \$base_dir,
	'sourcedir=s'       => \$fqsource_dir,
	'Project|project=s' => \$sge_project_opt,
	'config=s'          => \$config_file,
	'mainlibdir=s'		=> \$mainlib_dir, #Absolute path to main Library directory
	'reseqbindir=s'		=> \$reseq_bin_dir,
	'help'              => \$help,
);

my @dir_list = @ARGV;


if( $config_file ) {
	@dir_list = ru_add_dirs_from_config_file( $config_file, @dir_list );
}

if ( $help || !defined $base_dir || ! @dir_list ) {
	my $bin = basename $0;
	print STDERR "\n$bin [options] dirs*\n\n";
	print STDERR " -base   [bwa_dir]    alignment directory, sliced files will be base/fastq_dir\n";
	print STDERR " -source [fastq_dir]  subdirectory with source fastqs \n";
	print STDERR " -P       STRING      sge accounting project name for qsub\n";
	print STDERR " -config  FILE        config file ( libr.info format )\n";
	print STDERR " -mainlibdir FILE     absolute path to main Library directory\n";
	print STDERR "\n";

	exit;
}

$sge_project_opt = $sge_project_opt ? " -P $sge_project_opt " : '';

my $slots  = 2;    # default single
my $i      = 1;
my $script = "run_slice_${$}";

#print "dir_list is @dir_list\n";

@dir_list     = grep { -d $_ } @dir_list;
my %dreads    = ();
my %bad_links = ();

for my $dir ( @dir_list ) {
	# if your fastq files aren't gzipped, you are a bad person
	my @reads = glob("$mainlib_dir/$dir/$fqsource_dir/*fastq.gz $mainlib_dir/$dir/$fqsource_dir/*.fq.gz");
	if ( ! @reads ) { 
		warn "no gzipped fastq files found in $mainlib_dir/$dir/$fqsource_dir\nEnsure your fastq files are gzipped\n";
		next;
	}
	else {
		for my $read ( @reads ) {
			if ( ! -e $read ) {
				push @{ $bad_links{ $dir } }, $read;
			}
		}
		$dreads{ $dir } = \@reads;
	}
}
if ( %bad_links ) {
	print STDERR "\nERROR: some linked fastqs are missing ( restoring from tape perhaps? )\n\n";

	for my $dir ( sort keys %bad_links ) {
		for my $file ( @{ $bad_links{ $dir } } ) {
			print STDERR "missing: $file\n";
		}
	}
	print "\n";
	exit ;
}

for my $dir ( sort keys %dreads ) {
	
	if ( !-e "$mainlib_dir/$dir/$base_dir" ) {
		mkdir "$mainlib_dir/$dir/$base_dir" or die $!;
	}
	if ( !-e "$mainlib_dir/$dir/$base_dir/bam_dir" ) {
		mkdir "$mainlib_dir/$dir/$base_dir/bam_dir" or die $!;
	}
	if ( !-e "$mainlib_dir/$dir/$base_dir/$fastq_dir" ) {
		mkdir "$mainlib_dir/$dir/$base_dir/$fastq_dir" or die $!;
	}

	for my $read ( @{ $dreads{ $dir } }) {
		my $name  = basename $read;
		my $bname = $name;
		$bname =~ s/[.]fastq[.]gz$//;
		$bname =~ s/[.]fq[.]gz$//;

		my $type; # if undef it's guessed

		# if any paired are being sliced set pe to 4, to be kind.

		# TODO: to be efficient, should qsub separately single and paired

		if ( $bname =~ /^pe-/ ) {
			$type  = 'paired';
		}
		elsif ( $bname =~ /^se-/ ) {
			$type = 'single';
		}

		# check for empty fastq
		zopen( my $fh, $read );
		my $has_a_read = 0;
		while (<$fh>) {
			$has_a_read++;
			last if ( $has_a_read >= 3 );
		}
		close($fh);
		if ( !$has_a_read ) {
			warn "$dir has empty fastq $read, skipping";
			next;
		}
		my $guess = `$reseq_bin_dir/guess_qual_format $read`;
		chomp $guess;
		my $illoption = $guess eq 'sanger' ? '' : ' -ill2std ';
		my $cmd;
		if ( defined $type && $type eq 'paired' ) {
			print STDERR "type is paired\n";
			my $prefix = "$mainlib_dir/$dir/$base_dir/$fastq_dir/$bname";
			$cmd = "perl $reseq_bin_dir/fastq_slice.pl -verbose -gzip -num 8000000 $illoption -suffix fq -file $read -prefix $prefix";
		}
		elsif ( defined $type && $type eq 'single' ) {
			print STDERR "type is single end\n";
			my $prefix = "$mainlib_dir/$dir/$base_dir/$fastq_dir/$bname";
			$cmd = "perl $reseq_bin_dir/fastq_slice.pl -verbose -gzip -num 8000000 $illoption -suffix fq -file $read -prefix $prefix";
		}
		else {  # DEFAULT CASE, above could probably be removed as it's remnant from the maq->bwa switchover
			if ( !defined $type ) {
				($type) = split /\s+/, run_bt_or_die( "perl $reseq_bin_dir/fastq_single_or_paired.pl $read" );
			}
			if ( $type eq 'paired' ) {
				my $prefix = "$mainlib_dir/$dir/$base_dir/$fastq_dir/pe-$bname";
				# $prefix .= "perl $barcode" if ( $barcode );
				$cmd = "perl $reseq_bin_dir/fastq_slice.pl -verbose -split -gzip -num 8000000 $illoption -suffix fq -file $read -prefix $prefix";
			}
			elsif ( $type eq 'single' ) {
				my $prefix = "$mainlib_dir/$dir/$base_dir/$fastq_dir/se-$bname";
				# $prefix .= "perl $barcode" if ( $barcode );
				$cmd = "perl $reseq_bin_dir/fastq_slice.pl -verbose -gzip -num 8000000 $illoption -suffix fq -file $read -prefix $prefix";
			}
			else {
				warn "WARNING: skipping read $read don't understand what type it is $!";
				next;
			}
		}
		my $resultOfcmd = `$cmd`;
		if ( ${^CHILD_ERROR_NATIVE} ) {
			die "error running $cmd ( error = ${^CHILD_ERROR_NATIVE} )";
		}
	}
}
