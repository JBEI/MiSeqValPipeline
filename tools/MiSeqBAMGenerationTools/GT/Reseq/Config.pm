package GT::Reseq::Config;
use strict;
use warnings;
use Readonly;
use Exporter qw( import );

# config stuff specific to reseq pipelines

our @EXPORT = qw( 
  $SANGER_DATA_DIR
  $RESEQ_BIN_DIR
  $RESEQ_BIN_DIR_TEST
	get_users_email
	get_reseq_suffixes
	get_reseq_suffixes_alt
	get_file_suffixes_json
);

# TODO: remove Readonly
Readonly our $SANGER_DATA_DIR => '/home/nahituv/reseq/data';
Readonly our $RESEQ_BIN_DIR   => '/jgi/tools/groups/gentech/phoebetest/bin';
Readonly our $RESEQ_BIN_DIR_TEST => '/jgi/tools/groups/gentech/phoebetest/bin';
Readonly our $RESEQ_FILE_SUFFIX_URL => 'http://jel.jgi-psf.org/reseq/reseq_file_suffixes.json';

# bam pipeline standard suffixes and keys ( datatype in results.yaml )
my %reseq_suffixes = (
	'vcf',            '.bam.flt.anno.vcf.gz',
	'vcfidx',         '.bam.flt.anno.vcf.gz.tbi',
	'vcf-summary',    '.bam.flt.anno.vcf.snpEff_summary.html',
	'vcf-gene_stats', '.bam.flt.anno.vcf.snpEff_summary.genes.txt',
	'vcf-gff',        '.bam.flt.anno.vcf.gff.gz',
	'vcf-gffidx',     '.bam.flt.anno.vcf.gff.gz.tbi',
	'depth',          '.bam.depth.gz',
	'depthidx',       '.bam.depth.gz.tbi',
	'depthwig',       '.bam.depth.wig.gz',
	'depthU3',        '.bam.depth.under.3',
	'depthavg',       '.bam.depth_avg',
	'bam',            '.bam',
	'bamidx',         '.bam.bai',
	'breakdancer',    '.bd.max.filtered.xls',
	'nonstarters',    '.bam.nonstarters',
);

my %reseq_suffixes_alt = (  # common alternate suffixes
	'vcf',            '.bam.anno.vcf.gz',
	'vcfidx',         '.bam.anno.vcf.gz.tbi',
	'vcf-summary',    '.bam.anno.vcf.snpEff_summary.html',
	'vcf-gene_stats', '.bam.anno.vcf.snpEff_summary.genes.txt',
	'vcf-gff',        '.bam.anno.vcf.gff.gz',
	'vcf-gffidx',     '.bam.anno.vcf.gff.gz.tbi',
);

my %users_to_email = ( 
	'annau',          'alipzen@lbl.gov',
	'wendys',         'wsschackwitz@lbl.gov',
	'j_martin',       'j_martin@lbl.gov',
);

sub get_file_suffixes_json {
	use LWP::Simple;
	my $json = get( $RESEQ_FILE_SUFFIX_URL );
	if ( ! defined $json ) {
		die "ERROR: can't get json from $RESEQ_FILE_SUFFIX_URL $!";
	}
	return( $json );
}

sub get_users_email { 
	my $email;
	my $user = getlogin;

	if ( $user && exists $users_to_email{ $user } ) {
		$email = $users_to_email{ $user };
	}
	return( $email );
}

sub get_reseq_suffixes {
	return( %reseq_suffixes );
}
	
sub get_reseq_suffixes_alt {
	return( %reseq_suffixes_alt );
}

1;
