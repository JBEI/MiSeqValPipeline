package GT::BWA::Cleanup;
use Modern::Perl;
use File::Basename;
use Cwd;
use version;

=head1 NAME 
 
GT::BWA::Cleanup - 
  
=head1 VERSION
  
Version 0.01

=cut

our $VERSION = qv( 0.02 ); 

use base qw( Exporter );

our %EXPORT_TAGS = ( 'all' => [ qw( delete_bwa_inputs is_bam_corrupt ) ] );
 
our @EXPORT_OK   = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT      = @EXPORT_OK;
                   
=head1 SYNOPSIS
 
=head1 delete_bwa_inputs

=cut

  # - clean - files not needed
  # clear gz from base/fastq_dir
  # clear sai from base/sai_dir
  # clear split bam from base/bam_dir

sub delete_bwa_inputs {
	my ( $bwa_dir ) = @_;
	
	my %file_globs = ( fastq_dir => '*.fq.gz', sai_dir => '*.sai', bam_dir => '*.bam' );

	for my $sub_dir ( keys %file_globs ) {
		my $glob = "$bwa_dir/$sub_dir/$file_globs{ $sub_dir }";
		print STDERR "globbing $glob\n";
		my @files = glob( $glob );
		for my $file ( @files ) {
			print "unlink $file\n";
		}
	}
}

# should be in GT::BWA or GT::BAM
sub is_bam_corrupt {
	my ( $bam_file ) = @_;
	my $cmd = "samtools view -H $bam_file | head -0";
	my $is_corrupt = 0;
	system($cmd) && do {$is_corrupt = 1;};
	return( $is_corrupt );
} 



1;
