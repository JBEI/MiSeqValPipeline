package GT::Reseq::Utils;
use strict;
use warnings;
use version;
use parent qw( Exporter );

=head1 NAME

GT::Reseq::Utils - reseqish utilities

=head1 VERSION

Version 0.02

=cut

our $VERSION = qv( 0.02 );

our @EXPORT  = qw( ru_add_dirs_from_config_file );

=head1 ru_add_dirs_from_config_file()

 my @dirs = ru_add_dirs_from_config_file( 'lib.info', grep { -d } @ARGV );
 or
 my @dirs = ru_add_dirs_from_config_file( 'lib.info' );

=cut

sub ru_add_dirs_from_config_file {
	my ( $config_file, @dirs ) = @_;
	my $order    = 1;
	my %dir_list = ();
	for my $dir ( @dirs ) {
		$dir_list{ $dir } = $order++;
	}
	if( -e $config_file ) {
		open( my $fh, $config_file ) or die $!;
		while ( <$fh> ) {
			next if ( /^#/ );
			next unless ( /\S/ );
			my ( $pmo_info ) = split;
			$dir_list{$pmo_info} = $order++ if ( -d $pmo_info );
		}
		close( $fh ) or die $!;
	}
	return( sort { $dir_list{ $a } <=> $dir_list{ $b } } keys %dir_list );
}
