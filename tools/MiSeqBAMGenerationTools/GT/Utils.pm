package GT::Utils;
use warnings;
use strict;
use version;
use Data::Alias;
use Cwd qw( realpath );

use base qw( Exporter );

 
our @EXPORT      = qw( sort_by_last_num zopen bopen bfqopen run_bt_or_die );

=head1 SYNOPSIS

utilitarian functions

=cut

sub new {
  my ( $class, %args ) = @_;
  my $self = bless {%args}, $class;
  return $self;
}

=head1 sort_by_last_num 

@a = sort_by_last_num( @a );

for numeric sorting of lists like
qw( run_thing.1 run_thing.10 run_thing.2 )
-> run_thing.1, run_thing.2, run_thing.10 
qw( read.1@1.fq.gz read.2@1.fq.gz read.1@11000.fq.gz read.2@11000.fq.gz read.1@3000.fq.gz read.2@3000.fq.gz )
-> read.1@1.fq.gz read.2@1.fq.gz read.1@3000.fq.gz read.2@3000.fq.gz ead.1@11000.fq.gz read.2@11000.fq.gz 
  - note 'read.1' and 'read.2' NOT involved in the sorting 

=cut

sub sort_by_last_num {
 my @idx = ();
 for my $idx ( @_ ) {
   my $val;
   while ( $idx =~ /(\d+)/g ){ $val = $1 }
   push @idx, $val;
 }
 return( @_[ sort { $idx[ $a ] <=> $idx[ $b ] } 0 .. $#_ ] ); # ugly but efficient
}

=head1 zopen

zopen( my $fh, 'my.gz' ) or die $!;
while( <$fh> ) { print }
zopen( $fh, 'my.bz2' ) or die $!;
while( <$fh> ) { print }
close( $fh ) or die $!;

=cut

sub zopen{
  alias my ( $fh, $file ) = @_;

  my $open = $file =~ /[.]gz$/  ? "gzip  -cd $file |"
           : $file =~ /[.]bz2$/ ? "bzip2 -cd $file |"
           : $file;

  return( open( $fh, $open ) );
}

=cut

=head1 bopen

bopen( my $fh, 'my.bam' ) or die $!;
while( <$fh> ) { print }
bopen( $fh, 'my.sam' ) or die $!;
while( <$fh> ) { print }

=cut

sub bopen{
  alias my ( $fh, $file ) = @_;
  #
  # TODO: ? add uncompressed bam vs compressed bam 

  my $open = "samtools view $file |";

  return( open( $fh, $open ) );
}

=cut

=head1 bfqopen

bfqopen( my $fh, 'my.bfq' ) or die $!;
while( <$fh> ) { print }

=cut

sub bfqopen{
  alias my ( $fh, $file ) = @_;

  my $open = "maq bfq2fastq $file - |";

  return( open( $fh, $open ) );
}

=head1 run_bt_or_die

executes `$arg1` and returns the result
dies if $^CHILD_ERROR_NATIVE is set

my $result = run_bt_or_die('ls');

=cut

sub run_bt_or_die {
  my ( $cmd ) = @_;
 
  my $result = `$cmd`;
 
  if ( ${^CHILD_ERROR_NATIVE} ) {
    die "error running $cmd ( error = ${^CHILD_ERROR_NATIVE} )";
  }
 
  return $result;
}

1;
