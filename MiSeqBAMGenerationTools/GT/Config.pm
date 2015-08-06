package GT::Config;
use Modern::Perl;
use YAML::Syck;
use File::NFSLock;
use Fcntl qw ( LOCK_EX LOCK_SH LOCK_NB );

sub new {
  my ( $class, %args ) = @_;
  my $self = bless { 
    config_file => 'config.yml',
    ref_fasta   => '',
    %args,
  }, $class;

  return( $self );
}

sub load {
  my ( $self, $file ) = @_;

  if ( $file ) {
    $self->{config_file} = $file;
  }

  if ( ! -e $self->{config_file} ) {
    return undef;
  } 
  # locks fxored on gpfs 20140821 6:30pm
	my $lock = File::NFSLock->new($self->{config_file},LOCK_SH|LOCK_NB,10,30*60);
 
  my $config = LoadFile( $self->{config_file} );

  # locks fxored on gpfs 20140821 6:30pm
  $lock->unlock;

  while ( my ( $key, $val ) = each %$config ) {
    $self->{$key} = $val;
  }
  return;
}

sub ref_fasta {
  my ( $self, $arg ) = @_;
  $self->{ref_fasta} = $arg if defined $arg;
  return( $self->{ref_fasta} );
}

1;