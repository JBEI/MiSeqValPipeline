package GT::Clutils;
use Modern::Perl;
use threads;
use version;
use Cwd;
use File::Basename;
use Log::Log4perl;

=head1 NAME 
 
GT::Clutils - clusterish utilities
  
=head1 VERSION
  
Version 0.01
  
=cut

our $VERSION = qv( 0.02 ); 

my $log = Log::Log4perl::initialized() ? Log::Log4perl->get_logger() : clog();

sub clog {
  Log::Log4perl->init(\ qq{
    log4perl.rootLogger=INFO, Screen
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.layout = SimpleLayout
  });
  return Log::Log4perl->get_logger('');
}

use base qw( Exporter );

our %EXPORT_TAGS = ( 'all' => [ qw( write_array_cmds log_and_die write_array_run_script
                                    run_array_local
                              ) ] );
 
our @EXPORT_OK   = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT      = @EXPORT_OK;
                   
=head1 SYNOPSIS
 
cluster cruft

$i = 1;
@jobs = ( "mv this that", "blat that", "touch done", "rm lockfile" ); 
write_array_cmds(\@jobs, $i++, 'run_blat' );
@jobs = ( "mv what thus", "blat thus", "touch done", "rm lockfile" ); 
write_array_cmds(\@jobs, $i++, 'run_blat' );

write_array_run_script( 'run_blat.sh' );

$i--;
print `qsub -t 1-$i ./run_blat.sh`;
 
=cut

my $default_run_file_name = 'run_file';

sub new {
  my ( $class, %args ) = @_;
  my $self = bless {%args}, $class;
  return $self;
}

=head1 handle_suspense

installs $SIG{TSTP} and $SIG{CONT} handler

 use GT::Clutils;
 $GT::Clutils->handle_suspense();
 $GT::Clutils->handle_suspense( $log, 
   { preload => 1, no_cont => 1, no_stop => 1 } );

takes code ref as first argument or will print to STDERR 

 to pass methods as code refs try something like

 my $log  = get_logger();
 my $warn = sub { $log->warn( @_ ) }; 

 GT::Clutils->handle_suspense( $warn );

 or even 

 my $warn = sub { $log->warn( @_ ) }; 
                  shut_db if ( $_[0] =~ /^TSTP/ );
                  open_db if ( $_[0] =~ /^CONT/ );
            };

=cut

sub handle_suspense {
  my ( $self, $log, $opts ) = @_;

  if ( ! defined $log ) {
    $log = \&print_stderr;
  }

  $SIG{ TSTP } = sub { &$log( 'TSTP received' ); } unless( $opts->{ no_stop } );
  $SIG{ CONT } = sub { &$log( 'CONT received' ); } unless( $opts->{ no_cont } ); 

  if ( $opts->{ preload } ) {
    &$log( 'Handler installed' );
  }

  sub print_stderr {
    print STDERR @_, "\n";
  }
}
 
=head2 write_array_run_script
 
write_array_run_script( $run_file_name )
 
 run_file_name default name is 'run_file.sh'
 
 writes the script for sge array job submission
 ( the one that that calls the others )
 
=cut 

 sub write_array_run_script { 
  my ( $run_file_name, @task_ids ) = @_; 
                   
  $run_file_name ||= "$default_run_file_name.sh";
 
  my $run_file_base = basename($run_file_name, qw( .sh .pl .csh .bash .tcsh ) );
 
  open( my $fh, '>', $run_file_name ) or log_and_die(" can't open $run_file_name for writing $!" );
     
  print $fh <<EORF; 
#!/bin/bash 
   
umask 2
EORF
   for (my $i = 0; $i<scalar @task_ids; $i++) {
  print $fh <<EORF; 
  ./$run_file_base.$task_ids[$i]
EORF
  }
  close( $fh ) or log_and_die("can't close $run_file_name after writing $!");
  chmod 0775, $run_file_name or log_and_die("can't chmod 0755 $run_file_name after closing $!");
  return( $run_file_name );
} 

=head2 write_array_cmds {
 
write_array_cmds( $cmd, $job_index, $run_file_name );
 
$cmd can be scalar command or ref to array of commands
 
writes script files that contain the actual commands executed
expects job_index, optionally accepts base name for script
 
=cut
 
sub write_array_cmds {
  my ( $cmd, $job_index, $run_file_name, @bad_args ) = @_;
  if ( @bad_args ) {
    die "write_array_cmds called with too many options, check that array of command is passed by reference";
  }
 
  # $log->info("self: gone\tcmd: $cmd\tjob_index: $job_index\trun_file_name: $run_file_name\n");
 
  if ( ! defined $run_file_name ) { $run_file_name = $default_run_file_name; }
  if ( ! defined $job_index     ) { $job_index     = 0;                      }
 
  my $umask = umask;
  umask 2;
 
  open( my $fh, '>', "$run_file_name.$job_index" ) or die "can't open $run_file_name.$job_index for writing $!";
 
  umask $umask;
 
  my $cmds = ref $cmd ? join("\nexit_on_error \$?\n", @$cmd ) : $cmd;

  # -hold_jid* only wants to hold things that exit with code 100
 
  print $fh <<EORF;
#!/bin/sh
umask 2
exit_on_error()
{
  if [ "\${1}" -ne "0" ]; then
    echo exit on error code \${1}
    exit 100
  fi
}

$cmds
if [ \$? -eq 0 ]; then
  if [ ! -e scripts ]; then
   mkdir scripts
  fi
  mv \$0 scripts
else
  echo exit on error code \$?
  exit 100
fi
 
EORF
  close( $fh ) or die "can't close $run_file_name.$job_index after writing $!";
  chmod 0775, "$run_file_name.$job_index" or die "couldn't set $run_file_name.$job_index executable $!";
}
 
=head2 run_array_local
 
runs 'run_file.\d+$' in series
 
=cut
 
sub run_array_local {
  my ( $run_file_name, $max_threads ) = @_;
 
  my ( $run_file_base, $task_range, $tasks, $scripts ) = get_scripts( $run_file_name);
 
  my $MAX_THREADS = defined $max_threads          ? $max_threads
                  : exists $ENV{MAX_THREADS}      ? $ENV{MAX_THREADS}
                  : exists $ENV{OMP_NUM_THREADS}  ? $ENV{OMP_NUM_THREADS}
                  : 8;
 
  my @threads = ();
 
 
  for my $task ( @$tasks ) {
    while ( threads->list(threads::running) >= $MAX_THREADS ) {
      sleep 2;
    }
    # $log->info("executing ./$run_file_base.$task");
    push @threads, threads->create(\&run_cmd_or_die, "./$run_file_base.$task");
  }
  for my $thread( @threads ) {
    $thread->join();
  }
  return( 1 );  # didn't die so returns true.
}

 
=head1 run_cmd_or_die
 
system() wrapper
 
=cut
 
sub run_cmd_or_die {
  my ( $cmd ) = @_;
  my $from = ( caller( 1 ) ) [ 3 ];
  $log->info("($from) system( $cmd )");
  system( $cmd ) && log_and_die( "ERROR: failed executing $cmd" );
}

=head2 log_and_die {
 
logs to fatal and dies
 
=cut
 
sub log_and_die {
  my ( $message ) = @_;
 
  $message ||= 'no info supplied';
  my $from = ( caller( 1 ) ) [ 3 ];
  $log->fatal(" ( $from ) $message");
  croak($message);
}
 
 
=head1 get_scripts     
 
 takes optional script base name ( default 'run_script' )
 returns
   $run_file_base, $task_range, \@tasks, \@scripts
 
=cut
 
sub get_scripts {
  my ( $run_file_name ) = @_;
 
  $run_file_name ||= "$default_run_file_name.sh";

  my $run_file_base = basename( $run_file_name, qw( .sh .csh .bash .tcsh .pl ) );
 
  my @tasks = ();
 
  my @scripts = grep { /^$run_file_base.(\d+)$/ && do { push @tasks, $1 } } glob( "$run_file_base.[0-9]*" );
 
  @tasks = sort { $a <=> $b } @tasks;
 
  if ( ! defined  $tasks[ 0 ] || ! defined $tasks[ -1 ] ) {
    $log->warn( "tasks undeffed, run_file_name=$run_file_name\trun_file_bae=$run_file_base\tscripts=@scripts");
  }
  my $task_range = $tasks[ 0 ] . '-' . $tasks[ -1 ];
 
  return( $run_file_base, $task_range, \@tasks, \@scripts );
}

