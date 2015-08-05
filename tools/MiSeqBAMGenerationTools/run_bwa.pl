#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use File::Basename;
use Getopt::Long;
use Data::Alias;
use POSIX qw/ceil/;
use Cwd;
use GT::Clutils;
use GT::Utils;

# TODO: handle mix of single and paired end data in single call
# --- bs: -read1 x.fq.gz -read2 y.fq.gz -rg '@RG\tID:pe' -se z.fq.gz -rgse '@RG\tID:se'
#         if ( se )?

# TODO: make identifiable job name, last 2 digits of project_id and the job id?

$|=1;

# 
# avoid passing -V as it errors out nodes, 
# 


my ( $ref_fasta, @read1, @read2, @read_group_string, $reuse_rg, $out_bam, $help, $mem_in_g, $threads, 
	$picard_path, $bwa_path, $samtools_path, 
     $args_samse, $args_sampe, $args_aln, $args_mem, $args_sge, $sleep_time, $sge_project,
	   $algorithm, $verbose );

$mem_in_g       = 5;
$threads        = 4;     # default to 4 as cluster is annoying to schedule 8 thread anything on
$args_sampe     = '-a 800';
$args_mem       = '-M';
$args_aln       = '';
$args_sge       = '';
$sleep_time     = 10;
$algorithm      = 'backtrack';

#Location of Picard, BWA, and samtools commands
$picard_path = "";
$bwa_path = "";
$samtools_path = "";

GetOptions( 
  'ref=s'       => \$ref_fasta,
  'read1=s{,}'  => \@read1,
  'read2=s{,}'  => \@read2,
  'out=s'       => \$out_bam,
  'rg=s{,}'     => \@read_group_string,
  'reuse-rg'    => \$reuse_rg,
  'mem=i'       => \$mem_in_g,
  'help'        => \$help,
  'threads=i'   => \$threads,
  'picard_path=s'   => \$picard_path,
  'bwa_path=s'      => \$bwa_path,
  'samtools_path=s' => \$samtools_path,
  'aln-args=s'  => \$args_aln,
  'sge-args=s'  => \$args_sge,
  'sleep=i'     => \$sleep_time,
  'P=s'         => \$sge_project,
  'verbose'     => \$verbose,
  'samse-args'  => \$args_samse,
  'sampe-args'  => \$args_sampe,
  'mem-args'    => \$args_mem,
  'algorithm=s' => \$algorithm,

) or die;

sub usage {
  my $bin = basename $0;
  print "\nusage: $bin [options] -ref ref.fasta -read1 read_1.fastq [ -read2 read_2.fastq ] -out prefix [ -mem 10 ][ -t 8 ] \n";
  print "\n";
  print " -mem should be 1/3rd max mem used in sort step\n";
  print "\n";
  print " -ref           STRING   reference in fasta format, should have 'samtools index ref' and 'bwa index ref' run on it\n";
  print " -rg            STRING   read group string enclosed in quotes ( e.g. \'\@RG\\tID:id\\tSM:sample\\tLB:lib\' )\n";
  print " -reuse-rg               use same read group string for all read files\n";
  print " -nocluster              run locally\n";
  print " -aln-args      STRING   args for bwa aln ( e.g. -aln-args \'-M 4 -O 14\' )\n";
  print " -sampe-args    STRING   args for bwa sampe ( e.g. -aln-args \'-a 1200\' )\n";
	print "                         default: $args_sampe\n";
  print " -samse-args    STRING   args for bwa samse ( e.g. -samse \'-n 200000\' )\n";
	print " -nomodules              demodulate and use the source luke\n";
	print " -algorithm     STRING   mem or backtrack or bwasw ( backtrack is ( aln + ( samse or sampe ) ) ) and is default\n";
  print " -sge-args      STRING   args for sge e.g. -sge-args \"-l high.c\" ***USE THEM QUOTES***\n";
  print " -P             STRING   sge project ( e.g. jgi.p )\n";
  print " -sleep         INT      ( for testing, seconds to sleep between file stats )\n";
  print "\n";
	exit;
}

if ( $help || ! defined $ref_fasta || ! -e $ref_fasta  
  || ( ! @read1 || ! defined $out_bam ) 
  || ( @read2 &&  $#read2 != $#read1 ) ){
  usage();
}

for ( @read1, @read2 ) {
  if ( ! -e ) { 
    usage();
  }
}

my $pic_factor  = 0.9;
my $MEM_P_VALUE = int( $mem_in_g * $pic_factor );             # picard memory - fixmate is an offender, check if scales with readsize
my $MAX_RECORDS = ( $MEM_P_VALUE * 500000 ) - ( 250000 * $mem_in_g );

# not used -> my $mem_s_value = $mem_in_g * 2.5 * 1000000000; # samtools sort memory, it's quite off base

my ($cmd, $jid, %jids, %jids_aln);
  
my $job_prefix = 'bwa';

$sge_project = $sge_project ? " -P $sge_project " : '';


if ( ! -e "$ref_fasta.pac" ) {
	create_bwa_index( $ref_fasta, $job_prefix );
}

if ( ! -e "$ref_fasta.fai" ) {
  die "\nrun samtools faidx $ref_fasta first\n";
	# module load samtools ...
  $cmd = "$samtools_path/samtools faidx $ref_fasta";
  system( $cmd ) && die "ERROR: command failed ( $cmd ) $!";
}

# align read1

my $hold_option = '';

if ( keys %jids ) { 
  $hold_option = "-hold_jid " . join(",",sort keys %jids);
}

if ( $reuse_rg ) {
  @read_group_string = map { $read_group_string[ 0 ] } @read1;
}

my $is_paired   = @read2 ? 1 : 0;

  $algorithm eq 'backtrack' ? samXe() 
: $algorithm eq 'mem'       ? mem()
: $algorithm eq 'bwasw'     ? bwasw()
: die "unknown algorithm $algorithm";

######################
# # #    done    # # # 
######################


sub mem {
	my $mem_i      = 1;
	my $mem_script = "run_mem_${$}";
	my @bams       = ();
	my @uniqueIDs   = ();

	# generate alignment scripts
	for ( my $i = 0; $i <= $#read1; $i++ ) {
		alias my $r1 = $read1[ $i ];
		alias my $r2 = $read2[ $i ];
		alias my $rg = $read_group_string[ $i ];
		push @bams, setup_mem_reads( $r1, $r2, $rg, $args_mem, \$mem_i, $mem_script );
		push @uniqueIDs, ($i+1);
	}
	write_array_run_script( "$mem_script.sh", @uniqueIDs );
	$mem_i--;
	if ( ! $mem_i ) {
		warn "WARNING: no alignments to run!\n";
		exit;
	}
	my $ram_c_mem;
	my $h_vmem;
	$h_vmem = 1500;
	$ram_c_mem = ceil( $h_vmem / $threads );
	$h_vmem    .= 'M';
	$ram_c_mem .= 'M';
	
	$cmd = "./$mem_script.sh";
	# post process bams
	
	print STDERR "\nCMD: $cmd\n" if $verbose;

	my $mem_jid = run_or_die($cmd);
	chomp $mem_jid;
	( $mem_jid ) = split /[.]/, $mem_jid;

	my $last_jid = merge_and_index( $out_bam, $mem_jid, @bams );	

	print "Alignments finished";
}

sub merge_and_index {
	my ( $out_bam, $hold_jid, @bams ) = @_;

	my $cmd;

	if ( $#bams ) {
 		my @cmds = ();
		if ( -e "$out_bam.bai" || -l "$out_bam.bai" ) {
			push @cmds, "rm $out_bam.bai";
		}
		if ( -e "$out_bam.bam.bai" || -l "$out_bam.bam.bai") {
			push @cmds, "rm $out_bam.bam.bai";
		}
 		push @cmds, "java -jar $picard_path/picard.jar MergeSamFiles I=".join(" I=", @bams )." O=$out_bam.bam VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=$MAX_RECORDS AS=true SO=coordinate CREATE_INDEX=true";
 		push @cmds, "mv $out_bam.bai $out_bam.bam.bai";
		push @cmds, "$samtools_path/samtools flagstat $out_bam.bam > $out_bam.bam.flagstat";
		push @cmds, "rm @bams";

 		my $script_name = "run_merger_${$}";
 		write_array_cmds( \@cmds, 0, $script_name);

 		$cmd = "./$script_name.0";
	}
	else { 
		my @cmds = ();
		push @cmds, "mv $bams[0] $out_bam.bam"; 
		if ( -e "$out_bam.bai" || -l "$out_bam.bai" ) {
			push @cmds, "rm $out_bam.bai";
		}
		if ( -e "$out_bam.bam.bai" || -l "$out_bam.bam.bai") {
			push @cmds, "rm $out_bam.bam.bai";
		}
		push @cmds, "$samtools_path/samtools index $out_bam.bam"; 
		push @cmds, "$samtools_path/samtools flagstat $out_bam.bam > $out_bam.bam.flagstat";
	
 		my $script_name = $is_paired ? "run_pindex_${$}" : "run_sindex_${$}";  
 		write_array_cmds( \@cmds, 0, "$script_name" );
 		$cmd = "./$script_name.0";
	}

	print STDERR "\nCMD: $cmd\n" if $verbose;

	my $jid = run_or_die( $cmd );
	chomp $jid;
	( $jid ) = split /[.]/, $jid;

	return( $jid );
}

sub samXe { 
	my $args_samXe  = $is_paired ? $args_sampe : $args_samse; 

	# filename prefixes for the task scripts
	my $aln_script_pe   = "run_aln-pe_${$}";
	my $aln_script_se   = "run_aln-se_${$}";
	my $samXe_script_pe = "run_sampe_${$}";
	my $samXe_script_se = "run_samse_${$}";
	my $aln_script      = $is_paired ? $aln_script_pe   : $aln_script_se;
	my $samXe_script    = $is_paired ? $samXe_script_pe : $samXe_script_se;
		
	my @bams         = ();

	# indexes to count # of tasks
	my $samXe_i      = 1;
	my $aln_i        = 1;
	my @uniqueIDs    = ();

	for ( my $i = 0;  $i <= $#read1; $i++ ) {
 		alias my $r1 = $read1[ $i ];
 		alias my $r2 = $read2[ $i ];             # gets checked later in aln_backtrack_reads
 		alias my $rg = $read_group_string[ $i ];
		push @bams, setup_aln_backtrack_reads( $r1, $r2, $rg, $args_samXe, \$samXe_i, \$aln_i, $aln_script, $samXe_script  );
		push @uniqueIDs, ($i+1);
	}
	write_array_run_script( "$aln_script.sh", @uniqueIDs );
	write_array_run_script( "$samXe_script.sh", @uniqueIDs );

	$aln_i--;

	if ( ! $aln_i ) {
 	 warn "WARNING: no alignments to run!\n"; exit;
	}
	my $ram_c_mem;
	my $h_vmem   ;
	$h_vmem = 1500;
	$ram_c_mem = ceil( $h_vmem / $threads );
	$h_vmem    .= 'M';
	$ram_c_mem .= 'M';
	$cmd = "./$aln_script.sh";

	print STDERR "\nCMD: $cmd\n" if $verbose;

	my $aln_jid = run_or_die($cmd);
	chomp $aln_jid;

	( $aln_jid ) = split /[.]/, $aln_jid;

	# TODO: alter ram.c setting after 634346 completes.  either 1/2 it or add h_vmem=2*4mem_q_value
	#  *** request mem of ( sampe + view + SortSam -- alter SortSam request, it's well behaved? or how does max records scale with readlength?)

	my $sam_slots = 2;
	$ram_c_mem = ceil($mem_in_g/$sam_slots) . 'G';
	$h_vmem    = $mem_in_g . 'G';

	$cmd = "./$samXe_script.sh";

	print STDERR "\nCMD: $cmd\n" if $verbose;

	my $sX_jid = run_or_die( $cmd );;
	chomp $sX_jid;
	( $sX_jid ) = split /[.]/, $sX_jid;

	my $last_jid = merge_and_index ( $out_bam, $sX_jid, @bams );

	# TODO: add cmd -hold_jid $last_jid rm @bams;  rm sai files
	#      or add this all to a script for the last_jid

	print "$aln_jid\t$sX_jid\t$last_jid\n";

}

# command strings for picard the verbose *note file handle limit

sub picard_fmi_string{ 
	my ( $mem, $i, $o, $mr ) = @_;
	return( "java -jar $picard_path/picard.jar FixMateInformation I=$i O=$o MAX_RECORDS_IN_RAM=$mr VALIDATION_STRINGENCY=SILENT SO=coordinate" );
}

sub picard_md_string{
	my ( $mem, $i, $o, $mr ) = @_;
	return( "java -jar $picard_path/picard.jar MarkDuplicates I=$i O=$o MAX_RECORDS_IN_RAM=$mr M=$o.dupeMetrics VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=950 " );
}
sub picard_sc_string{
	my ( $mem, $i, $o, $mr ) = @_;
	return( "java -jar $picard_path/picard.jar SortSam I=$i O=$o SO=coordinate MAX_RECORDS_IN_RAM=$mr VALIDATION_STRINGENCY=SILENT " );
}
sub samtools_flagstat_string {
	my ( $mem, $i, $o, $mr ) = @_;
	return( "$samtools_path/samtools flagstat $i > $i.flagstat" );
}

sub is_interleaved {
	my ( $file ) = @_;
	zopen( my $fh, $file ) or die $!;
	my $name_one;
	for my $i ( 1 .. 4 ) {
		if ( ! defined $name_one ) { 
			$name_one = <$fh>;
		}
		else { 
			<$fh>;
		}
	}
	my $name_two = <$fh>;
	for my $name ( $name_one, $name_two ) {
		chomp $name;
		$name =~ s/\s.*//;
		$name =~ s/\/[12]$//;
	}
	my $answer = $name_one eq $name_two ? 1 : 0;
	return( $answer );
}

sub bam_chunk_name {
	my ( $read ) = @_;

  my $out_base = basename $read;
  $out_base =~ s/.gz$//;
  $out_base =~ s/.bz2$//;
  $out_base =~ s/.fq$//;
  $out_base =~ s/.fastq$//;

	return( $out_base );
}

sub setup_mem_reads {
	my ( $r1, $r2, $rg, $args_mem, $mem_i_ref, $mem_script ) = @_;
	alias my $mem_i = $$mem_i_ref;

	my $ob = bam_chunk_name( $r1 );

	if ( ! -e 'bam_dir' ) { mkdir 'bam_dir' or die "can't mkdir bam_dir $!" }

	$ob = "bam_dir/$ob";
	if ( -e "$ob.DONE" ) { unlink "$ob.DONE" }

	my @mem_cmds = ();

	
	push @mem_cmds, mem_cmd( $r1, $r2, $ref_fasta, $rg, $ob, $threads, $args_mem );
	
	if ( $is_paired ) {
    my $fixmate_records = ceil( $MAX_RECORDS * 0.8 );

    push @mem_cmds, picard_fmi_string( $MEM_P_VALUE, "$ob.bam", "$ob.fd.bam", $fixmate_records );
    push @mem_cmds, picard_md_string( $MEM_P_VALUE, "$ob.fd.bam", "$ob.bam", $MAX_RECORDS ) . " && rm $ob.fd.bam";

  }
  else {
		push @mem_cmds, picard_sc_string( $MEM_P_VALUE, "$ob.bam", "$ob.sc.bam", $MAX_RECORDS );
    push @mem_cmds, picard_md_string( $MEM_P_VALUE, "$ob.sc.bam", "$ob.bam", $MAX_RECORDS ) . " && rm $ob.sc.bam";
  }

	write_array_cmds( \@mem_cmds, $mem_i++, $mem_script );
	$mem_i_ref++;
	return( "$ob.bam" );
}

# setup scripts for bwa aln blah, bwa sampe blah ( or samse )

sub setup_aln_backtrack_reads {
	my ( $r1, $r2, $rg, $args_samXe,$Xe_i_ref, $aln_i_ref, $aln_script, $samXe_script ) = @_;
	alias my $samXe_i = $$Xe_i_ref;
	alias my $aln_i   = $$aln_i_ref;

	my $ob    = bam_chunk_name( $r1 );	
  my $name1 = basename $r1;
  my $name2 = basename $r2 if defined $r2;

  if ( ! -e 'sai_dir' ) { mkdir 'sai_dir' or die "can't mkdir sai_dir $!" }
  if ( ! -e 'bam_dir' ) { mkdir 'bam_dir' or die "can't mkdir bam_dir $!" }

  $ob    = "bam_dir/$ob";
  $name1 = "sai_dir/$name1";
  $name2 = "sai_dir/$name2" if ( defined $name2 );

  if ( -e "$name1.sai.log.DONE" ) { unlink "$name1.sai.log.DONE" }
  if ( $name2 && -e "$name2.sai.log.DONE" ) { unlink "$name1.sai.log.DONE" }
  if ( -e "$ob.DONE" ) { unlink "$ob.DONE" } 

	my @aln_cmds = ();
	push @aln_cmds,  align_cmd( $r1, $name1, $ref_fasta, $threads, $args_aln );

  write_array_cmds( \@aln_cmds, $aln_i++, $aln_script );

  if ( $r2 ) { 
		pop @aln_cmds;
		push @aln_cmds,  align_cmd( $r2, $name2, $ref_fasta, $threads, $args_aln );
    write_array_cmds( \@aln_cmds, $aln_i++, $aln_script );
  }
  else { 
    write_array_cmds( 'echo dummy cmd for samse', $aln_i++, $aln_script );  # for keeping job dependencies synced during the occasional SE reads
  }

  my @samXe_cmds = ();

  push @samXe_cmds, samXe_cmd( $r1, $r2, $name1, $name2, $ref_fasta, $rg, $ob, $MEM_P_VALUE, $MAX_RECORDS, $args_samXe ); 

	# 
	# TODO: aka TOREMEMBER max records why fixmate = fixmate_records, maxduplicates = MAX_RECORDS ?
	#

  if ( $is_paired ) {
    # oh FixMate why do you vex me so.
    my $fixmate_records = ceil( $MAX_RECORDS * 0.8 );

		push @samXe_cmds, picard_fmi_string( $MEM_P_VALUE, "$ob.bam", "$ob.fd.bam", $fixmate_records );
		push @samXe_cmds, picard_md_string( $MEM_P_VALUE, "$ob.fd.bam", "$ob.bam", $MAX_RECORDS ) . " && rm $ob.fd.bam";

  }
  else {
		push @samXe_cmds, picard_md_string( $MEM_P_VALUE, "$ob.bam", "$ob.md.bam", $MAX_RECORDS ) . " AS=true && rm $ob.bam && mv $ob.md.bam $ob.bam";
  }

  write_array_cmds( \@samXe_cmds, $samXe_i, $samXe_script );

  $samXe_i += 2;

	return( "$ob.bam" );
  #push @bams, "$ob.bam";
}

sub mem_cmd {
	my ( $read1, $read2, $ref, $rg, $out_bam, $threads, $args_mem ) = @_;

	if ( ! defined $read2 && is_interleaved( $read1 ) ) {
		$args_mem .= " -p";
	}
	my $thread_arg =  $threads      ? " -t $threads " : '';
	my $r2_arg     = defined $read2 ? $read2          : '';

	my $cmd = "$bwa_path/bwa mem $args_mem $thread_arg -R \'$rg\' $ref $read1 $r2_arg  | $samtools_path/samtools view -bt $ref.fai -o $out_bam.bam - ";

	return( $cmd );
}

sub samXe_cmd {
  my ( $read1, $read2, $name1, $name2, $ref, $rg, $out_bam, $MEM_P_VALUE, $MAX_RECORDS, $args ) = @_;
	$args ||= ''; # samse often doesn't have args
  my $bcmd;
  if ( defined $read2 ) { # sampe
    $bcmd = "$bwa_path/bwa sampe -r \'$rg\' $args -P $ref $name1.sai $name2.sai $read1 $read2 | $samtools_path/samtools view -bt $ref.fai -o $out_bam.bam - ";
  }
  else {                  # samse
    $bcmd = "$bwa_path/bwa samse -r \'$rg\' $args $ref $name1.sai $read1 | $samtools_path/samtools view -ut $ref.fai -o - - ";
    $bcmd .= "| java -jar $picard_path/picard.jar SortSam I=/dev/stdin O=$out_bam.bam SO=coordinate MAX_RECORDS_IN_RAM=$MAX_RECORDS VALIDATION_STRINGENCY=SILENT";
  }

  return( $bcmd );
}

sub align_cmd {
  my ( $read, $name, $ref, $threads, $args ) = @_;
  my $cmd = "$bwa_path/bwa aln -t $threads $args -f $name.sai $ref $read 2>&1 > $name.sai.log && touch $name.sai.log.DONE";
  return( $cmd );
}

sub create_bwa_index {
	my ( $ref_fasta, $job_prefix, @stuff ) = @_;
  die "\n\nrun bwa index on the reference first\n\n";
  
  my $cmd  = "$bwa_path/bwa index $ref_fasta 2>&1 > $ref_fasta.idx.log && touch $ref_fasta.idx.DONE";

  print `$cmd`;
  until( -e "$ref_fasta.idx.DONE" ) {
    sleep( $sleep_time );
  } 
  sleep(2);
  unlink( "$ref_fasta.idx.DONE" );
}

sub run_or_die {
  my ( $cmd ) = @_;
  my $retval = `$cmd`;
  if ( ${^CHILD_ERROR_NATIVE} ) {
    die "FATAL: command failed ( $cmd ) ${^CHILD_ERROR_NATIVE} $!";
  }
  return $retval;
}