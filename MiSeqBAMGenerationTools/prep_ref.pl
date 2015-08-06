#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw/realpath/;
use File::stat;
use GT::Sequtils;

my $gt = GT::Sequtils->new();

my ( $ref_file_in, $ref_file_out, $change_bad_to_n, $scrub_ids, 
     $picard_path, $bwa_path, $samtools_path, 
     $do_bwa_index, $do_dict, $do_cnvnator, $help ) ;

$scrub_ids    = 1;
$do_dict      = 1;
$do_bwa_index = 0;
$do_cnvnator  = 0;

$picard_path = "";
$bwa_path = "";
$samtools_path = "";

GetOptions( 
  'file|ref=s'      => \$ref_file_in,
  'out=s'           => \$ref_file_out,
  'bad_to_n'        => \$change_bad_to_n,
  'scrub_ids!'      => \$scrub_ids,
	'index!'          => \$do_bwa_index,
	'dict!'           => \$do_dict,
	'cnvnator!'       => \$do_cnvnator,
  'picard_path=s'   => \$picard_path,
  'bwa_path=s'      => \$bwa_path,
  'samtools_path=s' => \$samtools_path,
  'help'            => \$help,
) or die $!;

if ( ! defined $ref_file_in && @ARGV ) {
  $ref_file_in = shift;
}

if ( $help || ! defined $ref_file_in || ! -e $ref_file_in ) {
	my $bin = basename $0;

print STDERR<<EOH;

   checks fasta for odd characters and reformats to same line length
   runs samtools faidx 
   runs picard CreateSequenceDictionary 
   optionally runs bwa index 
   Dependencies: samtools, picard, bwa

 Usage: $bin [options] input.fasta

   -bad_to_n     changes non ACGTNacgtn characters to N, probably fine if they are iupac ambiguity codes

   -noscrub_ids  do not change id entry at all 
                 default(scrub_ids) is to make ncbi ids easier to handle for downstream tools
                 >gi|12|thing blah blah => >blah gi|12|thing blah blah  
                 >gi|12|thing.1|        => >thing.1 gi|12|thing.1|

   -index        run bwa index on the fasta         [ default is SKIP THIS ]
   -cnvnator     process fasta for CNVNator         [ default is SKIP THIS ]
   -dict         run picard CreatSequenceDictionary [ default is DO THIS   ]
   -out          specify a particular output file name instead of replacing input file

   -picard_path  Path to Picard jars
   -bwa_path     Path to BWA
   -samtools_path Path to Samtools

   by default the new fasta is the input.fasta name, and original is saved as input.fasta.orig

EOH

  exit;
}

my $swap_filenames = defined $ref_file_out ? 0 : 1;

if ( $swap_filenames ) {
  $ref_file_out = "$ref_file_in.temp";
}
else {
	if ( ! ( $ref_file_out =~ /.fa$|.fasta$/ ) ) {             # broad tools insist on ref with .fa
		$ref_file_out .= '.fa';
	}
}


open( my $infh, $ref_file_in ) or die $!;

# fix for fastas that start with blank lines, will wonders never cease.

my $offset = 0;

while ( <$infh> ) {
  last if (/^>/);
  $offset = tell $infh; 
}
seek( $infh, $offset, 0 );

open( my $ofh, '>', $ref_file_out ) or die $!;

my $warning_given   = 0;
my $bases_converted = 0;

while ( my ( $id, $def, $seq ) = $gt->next_seq( $infh ) ) {

  # scrub weird line endings from defline, split will catch the rest.
  $def =~ s/\r$//;

  my $full_defline = '>';

  # shorten ncbi ids
  if ( $scrub_ids ) {
    if ( $id =~ /^gi\|\d+/ ) {
      $id =~ s/\|$//;
      $id = ( split /\|/, $id ) [ -1 ];
    } 

  # uhm, no id + actual def would be '> blah', a
    if ( ! defined $id ) { warn "couldn't pick out an id from $def -*- check the deflines please -*- ok to manually fix -*-\n"; }

    my ( $pre_mangled_id ) = split /\s+/, $def;

    $full_defline .= defined $id ? $id eq $pre_mangled_id ? '' : $id : '';
    $full_defline .= defined $def ? " $def\n" : "\n";
    $full_defline =~ s/^>\s+/>/;
  }
  else {
    $full_defline  .= "$def\n";
  }

  if ( $seq =~ /[^ACGNTacgnt]/  ) {
    warn "non-[ACGTN] characters found in sequence! \n" unless( $bases_converted );
    if ( $change_bad_to_n ) { 
      warn "Changing them to N\n" unless ( $bases_converted );
      $bases_converted += $seq =~ s/[^ACGNTacgnt]/N/g;
    }
    else {
      my ( $first_bad ) = $seq =~ /([^ACGNTacgnt])/;
      warn "\nrun with -bad_to_n to change them to N's\n\n,probably fine if they were iupac ambiguity codes\n\nFirst bad char was ( $first_bad )\n\n"; 
    }
    if ( ! $seq ) { 
      warn "No valid sequence left for ( seq $def )";
      next;
    }
  }

  
#  my $full_defline = '>' . { defined $id ? $id : '' } . { defined $def ? " $def\n" : "\n" };

  print $ofh $full_defline;

  while ( $seq ) {
    print $ofh uc( substr( $seq, 0, 70, '' ) ), "\n";
  }
}

close( $infh ) or die $!;
close( $ofh  ) or die $!;

if ( $swap_filenames ) {
  rename( $ref_file_in, "$ref_file_in.bak" ) or die $!;
  rename( $ref_file_out, $ref_file_in      ) or die $!; 
	$ref_file_out = $ref_file_in;
}

my $cmd;
$cmd = "$samtools_path/samtools faidx $ref_file_out";
system( $cmd ) && die "ERROR: $cmd failed $!";

## picard chokes if file isn't .fa or .fasta
my $dict_name = $ref_file_out;

if ( $do_dict  && ( $dict_name =~ s/.fa$// || $dict_name =~ s/.fasta$// ) ) { 
	$dict_name .= '.dict';

	if ( -e $dict_name ) {
		unlink( $dict_name );
	}
	print STDERR "\ncreating dictionary file with picard\n\n";


	$cmd = "java -jar $picard_path/picard.jar CreateSequenceDictionary R=$ref_file_out O=$dict_name > /dev/null";
	system( $cmd ) && warn "ERROR: $cmd failed $!";
}
else {
	print STDERR "skipping picard .dict creation ( requires ref to end in .fa or possibly .fasta )\n";
}

if ( $do_bwa_index ) {
	print "\nindexing with bwa -path $bwa_path \n";
  $cmd = "$bwa_path/bwa index $ref_file_out";

	system( $cmd ) 
}
if ( $bases_converted ) {
  my $bases = $bases_converted == 1 ? 'base' : 'bases';
  my $ns    = $bases_converted == 1 ? 'N'    : 'Ns';
  print STDERR "$bases_converted $bases converted to $ns\n";
}

if ( $do_cnvnator ) {

	my $fa_dir  = dirname $ref_file_out;
	my $fa_name = realpath $ref_file_out;

	my $cnv_dir = "$fa_dir/cnvnator_dir";

	if ( ! -e $cnv_dir ) {
 	 mkdir "$cnv_dir" or die "ERROR: can't create $cnv_dir $!";
	}
	else {
		print STDERR "unlinking previous contents of cnvnator_dir\n";
		unlink( glob("$cnv_dir/*.fa") );
	}

	chdir "$cnv_dir" or die $!;

	my $cmd = "fasta_split $fa_name";
	print STDERR "cmd: $cmd\n";
	system( $cmd ) && die "ERROR: $cmd failed $!";
}
