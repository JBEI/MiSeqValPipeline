package GT::Sequtils;

use Data::Alias;

# grab bag of junk

sub new {
  my ( $class, %args ) = @_;
  my $self = bless { %args }, $class;
}

sub next_seq {

    my ( $self, $fh ) = @_;
 
    local $/ = "\n>";
 
    return unless( defined( my $entry = <$fh>));
 
    chomp($entry);
 
    $entry =~ s/^>//;
 
    my ($defline, $sequence) = split(/\n/,$entry,2);
 
    my ( $id, $fulldesc) = split /\s+/, $defline, 2;
 
    if ( defined $id && $id eq '') {   # in case of space '> blah'
      $id = $fulldesc;
    }
    elsif ( ! defined $fulldesc ) {
      $fulldesc = '';
    }

# remove spaces and all manner of newlines ( mac, unix, dos )
 
    defined $sequence && $sequence =~ tr/ \t\n\r//d;
#    return( $defline, $sequence ) ;
    return( $id, $defline, $sequence );
 
}

sub next_qual {
  my ( $self, $fh ) = @_;
  local $/ = "\n>";
  return unless( defined( my $entry = <$fh> ) );
  chomp( $entry );
  $entry =~ s/^>//;
  
  my ( $defline, $sequence ) = split( /\n/, $entry, 2 );
  my ( $id, $fulldesc) = split /\s+/, $defline, 2;

  if ( defined $id && $id eq '') {   # in case of space '> blah'
    $id = $fulldesc;
  }
  elsif ( ! defined $fulldesc ) {
    $fulldesc = '';
  }

# remove spaces and all manner of newlines ( mac, unix, dos )
  
  return( $id, $defline, $sequence );
}

sub print_fasta {  # arguments in order of what's most required, next_seq should switch to this.

  my ( $self, $seq, $id, $defline ) = @_;
  
  my $full_defline = '>' . ( defined $id ? $id : '' ) . ( defined $defline ? " $defline\n" : "\n" );
  
  print $full_defline;
  
  while ( $seq ) {
     print substr( $seq, 0, 50, '' ), "\n";
  }
}

#destructive, alters argument
# if given a ref, operate on arg, else operate on copy
sub revc {
  my $is_ref = ref $_[1];

# lame, alias $seq_arg = $_[1]; alias my $seq = $is_ref ? $$seq_arg : $seq_arg; is easier to read

  alias my $seq = $is_ref ? ${$_[1]} : $_[1]; 

  # ** note N,S,W skipped as they are they vain and self-complementary

  $seq =  reverse($seq);
  $seq =~ tr/acgtrymkhbvdACGTRYMKHBVD/tgcayrkmdvbhTGCAYRKMDVBH/;

  return $seq unless $is_ref; 

}

1;
