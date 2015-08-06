#!/usr/bin/env perl
use strict;
use warnings;
use Data::Alias;

sub fqopen {
  alias my ( $fh, $file ) = @_;
 
  my $result;
  my $open = $file =~ /[.]gz$/  ? "gzip  -cd $file |"
           : $file =~ /[.]bz2$/ ? "bzip2 -cd $file |" 
           : $file;
 
  return( open( $fh, $open ) );
}

my $fh;

for my $file ( @ARGV ) {
  fqopen( $fh, $file ) or die "can't open $file $!";
  my ( $n1 ) = split /\s+/, <$fh> ;
  for ( 0 .. 2 ){
    <$fh>;
  }
  my ( $n2 ) = split /\s+/, <$fh>;
  $n1 =~ s/\/[12]//;
  $n2 =~ s/\/[12]//;
  if ( $n1 eq $n2 ) {
    warn "PAIRED";
    print "paired $file\n";
  }
  else {
    print "single $file\n";;
  }
}