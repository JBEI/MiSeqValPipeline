#!perl -T

use Test::More tests => 4;
use strict;
use warnings;

BEGIN
{
  use_ok( 'Bio::DNAssemble' )               || print "Can't use DNAssemble\n";
  use_ok( 'Bio::DNAssemble::Blast' )        || print "Can't use Blast\n";
  use_ok( 'Bio::DNAssemble::Oligo' )        || print "Can't use Oligo\n";
  use_ok( 'Bio::DNAssemble::Vector' )       || print "Can't use Vector\n";
}

