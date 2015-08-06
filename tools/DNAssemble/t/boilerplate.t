#!perl -T

use 5.006;
use strict;
use warnings;
use Test::More tests => 8;

sub not_in_file_ok
{
  my ($filename, %regex) = @_;

  open(my $fh, '<', $filename) || die "couldn't open $filename for reading: $!";
  my $ref = do {local $/; <$fh>};
  close $fh;
  my @data = split(m{\n}msx, $ref);

  my %violated;

  foreach (my $line = @data)
  {
    while (my ($desc, $regex) = each %regex)
    {
      if ($line =~ $regex)
      {
        push @{$violated{$desc}||=[]}, $.;
      }
    }
  }
  close $fh;

  if (%violated)
  {
    fail("$filename contains boilerplate text");
    diag "$_ appears on lines @{$violated{$_}}" for keys %violated;
  }
  else
  {
    pass("$filename contains no boilerplate text");
  }
  return;
}

sub module_boilerplate_ok
{
  my ($module) = @_;
  not_in_file_ok($module =>
      'the great new $MODULENAME'   => qr/ - The great new /,
      'boilerplate description'     => qr/Quick summary of what the module/,
      'stub function definition'    => qr/function[12]/,
  );
  return;
}

TODO:
{
  local $TODO = "Need to replace the boilerplate text";

  not_in_file_ok(README =>
    "The README is used..."       => qr/The README is used/,
    "'version information here'"  => qr/to provide version information/,
  );

  not_in_file_ok(Changes =>
    "placeholder date/time"       => qr(Date/time)
  );

  module_boilerplate_ok('lib/Bio/DNAssemble.pm');
  module_boilerplate_ok('lib/Bio/DNAssemble/Blast.pm');
  module_boilerplate_ok('lib/Bio/DNAssemble/Oligo.pm');
  module_boilerplate_ok('lib/Bio/DNAssemble/Palindrome.pm');
  module_boilerplate_ok('lib/Bio/DNAssemble/Vmatch.pm');
  module_boilerplate_ok('lib/Bio/DNAssemble/Vector.pm');
}

