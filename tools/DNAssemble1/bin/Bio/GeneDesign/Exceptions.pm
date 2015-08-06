#
# GeneDesign exceptions
#

=head1 NAME

Bio::GeneDesign::Exceptions

=head1 VERSION

Version 5.52

=head1 DESCRIPTION

GeneDesign is a library for the computer-assisted design of synthetic genes

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>

=cut

package Bio::GeneDesign::Exceptions;

use strict;
use warnings;

use Exception::Class
(
  
  "Bio::GeneDesign::Exception::UnBBable" =>
  {
    description => 'No building blocks could be carved from this chunk'
  },
  
  "Bio::GeneDesign::Exception::UnOLable" =>
  {
    description => 'No oligos could be chopped from this building block'
  },
);

1;