package ESPT::Aprmcrd;

use base qw(ESPT::ESSfile);
use strict;
use warnings;

=head1 NAME

ESPT::Aprmcrd - AMBER prmcrd file object.

=head1 SYNOPSIS

	use ESPT::Aprmcrd;

	my $file = Aprmcrd->new();

=head1 DESCRIPTION

This module provides methods to quickly access data contianed in an AMBER prmcrd file
object.  AMBER prmcrd files can only be read currently.

=cut

our $VERSION = '0.01';

### Version History ###
# 0.01	digest prmcrd files from Amber9

=head1 ATTRIBUTES

All attributes are currently read-only and get populated by reading the assigned ESS file.  Attribute values are
accessible through the B<get> method.

=over 15

=item  CARTCOORD

Cartesian coordinates stored as a NATOM x 3 matrix

=back

=head1 METHODS

Method parameters denoted in [] are optional.

=over 15

=item B<$file-E<gt>new()>

Creates a new Aprmcrd object

=cut

## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $prmcrd = ESPT::ESSfile->new();

	$prmcrd->{PROGRAM} = "AMBER";
	$prmcrd->{TYPE} = "prmcrd";

	# molecular info
	$prmcrd->{CARTCOORD} = [];		# Current cartesian coordinates

	bless($prmcrd, $class);
	return $prmcrd;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>

Analyze the spin results in file called filename.  Spin defaults to Alpha.

=cut

# set filename & spin then digest the file
sub analyze : method {
	my $prmcrd = shift;
	$prmcrd->prepare(@_);
	$prmcrd->digest();
	return;
}


## subroutines ##

sub digest {

my $prmcrd = shift;

# flags & counters
my $counter = 0;
my $Titleflag = 1;

# open filename for reading or display error
open(PRMCRDFILE,$prmcrd->{FILENAME}) || die "Could not read $prmcrd->{FILENAME}\n$!\n";

# grab everything which may be useful
while (<PRMCRDFILE>){
	# skip blank lines
	next if /^$/;

	# title; first line of text
	if ( $Titleflag == 1 && /^[\w\d\-\(\)]+/ ) {
		chomp($_);
		s/\s+$//;
		$prmcrd->{TITLE} = $_;
		$Titleflag = 0;
		next;
	}
	# number of atoms
	if ( $Titleflag == 0 && /^\s+(\d+)$/ ) {
		$prmcrd->{NATOMS} = $1;
		next;
	}
        # current cartesian coordinates
        # store in an N x 3 array for the time being
        # switch to PerlMol objects in the future
        if ( /^\s+((?:-*\d+\.\d+\s+){1,6})/ ) {
		my @carts = split /\s+/, $1;
                for (my $i=0; $i<scalar(@carts); $i++) {
                        push @{ $prmcrd->{CARTCOORD} [$counter] }, $carts[$i];
                        $counter++ if $#{$prmcrd->{CARTCOORD} [$counter]}  == 2;
                }
                next;
        }

}
}


1;
__END__
# Below is the documentation for this module.

=head1 VERSION

0.01

=head1 SEE ALSO

F<ESSfile>

=head1 AUTHOR

Jason L. Sonnenberg, E<lt>sonnenberg.11@osu.eduE<gt>

