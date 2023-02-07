package ESPT::Gfchk;

use base qw(ESPT::ESSfile);
use ESPT::Glib 0.01;
use strict;
use warnings;

=head1 NAME

ESPT::Gfchk - Gaussian formatted checkpoint file object.

=head1 SYNOPSIS

   use ESPT::Gfchk;

   my $fchk = Gfchk->new();

=head1 DESCRIPTION

This module provides methods to quickly access data contianed in a Gaussian 
formatted checkpoint fileobject.  Guassian formatted checkpoint files can 
only be read currently.

=cut

our $VERSION = '0.02';

### Version History ###
# 0.01	digest opt freq fchk files
# 0.02	use ESPT namespace
#	use ESPT::Glib module; redundant internal
#	coordinates, gradient and H; SCF and total energies
#	

##### NOTE #####
# Gaussian stores all values in a.u. and standard orientation in the fchk file
#
# Distance = bohrs
# Energy = hartree
# Gradient = hartree/bohr
# Hessian = hartree/bohr^2
################

### To Do ###
# -Store only lower triangle of square matrices and
#  have the get or accessor method return upper triangle as requested.


=head1 ATTRIBUTES

All attributes are currently read-only and get populated by reading the assigned ESS file.  Attribute values 
are accessible through the B<$Gfchk-E<gt>get()> method.

=over 15

=item C

Coefficient matrix, NBASIS x NBASIS. The coefficients correspond to Alpha or
Beta depending upon what spin was passesd to B<$Gfchk-E<gt>analyze()>.

=item CARTCOORD

NATOMS x 3 matrix containing the current cartesian coordinates

=item CHARGE

Total molecular charge.

=item EELEC

Electronic energy.

=item ESCF

SCF energy. This will be either the Hartree-Fock or the DFT energy. See Gaussian documentation
for more information.

=item EIGEN

Array of length NBASIS, containing the eigenvalues.  The eigenvalues correspond to Alpha or 
Beta depending upon what spin was passesd to B<$Gfchk-E<gt>analyze()>.

=item EINFO

Text description of the energy contained in ENERGY.

=item ENERGY

Molecular energy as described by EINFO.  Defaults to SCF electronic energy.

=item FUNCTIONAL

DFT functional utlized in this job.

=item GRADIENT

Array containing the Cartesian gradients.

=item HESSIAN

Lower-triangular matrix containing the Cartesian Hessian.

=item HOMO

Number corresponding to the highest occupied molecular orbital. The value corresponds
to either Alpha or Beta electrons depending upon what spin was passesd to 
B<$Gfchk-E<gt>analyze()>.

=item IRCCOORD

A rank three tensor containing Cartesian coordinates for each IRC geometry.

=item IRCENERGY

Array, with length equal to IRCPOINTS, containing the electronic energy at each IRC geometry.

=item IRCGRADIENT

Cartesian gradients for each IRC geometry stored as a rank two tensor.

=item IRCPOINTS

Total number of IRC steps.

=item IRCSTEP

Array of reaction coordinate values for each IRC step.

=item KEYWORDS

Array containing keywords used in this job.

=item MASS

Array of length NATOMS, containing the atomic masses.

=item MULTIPLICITY

Molecular spin multiplicity, 2S+1.

=item NRINT

Total number of redundant internal coordinates.

=item RINTCOORD

N x 4 matrix containing the redundant internal coordinates. Each coordinate
is defined by four integers corresponding to the atom numbers.  Bond coordinates
have zeros in columns 3 & 4.  Bond angle coordinates have a zero in column 4.

=item ROUTE

Gaussian route line

=item SSQUARED

<S**2> expectation value.

=back

=head1 METHODS

Method parameters denoted in [] are optional.

=over 15

=item B<$file-E<gt>new()>

Creates a new Gfchk object

=cut

## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $fchk = ESPT::ESSfile->new();

	$fchk->{PROGRAM} = "Gaussian";
	$fchk->{TYPE} = "fchk";

	# Link 0 & Route commands
	$fchk->{ROUTE} = undef;
	$fchk->{KEYWORDS} = [];

	# calc info
	$fchk->{FUNCTIONAL} = undef;

	# IRC data
	$fchk->{IRCCOORD} = [];		# Coordinates for each IRC Geometry
	$fchk->{IRCENERGY} = [];	# Energy at each IRC Geometry 
	$fchk->{IRCSTEP} = [];		# Reaction coordinate value at each IRC step
	$fchk->{IRCGRADIENT} = [];	# Gradient at each IRC Geometry
	$fchk->{IRCPOINTS} = 0;		# Total number of IRC Geometries

	# molecular info
	$fchk->{C} = [];		# coefficient matrix
        $fchk->{CARTCOORD} = [];        # Current cartesian coordinates
	$fchk->{CHARGE} = undef;
	$fchk->{EIGEN} = [];
	$fchk->{EELEC} = undef;		# electronic energy
	$fchk->{ESCF} = undef;		# SCF energy
	$fchk->{ENERGY} = undef; 	# total energy 
	$fchk->{EINFO} = "E(elec)";	# total energy description
	$fchk->{GRADIENT} = [];
	$fchk->{HESSIAN} = [];		# lower triangle only
	$fchk->{HOMO} = undef;
	$fchk->{MULTIPLICITY} = undef;
	$fchk->{MASS} = undef;		# atomic masses
	$fchk->{NRINT} = 0;		# number of red. internals
	$fchk->{RINTCOORD} = [];	# Redundant internal coordinates
	$fchk->{SSQUARED} = [];		# <S**2>

	bless($fchk, $class);
	return $fchk;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>
    
Analyze the spin results in file called filename.  Spin defaults to Alpha.

=cut

# set filename & spin then digest the file
sub analyze : method {
	my $fchk = shift;
	$fchk->prepare(@_);
	$fchk->digest();
	return;
}


## subroutines ##

sub digest {

my $fchk = shift;

# flags & counters
my $atomflag = 0;
my $atomtot = 0;
my $cartflag = 0;
my $carttot = 0;
my $col = 0;
my $counter = 0;
my $eigenflag = 0;
my $eigentot = 0;
my $geomcount = 0;
my $gradientflag = 0;
my $hessflag = 0;
my $hesstot = 0;
my $irccoord = 0;
my $ircdata = 0;
my $ircflag = 0;
my $ircgeomflag = 0;
my $ircgradientflag = 0;
my $ircresults = 0;
my $ircresultsflag = 0;
my $redinttot = 0;
my $redintflag = 0;
my $row = 0;
my $rparsed = 0;
my $Titleflag = 0;
my $Cflag = 0;
my $weightflag = 0;

# open filename for reading or display error
open(FCHKFILE,$fchk->{FILENAME}) || die "Could not read $fchk->{FILENAME}\n$!\n";

# quick check for IRC data since IRC data ends up at 
# the end of the fchk file.  This is not the most 
# elegant way to do this, but works for now. -JLS
while (<FCHKFILE>) {
        # skip blank lines
        next if /^$/;

	# check for IRC data and get # of IRC points
	# grow IRC arrays according to results
	# remember that arrays indices start at 0
	if ( /^IRC\sNumber\sof\sgeometries\s+I\s+N=\s+\d+/ ) {
		$ircflag = 1;
		next;
	}
	if ( $ircflag == 1 && /^\s+(\d+)/ ) {
		$ircflag = 0;
		$fchk->{IRCPOINTS} = $1;
		$#{$fchk->{IRCCOORD}} = $1 - 1;
		$#{$fchk->{IRCGRADIENT}} = $1 - 1;
		last;
	}
}

# rewind file
seek(FCHKFILE, 0, 0);		

# grab everything which may be useful
while (<FCHKFILE>){
	# skip blank lines
	next if /^$/;

	# title which is the first line
	# only the first 72 characters are present as of G03
	if ( $Titleflag == 0 ) {
		chomp($_);
		s/\s+$//;
		$fchk->{TITLE} = $_;
		$Titleflag = 1;
		next;
	}
	# Job Type, Method, & Basis
	if ( $rparsed == 0 && /^[SPOFIRCMADBVGa-z=]+/ ){
		chomp($fchk->{ROUTE} = lc($_));
		rparser($fchk);
		$rparsed = 1;
		next;
	}
        # Number of Atoms
        if ( /^Number\s+of\s+atoms\s+I\s+(\d+)/ ) {   
                $fchk->{NATOMS} = $1;
		next;
        }
        # charge
        if ( /^Charge\s+I\s+(-*\d+)/ ) {   
                $fchk->{CHARGE} = $1;
		next;
        }
        # multiplicity 
        if ( /^Multiplicity\s+I\s+(\d+)/ ) {   
                $fchk->{MULTIPLICITY} = $1;
		next;
        }
        # electrons
	# figure HOMO & LUMO, alphas fill first
        if ( /^Number\s+of\s+alpha\s+electrons\s+I\s+(\d+)/ ) {
                $fchk->{ALPHA} = $1;
		next;
        }
        if ( /^Number\s+of\s+beta\s+electrons\s+I\s+(\d+)/ ) {
                $fchk->{BETA} = $1;
		next;
        }
        # basis functions 
        if ( /^Number\s+of\s+basis\s+functions\s+I\s+(\d+)/ ) {
		$fchk->{NBASIS} = $1;
		next;
        }
	# SCF energy 
	if ( /^SCF\s+Energy\s+R\s+(-*\d+\.\d+E-*\+*\d{2,})/ ) {
		$fchk->{ESCF} = $1;
		next;
	}
	# Total electronic energy
	if ( /^Total\s+Energy\s+R\s+(-*\d+\.\d+E-*\+*\d{2,})/ ) {
		$fchk->{EELEC} = $1;
		$fchk->{ENERGY} = $1;
		next;
	}
	# <S**2>; use 2D array to conform with other ESPT modules
	if ( /^S\*\*2\s+R\s+(\d\.\d+E\+\d{2,})/ ) {
		$fchk->{SSQUARED} [0] = $1;
		next;
	}
	# Atoms; stored as atomic numbers 
	if ( /^Atomic\s+numbers\s+I\s+N=\s+(\d+)/ ) {
		$atomtot = $1;
		$atomflag = 1;
		$counter = 0;
		next;
	}
	if ( $atomflag == 1 && /^\s+((?:\d+\s+){1,6})/ ) {
		my @atomnum = split /\s+/, $1;
		for (my $i=0; $i<scalar(@atomnum); $i++) {
			$fchk->{ATOMS} [$counter] =  $fchk->atomconvert($atomnum[$i]);
			$counter++;
			$atomflag = 0 if $counter == $atomtot;
		}
		next;
	}
	# Redundant internal coordinates
	# not present for all calculations
	if ( /^Number\sof\sredundant\sinternal\sbonds\s+I\s+(\d+)/ ) {
		$fchk->{NRINT} = $1;
		next;
	}
        if ( /^Number\sof\sredundant\sinternal\sangles\s+I\s+(\d+)/ ) {
                $fchk->{NRINT} = $fchk->{NRINT} + $1;
                next;
        }
        if ( /^Number\sof\sredundant\sinternal\sdihedrals\s+I\s+(\d+)/ ) {
                $fchk->{NRINT} = $fchk->{NRINT} + $1;
                next;
        }
	if ( /^Redundant\sinternal\scoordinate\sindices\s+I\s+N=\s+(\d+)/ ) {
		$redinttot= $1;
		$redintflag = 1;
		$counter=0;
		next;
	}
	if ( $redintflag == 1 && /^\s+((?:-*\d+\s+){1,6})/ ) {
		my @ints = split /\s+/, $1;
		for (my $i=0; $i<scalar(@ints); $i++) {
			push @{$fchk->{RINTCOORD} [$counter] }, $ints[$i];
			$counter++ if $#{$fchk->{RINTCOORD} [$counter]} == 3;
			$redintflag = 0 if $counter*4 == $redinttot;
		}
		next;
	}
	# current cartesian coordinates
	# store in an N x 3 array
	if ( /^Current\s+cartesian\s+coordinates\s+R\s+N=\s+(\d+)/ ) {
		$carttot = $1;
		$cartflag = 1;
		$counter = 0;
		next;
	}
	if ( $cartflag == 1 && /^\s+((?:-*\d\.\d+E[-\+]\d{2,}\s+){1,5})/ ) {
		my @carts = split /\s+/, $1;
		for (my $i=0; $i<scalar(@carts); $i++) {
			push @{ $fchk->{CARTCOORD} [$counter] }, $carts[$i];
			$counter++ if $#{$fchk->{CARTCOORD} [$counter]}  == 2;
			$cartflag = 0 if $counter*3 == $carttot;
		}
		next;
	} 
	# Real atomic weights
	if ( /Real\s+atomic\s+weights\s+R\s+N=\s+(\d+)$/ ) {
		$weightflag = 1;
		$counter = 0;
		next;
	}
	if ( $weightflag ==1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/ ) {
		my @weights = split /\s+/, $1;
		for (my $i=0; $i<scalar(@weights); $i++) {
			$fchk->{MASS} [$counter] = $weights[$i];
			$counter++;
			$weightflag = 0 if $counter == $fchk->{NATOMS};
		}
		next;
	}
        # Eigenvalues; occur only once per spin
	# must still use 2D array to stay in line with other ESPT modules
        if ( /^$fchk->{SPIN}\s+Orbital\s+Energies\s+R\s+N=\s+(\d+)$/ ) {
		$eigentot = $1;
		$eigenflag = 1;
		$counter = 0;
		next;
	}
	if ( $eigenflag == 1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/) {
		my @eig = split /\s+/, $1;
                for (my $i=0; $i<scalar(@eig); $i++) {
			$fchk->{EIGEN} [0] [$counter] = $eig[$i];
			$counter++;
			$eigenflag = 0 if $counter == $eigentot;
                }
		next;
        }
#Fix	# MO coeffients (square matrix, multiple occurances)
#	if ( /\s+($fchk->{SPIN})*Molecular Orbital Coefficients/ ) {
#		$Cflag = 1;
#		$fchk->{C} = undef;
#		$counter = 0;
#		next;
#	}
#	if ( $Cflag == 1 && /\s*(\d+)\s(\d+)\s*(\w+)\s+(\d+[A-Z]+\s?\-?\+?[0-9]*)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*/ ) {
#		$fchk->{BASISLABELS}[$counter] = [$1, $2, $3, $4];
#		push @{ $fchk->{C}[$counter] }, $5, $6, $7, $8, $9;
#		$counter++;
#		$counter = 0 if $counter == $fchk->{NBASIS};
#		next;
#	} elsif ( $Cflag == 1 && /\s*(\d+)\s*(\d+[A-Z]+\s?\-?\+?[0-9]*)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*/ ) {
#             	$fchk->{BASISLABELS}[$counter] = [$1, $fchk->{BASISLABELS}[$counter - 1] [1], $fchk->{BASISLABELS}[$counter - 1] [2], $2];
#		push @{ $fchk->{C}[$counter] }, $3, $4, $5, $6, $7;
#              	$counter++;
#              	$counter = 0 if $counter == $fchk->{NBASIS};
#		next;
#	}
	# Cartesian Gradient (3*NATOMS x 1 vector)
        if ( /Cartesian\s+Gradient\s+R\s+N=\s+(\d+)$/ ) {
                $gradientflag = 1;
                $counter = 0;
                next;
        }
        if ( $gradientflag == 1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/ ) {
                my @gradients = split /\s+/, $1;
                for (my $i=0; $i<scalar(@gradients); $i++) {
                        $fchk->{GRADIENT} [$counter] = $gradients[$i];
                        $counter++;
                        $gradientflag = 0 if $counter == 3*$fchk->{NATOMS};
                }
                next;
        }
	# Cartesian Hessian
	# 3*NATOMS x 3*NATOMS matrix delivered as a lower
	# triangular matrix (3*NATOMS)(3*NATOMS+1)(1/2) elements
        if ( /Cartesian\s+Force\s+Constants\s+R\s+N=\s+(\d+)$/ ) {
		$hesstot = $1;
                $hessflag = 1;
                $counter = $row = $col = 0;
                next;
        }
        if ( $hessflag == 1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/ ) {
                my @hessian = split /\s+/, $1;
                for (my $i=0; $i<scalar(@hessian); $i++) {
                        $fchk->{HESSIAN} [$row] [$col] = $hessian[$i];
                        $fchk->{HESSIAN} [$col] [$row] = $hessian[$i] unless $row == $col;
                        $counter++;
			if ( $row == $col ) {
				$col++;
				$row = -1;
			}
			$row++;
                        $hessflag = 0 if $counter == $hesstot;
                }
                next;
        }
	# IRC data per geom
	if ( /^IRC\sNum\sresults\sper\sgeometry\s+I\s+(\d+)/ ) {
		$ircresults = $1;
		next;
	}
	# IRC # of coordinates
	if ( /^IRC\sNum\sgeometry\svariables\s+I\s+(\d+)/ ) {
		$irccoord = $1;
		next;
	}
	# IRC energy and step value
	# These steps are sequential and proceed either 
	# forward or backward from the TS which is 0.0 
	# along the IRC. Energies come first followed by 
	# the cooresponging IRC value. Currently all IRC
	# values are given as positive regardless of whether
	# the displacement is towards products or reactants.
	if ( /^IRC\spoint\s+\d+\sResults\sfor\seach\sgeom.*\s+R\s+N=\s+\d+/ ) {
		$ircresultsflag = 1;
		$counter = 0;
		next;
	}
	if ( $ircresultsflag == 1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/ ) {
                my @results = split /\s+/, $1;
                for (my $i=0; $i<scalar(@results); $i++) {
                        $counter++;
			# use modulo math to determine if this is an energy or step value
			if ( $counter % 2 ) {
				push(@{$fchk->{IRCENERGY}}, $results[$i]);
			} else {
				push(@{$fchk->{IRCSTEP}}, $results[$i]);

			}
                }
		$ircresultsflag = 0 if $counter == $ircresults * $fchk->{IRCPOINTS};
                next;
	}	
	# IRC geometries
	if ( /^IRC\spoint\s+\d+\sGeometries\s+R\s+N=\s+\d+/ ) {
		$ircgeomflag = 1;
		$geomcount = 0;
		$counter = 0;
		next;
	}
	if ( $ircgeomflag == 1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/ ) {
                my @coords = split /\s+/, $1;
                for (my $i=0; $i<scalar(@coords); $i++) {
                        push @{ $fchk->{IRCCOORD} [$geomcount] [$counter] }, $coords[$i];
                        $counter++ if $#{$fchk->{IRCCOORD}[$geomcount] [$counter]}  == 2;
			if ( $counter == $fchk->{NATOMS} ) {
				$geomcount++;
				$counter = 0;
			}
                        $ircgeomflag = 0 if $geomcount == $fchk->{IRCPOINTS};
                }
                next;
        }
	# IRC gradients (3*NATOMS x 1 vector)
        if ( /IRC\spoint\s+\d+\sGradient\sat\seach\sgeom.*\s+R\s+N=\s+\d+/ ) {
                $ircgradientflag = 1;
		$geomcount = 0;
                $counter = 0;
                next;
        }
        if ( $ircgradientflag == 1 && /^\s+((?:-*\d\.\d+E-*\+*\d{2,}\s+){1,5})/ ) {
                my @gradients = split /\s+/, $1;
                for (my $i=0; $i<scalar(@gradients); $i++) {
                        $fchk->{IRCGRADIENT} [$geomcount] [$counter] = $gradients[$i];
                        $counter++;
			if ( $counter == 3*$fchk->{NATOMS} ) {
				$geomcount++;
				$counter = 0;
			} 
                        $ircgradientflag = 0 if $geomcount == $fchk->{IRCPOINTS};
                }
                next;
        }
}		

# set HOMO
	$fchk->{HOMO} = $fchk->{uc($fchk->{SPIN})};
}


# convert Scientific notation to decimals
sub sci2dec {
	my $value = shift;
	return unless defined $value;
	
	my $dec = 1*$value;
	return $dec;
}

1;
__END__

=head1 VERSION

0.02

=head1 SEE ALSO

F<ESPT::ESSfile>

=head1 AUTHOR

Jason L. Sonnenberg, E<lt>sonnenberg.11@osu.edu<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2006 by Jason L. Sonnenberg

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. I would like to hear of any
suggestions for improvement.

=cut
