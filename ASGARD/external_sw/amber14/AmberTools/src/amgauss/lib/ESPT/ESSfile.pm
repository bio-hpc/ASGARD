package	ESPT::ESSfile;

use strict;
use warnings;

=head1 NAME

ESPT::ESSfile - Generic Electronic Structure Suite (ESS) file object.

=head1 SYNOPSIS

	package ESPT::MyFile;
	use base qw(ESPT::ESSfile);

	package Main;
	my $object = ESPT::MyFile->new();
	$object->prepare(filename);

=head1 DESCRIPTION

This module is the base class for all of the Electronic Structure Perl Toolkit (ESPT) objects.  It provides the generic 
attributes, methods, and subroutines common to all ESPT objects.

=cut

our $VERSION = '0.03';

### Version History ###
# 0.01	consolidate common properties & methods
# 0.02	moved to ESPT namespace, Time property
# 0.03	generalized get method for Rank N tensors

### To Do ###

=head1 ATTRIBUTES

All attributes are currently read-only and get populated by reading the assigned ESS file.  Attribute values are 
accessible through the B<get()> method unless otherwise noted.

=over 10

=item ALPHA

Total number of alpha electrons.

=item ATOMS

Array of atoms stored as atomic symbols. Array length equals NATOMS.

=item BASIS

Basis set employed in the calculation.

=item BETA

Total number of beta electrons.

=item COMPLETE

Flag indicating ESS job completion (1). Defaults to 0.

=item DEBUG

Debug flag enabling or disabling (0) verbose output. Useful when writing new code.
Defaults to 0. Accessible via the B<debug()> method.

=item FILENAME

Full name of the file, including path if passed, assigned to the object.

=item JOBTYPE

ESS job type stored as capitalized keyword(s). Current keywords are: 

=back

=over 10

=over 5

=over 10

=item SP

Single Point

=item OPT

Optimization

=item OPT FREQ

Optimization & Frequency

=item OPT SP

Optimization & Single Point

=item FREQ

Frequency

=back 

=back

=back

=over 10

=item NATOMS

Total number of atoms.

=item NBASIS

Total number of basis functions in the basis set.

=item PROGRAM

Full name of the ESS program 

=item SPIN

Type of electrons, Alpha or Beta, to analyse. 

=item THEORY

Theory level used in the calculation.

=item TITLE

ESS job title.

=item TIME

Total time for all calculations perfomed in the ESS file. Stored as an array with four elements, [days, hours, minutes, seconds].

=item TYPE

File type. This usually eqaul to the file's extension such as log, out, fchk, etc.

=back

=head1 METHODS

Method parameters denoted in [ ] are optional.

=over 15

=item B<$obj-E<gt>new()>

Creates a new ESSfile object.

=cut

## the object constructor **
        
sub new {
        my $invocant = shift;
        my $class = ref($invocant) || $invocant;
        my $ESSfile = {};

	# Generating ESS
	$ESSfile->{PROGRAM} = undef;

	# File type (log, out, fchk, etc.)
	$ESSfile->{TYPE} = undef;

	# Analysis info
	$ESSfile->{DEBUG} = 0;
        $ESSfile->{FILENAME} = undef;
        $ESSfile->{SPIN} = undef;

        # Calculation info
	$ESSfile->{BASIS} = undef;
	$ESSfile->{COMPLETE} = 0;		# flag indicating job completion
	$ESSfile->{JOBTYPE} = undef;		# type of calculation performed
	$ESSfile->{NBASIS} = undef;
	$ESSfile->{THEORY} = undef;
	$ESSfile->{TIME} = [ 0, 0, 0, 0 ]; 	# Total computation time in days [0], hours [1], minutes [2], seconds [3]
	$ESSfile->{TITLE} = undef;	

	# Molecular info
	$ESSfile->{ALPHA} = undef;
	$ESSfile->{ATOMS} = undef;
	$ESSfile->{BETA} = undef;
	$ESSfile->{NATOMS} = undef;
	
	
	bless($ESSfile, $class);
	return $ESSfile;
}


### Methods ###

=item B<$obj-E<gt>prepare(filename [spin])>

Set FILENAME, and SPIN.  SPIN defaults to Alpha.

=cut

# set filename, spin, & debug options 
sub prepare : method {
        my $ESSfile = shift;
        $ESSfile->{FILENAME} = shift;
        $ESSfile->{SPIN} = shift;
        $ESSfile->{SPIN} ||= "Alpha";
        return;
}



=item B<$obj-E<gt>atomconvert(atom)>

Convert atomic symbols to atomic numbers and vice versa. Atom must 
be a valid atomic symbol or number.

=cut


# Convert atomic numbers to atomic symbols
sub atomconvert : method {
	shift; # remove object reference
	my (@num2sym, %sym2num);
	
	#populate atomic symbols
	@num2sym = ("H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne");
	push @num2sym, ("Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar");
	push @num2sym, ("K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr");
	push @num2sym, ("Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe");
	push @num2sym, ("Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu");
	push @num2sym, ("Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn");
	push @num2sym, ("Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr");
	push @num2sym, ("Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub", "Uut", "Uuq", "Uup");
	
	for (my $i=0; $i<scalar(@num2sym); $i++){
		$sym2num{$num2sym[$i]} = $i+1;
	}
	$sym2num{D} = $sym2num{T} = 1;

	my $value = shift;
	
	# determine what was passed and return the opposite
	if ( $value =~ /^[a-zA-Z]+\Z/ ){
		return $sym2num{$value};
	} else {
		if ( $value <= $#num2sym ) {
			return $num2sym[$value - 1];
		} else {
			return $value;
		}
	}
}

=item B<$obj-E<gt>debug([debuglevel])>

Set or retrieve the DEBUG attribute. Passing a debuglevel equal to 1 enables standard debug 
printing.  Verbose debug printing can be enabled by passing integers greater than one.

=cut

sub debug : method {
	my $ESSfile = shift;

	my $debug = shift;

	if ($debug =~ /\d+/ ) {
		$ESSfile->{DEBUG} = $debug;
		return;
	} else {
		return $ESSfile->{DEBUG};
	}
}

=item B<$obj-E<gt>get(attribute [index1] [index2] ... [indexN])>

Get attribute data stored in an N dimensional tensor. If the tensor indicies are
not passed, then the last value for that attribute will be passesd. If a requested
attribute is not present, then a null string is returned. If the requested datum is
not present then undef will be returned.

=cut

# retrieve values. for array data the last element is returned
# unless otherwise requested
sub get : method {
        my $ESSfile = shift;

        # required parameters : value
        my $value = shift;
	print "Get method called for ", $value, "\n" if $ESSfile->{DEBUG} >= 4;

        # optional parameters : index values
	my @index = @_;

        # warn if requested property is not coded & return empty string
        my $valid = 0;
        foreach my $p ( keys %$ESSfile ) {
                if ( $value eq $p ) {
			$valid = 1;
	        	last;
		} elsif ( $value eq "DEBUG" ) {
			$valid = 2;
			last;
		}
			
        }
        if ( $valid == 0 ) {
                warn "Requested property, $value, is not currently coded.\n$!\n" if $ESSfile->{DEBUG} == 1;
                return "";
        } elsif ( $valid ==2 ) {
		warn "Please use the debug() method to access the DEBUG attribute.\n$!\n";
		return "";
	}

	# build hard reference 
	my $reference = \$ESSfile->{$value};

        # undef if no scalar value is present
        return "undef" unless ( defined $ESSfile->{$value} );

        # check for  multidimensional array in a nested
        # manner to avoid autovivification problems.
        my $i =0;
	while ( ref $$reference eq "ARRAY" ) {
	  # get array length
          $index [$i] ||= scalar(@{$$reference}) - 1 unless ( defined $index[$i] );
       	  return "undef" unless ( defined @{$$reference} [$index[$i]] );

	  # update reference 
	  $reference = \@{$$reference} [$index[$i]];       
	  return $$reference unless ( ref $$reference eq "ARRAY" );
	  $i++;
      	}

	# scalar values	
	return $$reference;

}

=item B<$obj-E<gt>MOdecoder(MO)>

Return the molecular orbital number for the requested MO. MO may be HOMO, LUMO or SHOMO.

=cut

# decode homo, shomo & lumo 
# correct for arrays starting with 0 
sub MOdecoder : method  {
        my $ESSfile = shift;
        my $MO = lc(shift);

        if ($MO eq "homo") {
                return $ESSfile->{HOMO} - 1;
        } elsif ($MO eq "lumo") {
                return $ESSfile->{HOMO};
        } elsif ($MO eq "shomo") {
                return $ESSfile->{HOMO} - 2;
        } else {
                return $MO - 1;
        }
}

### Utilities ###

sub printattributes : method {
	# print out object attributes
	my $ESSfile = shift;
        foreach my $p ( sort keys %$ESSfile ) {
		print "$p\n";
        }

}


1;
__END__

=head1 VERSION

0.03

=head1 AUTHOR

Jason L. Sonnenberg, E<lt>sonnenberg.11@osu.eduE<gt>

=head1 COPYRIGHT

Copyright 2006 by Jason L. Sonnenberg

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. I would like to hear of any
suggestions for improvement.

=cut

