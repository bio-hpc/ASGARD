package ESPT::Aprmtop;

use base qw(ESPT::ESSfile);
use strict;
use warnings;

=head1 NAME

ESPT::Aprmtop - AMBER prmtop file object.

=head1 SYNOPSIS

	use ESPT::Aprmtop;

	my $file = Aprmtop->new();

=head1 DESCRIPTION

This module provides methods to quickly access data contained in an AMBER prmtop file
object.  AMBER prmtop files can only be read currently.

=cut

our $VERSION = '0.02';

### Version History ###
# 0.01	digest prmtop files from Amber9
#
# 0.02	use integer atomic masses for atom 
# 	type determination

=head1 ATTRIBUTES

All attributes are currently read-only and get populated by reading the assigned ESS file.  Attribute values are
accessible through the B<get> method.

=over 15

=item CHARGE

Total molecular charge given by (total number of protons) - (total Molecular Mechanics charge).  
Values rounded to nearest integer.  

=item MMCHARGE

Molecular mechanic atomic charges as determined by AMBER stored as an NATOM array.

=item MULTIPLICITY

Molecular spin multiplicity. Defaults to 2S+1.

=back

=head1 METHODS

Method parameters denoted in [] are optional.

=over 15

=item B<$file-E<gt>new()>

Creates a new Aprmtop object

=cut

## the object constructor **

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $prmtop = ESPT::ESSfile->new();

	$prmtop->{PROGRAM} = "AMBER";
	$prmtop->{TYPE} = "prmtop";

	# molecular info
	$prmtop->{CHARGE} = undef;
	$prmtop->{MMCHARGE} = [];
	$prmtop->{MULTIPLICITY} = undef;

	bless($prmtop, $class);
	return $prmtop;
}


## methods ##

=item B<$file-E<gt>analyze(filename [spin])>
    
Analyze the spin results in file called filename.  Spin defaults to Alpha.

=cut

# set filename & spin then digest the file
sub analyze : method {
	my $prmtop = shift;
	$prmtop->prepare(@_);
	$prmtop->digest();
	return;
}


## subroutines ##

sub digest {

my $prmtop = shift;

# flags & counters
my $atomflag = 0;
my @atomicsymbol;
my $cartflag = 0;
my $carttot = 0;
my $chargeflag = 0;
my $counter = 0;
my $massflag = 0;
my $Titleflag = 0;

# open filename for reading or display error
open(PRMTOPFILE,$prmtop->{FILENAME}) || die "Could not read $prmtop->{FILENAME}\n$!\n";

# grab everything which may be useful
while (<PRMTOPFILE>){
	# skip blank lines
	next if /^$/;

	# title (assume only one line long)
	if ( /^\%FLAG\s+TITLE/ ) {
		$Titleflag = 1;
		next;
	}
	if ( $Titleflag == 1 && /^[\w\d\-\(\)]+/ ) {
		chomp($_);
		s/\s+$//;
		$prmtop->{TITLE} = $_;
		$Titleflag = 0;
		next;
	}
	# Atoms; 
	# Hopefully stored as atomic symbol + an optional id # 
	# seperated by a space(s). Ignore everything 
	# that does not follow this syntax.
        if ( /^%FLAG\s+ATOM_NAME/ ) {   
                $atomflag = 1;   
                $counter = 0;
               next;
        }
        if ( $atomflag == 1 && /^((?:[a-zA-Z]{1,2}\d*\s+){1,20})/ ) {
                my @atomsym = split /\s+/, $1, 20;
                for (my $i=0; $i<scalar(@atomsym); $i++) {
			# drop the trailing numbers and spaces
			$atomsym[$i] =~ s/\d+//;
			$atomsym[$i] =~ s/\s+//;
			$atomicsymbol[$counter] = $atomsym[$i];
                        $counter++;
                }
                next;
        }
        if ( $atomflag == 1 && /^%FLAG/ ) {
                $atomflag = 0;
        }
	# Atomic charges - must be totaled to know molecular charge
	if ( /^%FLAG\s+CHARGE/ ) {
		$chargeflag = 1;
		$counter = 0;
		next;
	}
	if ( $chargeflag == 1 && /^\s+((?:-*\d\.\d+E[-+]\d{1,2}\s+){1,5})/ ) {
		my @charges = split /\s+/, $1;
		for (my $i=0; $i<scalar(@charges); $i++) {
			$prmtop->{MMCHARGE} [$counter] = $charges[$i];
			$counter++;
			
		}
		next;
	}
        if ( $chargeflag == 1 && /^%FLAG/ ) {
                $chargeflag = 0;
        }
	# Pick up Atoms from mass section
	# use integer atomic masses to determine atom
	if ( /^%FLAG\s+MASS/ ) {
		$massflag = 1;
		$counter = 0;
		next;
	}
	if ( $massflag == 1 && /^\s+((?:-*\d\.\d+E[-+]\d{1,2}\s+){1,5})/ ) {
		my @mass = split /\s+/, $1;
		for (my $i=0; $i<scalar(@mass); $i++) {
			$prmtop->{ATOMS} [$counter] = mass2sym(sprintf("%.0f", $mass[$i]));
			$counter++;
			# set NATOMS
			$prmtop->{NATOMS} = $counter;
		}
		next;
	}	
	if ( $massflag == 1 && /^%FLAG/ ) {
		$massflag = 0;
	}
}

# ensure that all atoms have been properly 
# assigned an atomic symbol.
for (my $i=0; $i<$prmtop->{NATOMS}; $i++) {
	$_ = $prmtop->{ATOMS} [$i];
	if ( /ArCa|CoNi|BiPo|CmBk/ ) {
		if ( exists($atomicsymbol[$i]) ) {
			$prmtop->{ATOMS} [$i] = $atomicsymbol[$i];
		}
	}
}

# determine molecular charge
my $chargetot = 0;
my $protontot = 0;
for (my $i=0; $i<$prmtop->{NATOMS}; $i++) {
	$chargetot = $chargetot + $prmtop->{MMCHARGE} [$i];
	$protontot = $protontot + $prmtop->atomconvert($prmtop->{ATOMS} [$i]);
}

# convert to integer units of electrons
$prmtop->{CHARGE} = sprintf("%.0f", $chargetot/18.2223);

# correct for negative zeros
$prmtop->{CHARGE} = 0 if $prmtop->{CHARGE} eq "-0";

# set multiplicity, 2S+1, assuming highest spin state
my $electrons = $protontot - $prmtop->{CHARGE};
$prmtop->{MULTIPLICITY} = $electrons - 2*(int($electrons/2)) + 1;

}

# utility subroutines
sub mass2sym {
        my (%imass, $result);

        %imass = (
                1 => "ArCa",
                4 => "He",
                7 => "Li",
                9 => "Be",
                11 => "B",
                12 => "C",
                14 => "N",
                16 => "O",
                19 => "F",
                20 => "Ne",
                23 => "Na",
                24 => "Mg",
                27 => "Al",
                28 => "Si",
                31 => "P",
                32 => "S",
                35 => "Cl",
                40 => "ArCa", 
                39 => "K",
		45 => "Sc",
		48 => "Ti",
		51 => "V",
		52 => "Cr",
		55 => "Mn",
		56 => "Fe",
		59 => "CoNi",
		64 => "Cu",
		65 => "Zn",
		70 => "Ga",
		73 => "As",
		79 => "Se",
		80 => "Br",
		84 => "Kr",
		85 => "Rb",
		88 => "Sr",
		89 => "Y",
		91 => "Zr",
		93 => "Nb",
		96 => "Mo",
		98 => "Tc",
		101 => "Ru",
		103 => "Rh",
		106 => "Pd",
		108 => "Ag",
		112 => "Cd",
		115 => "In",
		119 => "Sn",
		122 => "Sb",
		128 => "Te",
		127 => "I",
		131 => "Xe",
		133 => "Cs",
		137 => "Ba",
		139 => "La",
		140 => "Ce",
		141 => "Pr",
		144 => "Nd",
		145 => "Pm",
		150 => "Sm",
		152 => "Eu",
		157 => "Gd",
		159 => "Tb",
		162 => "Dy",
		165 => "Ho",
		167 => "Er",
		169 => "Tm",
		173 => "Yb",
		175 => "Lu",
		178 => "Hf",
		181 => "Ta",
		184 => "W",
		186 => "Re",
		190 => "Os",
		192 => "Ir",
		195 => "Pt",
		197 => "Au",
		201 => "Hg",
		204 => "Ti",
		207 => "Pb",
		209 => "BiPo",
		210 => "At",
		222 => "Rn",
		223 => "Fr",
		226 => "Ra",
		227 => "Ac",
		232 => "Th",
		231 => "Pa",
		238 => "U",
		237 => "Np",
		244 => "Np",
		243 => "Am",
		247 => "CmBk",
		251 => "Cf",
		252 => "Es",
		257 => "Fm",
		258 => "Md",
		259 => "No",
		260 => "Lr",
		261 => "Rf",
		262 => "Ha",
        );

        my $value = shift;
        $result = $imass{$value};
	if ( defined($result) ) {
		return $result;

	} else {
		return "X";
	}
}

1;
__END__

=head1 VERSION

0.02

=head1 SEE ALSO

F<ESPT::ESSfile>

=head1 AUTHOR

Jason L. Sonnenberg, E<lt>sonnenberg.11@osu.eduE<gt>

=cut

