# Tests for geometry-related methods such as distance, angle, and dihedral

use strict;
use warnings;
use Test::More;
use Chemistry::File::Dumper;
use Math::VectorReal;

#plan 'no_plan';
plan tests => 10;

my $mol = Chemistry::Mol->read("t/mol.pl");
isa_ok( $mol, 'Chemistry::Mol' );

my (@a);

# angle
@a = $mol->atoms(1,2,3);
is( scalar @a, 3, 'three atoms');
is_float( Chemistry::Atom::angle_deg(@a), 110.7, 0.1, "angle" );

# dihedral
@a = $mol->atoms(1,2,3,4);
is_float( Chemistry::Atom::dihedral_deg(@a), -85.6,  0.1, "dihedral" );

# ill-defined angles and dihedrals

my $v0 = vector(0,0,0);
my $v1 = vector(1,0,0);
my $v2 = vector(2,0,0);
my $v3 = vector(3,0,0);

is( Chemistry::Atom::angle($v0, $v1, $v0), 0, "zero angle" );
is( Chemistry::Atom::angle($v0, $v0, $v0), 0, "bad angle" );
is( Chemistry::Atom::angle($v0, $v0, $v1), 0, "bad angle" );
is_float( Chemistry::Atom::angle_deg($v0, $v1, $v2), 180, 0.1, "linear angle" );
is( Chemistry::Atom::dihedral($v0, $v0, $v0, $v0), 0, "bad dihedral" );
is( Chemistry::Atom::dihedral($v0, $v1, $v2, $v3), 0, "bad dihedral" );


#############

sub is_float {
    my ($got, $expected, $tol, $name) = @_;
    ok( abs($got - $expected) < $tol, $name ) 
        or diag "got $got, expected $expected";
}

