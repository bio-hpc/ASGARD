use Test::More;

# These tests try to make sure that objects are destroyed when they
# fall out of scope; these requires avoiding circular strong references

use strict;
use warnings;
use Chemistry::Ring;

plan tests => 10;
#plan 'no_plan';

my $dead_atoms  = 0;
my $dead_bonds  = 0;
my $dead_mols   = 0;
my $dead_rings  = 0;

my $badatom;

{
    # make a cyclopropane
    my $mol = Chemistry::Mol->new;
    $mol->new_atom(symbol => 'C');
    $mol->new_atom(symbol => 'C');
    $mol->new_atom(symbol => 'C');
    $mol->new_bond(atoms => [$mol->atoms(1,2)]);
    $mol->new_bond(atoms => [$mol->atoms(2,3)]);
    $mol->new_bond(atoms => [$mol->atoms(3,1)]);

    isa_ok( $mol, 'Chemistry::Mol' );
    is( scalar $mol->atoms, 3,   'atoms before');

    Chemistry::Ring::aromatize_mol($mol);
    is( $dead_atoms,    0,  "before gc - atoms" );
    is( $dead_bonds,    0,  "before gc - bonds" );
    is( $dead_mols,     0,  "before gc - mols" );
    is( $dead_rings,    0,  "before gc - rings" );
}

is( $dead_atoms,    3,  "after gc - atoms" );
is( $dead_bonds,    3,  "after gc - bonds" );
is( $dead_mols,     1,  "after gc - mols" );
is( $dead_rings,    1,  "after gc - rings" );

sub Chemistry::Mol::DESTROY  { $dead_mols++ }
sub Chemistry::Atom::DESTROY { $dead_atoms++ }
sub Chemistry::Bond::DESTROY { $dead_bonds++ }
sub Chemistry::Ring::DESTROY { $dead_rings++ }

