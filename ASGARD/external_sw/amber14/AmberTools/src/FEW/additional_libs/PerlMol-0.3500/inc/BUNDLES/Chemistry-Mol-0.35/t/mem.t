use Test::More;

# These tests try to make sure that objects are destroyed when they
# fall out of scope; these requires avoiding circular strong references

use strict;
use warnings;

#plan 'no_plan';
plan tests => 8;

use Chemistry::File::Dumper;
my $dead_atoms = 0;
my $dead_bonds = 0;
my $dead_mols = 0;

{
    my $mol = Chemistry::Mol->read("t/mol.pl");
    isa_ok( $mol, 'Chemistry::Mol' );
    is( scalar $mol->atoms, 8,   'atoms before');

    # make sure cloned molecules are also gc'ed
    my $mol2 = $mol->clone;

    # atom deletion garbage collection test
    $mol->atoms(2)->delete;
    is( $dead_atoms,    1,  "delete one atom - atoms" );
    is( $dead_bonds,    4,  "delete one atom - bonds" );
    is( $dead_mols,     0,  "delete one atom - mols" );
}

is( $dead_atoms,    16,  "out of scope - atoms" );
is( $dead_bonds,    14,  "out of scope - bonds" );
is( $dead_mols,     2,  "out of scope - mols" );

sub Chemistry::Mol::DESTROY { $dead_mols++ }
sub Chemistry::Atom::DESTROY { $dead_atoms++ }
sub Chemistry::Bond::DESTROY { $dead_bonds++ }

