use Test::More;

# These tests try to make sure that objects are destroyed when they
# fall out of scope; these requires avoiding circular strong references

use strict;
use warnings;
use Chemistry::File::SMILES;

plan tests => 8;
#plan 'no_plan';

my $dead_atoms  = 0;
my $dead_bonds  = 0;
my $dead_mols   = 0;

my $badatom;

{
    my $mol = Chemistry::Mol->parse('CC', format => 'smiles');

    isa_ok( $mol, 'Chemistry::Mol' );
    is( scalar $mol->atoms, 2,   'atoms before');

    is( $dead_atoms,    0,  "before gc - atoms" );
    is( $dead_bonds,    0,  "before gc - bonds" );
    is( $dead_mols,     0,  "before gc - mols" );
}

is( $dead_atoms,    2,  "after gc - atoms" );
is( $dead_bonds,    1,  "after gc - bonds" );
is( $dead_mols,     1,  "after gc - mols" );

sub Chemistry::Mol::DESTROY  { $dead_mols++ }
sub Chemistry::Atom::DESTROY { $dead_atoms++ }
sub Chemistry::Bond::DESTROY { $dead_bonds++ }

