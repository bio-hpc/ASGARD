use Test::More;

# this tests make sure that safe_clone actually increments the ids for the
# objects in the new (cloned) molecule

use strict;
use warnings;

#plan 'no_plan';
plan tests => 9;

use Chemistry::File::Dumper;

my $mol = Chemistry::Mol->new;
$mol->new_atom(symbol => 'C');
$mol->new_atom(symbol => 'C');
$mol->new_bond(atoms => [$mol->atoms(1,2)]);

is( $mol->atoms(1)->id,   'a1',   'atom(1) before clone');
is( $mol->bonds(1)->id,   'b1',   'bond(1) before clone');
is( $mol->id,             'mol1', 'mol->id before clone');

my $mol2 = $mol->clone;
is( $mol2->atoms(1)->id,   'a1',   'atom(1) after clone');
is( $mol2->bonds(1)->id,   'b1',   'bond(1) after clone');
is( $mol2->id,             'mol1', 'mol->id after clone');

$mol2 = $mol->safe_clone;
is( $mol2->atoms(1)->id,   'a3',   'atom(1) after safe_clone');
is( $mol2->bonds(1)->id,   'b2',   'bond(1) after safe_clone');
is( $mol2->id,             'mol2', 'mol->id after safe clone');

