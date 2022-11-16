use strict;
use warnings;
use Test::More;
use Chemistry::Mol;

plan 'no_plan';
#plan tests => 8;

# low-level test
# some typical cases
my @tests = (
    #symbol, explicit valence, charge, radical, expected result
    [ C => 1, 0, 0, 3 ],
    [ C => 4, 0, 0, 0 ],
    [ C => 1, 1, 0, 2 ],
    [ C => 1, -1, 0, 2 ],
    [ C => 1, 0, 1, 1 ],
    [ C => 1, 0, 3, 1 ],
    [ C => 1, 0, 2, 2 ],
    [ O => 1, -1, 0, 0 ],
    [ O => 1, 1, 0, 2 ],
    [ N => 4, 1, 0, 0 ],
    [ N => 2, 1, 0, 2 ],
    [ N => 2, -1, 0, 0 ],
    [ N => 1, 0, 3, 0 ],
    [ B => 1, -1, 0, 3 ],
);

for my $test (@tests) {
    my $expected = pop @$test;
    my $got = Chemistry::Atom->_calc_implicit_hydrogens(@$test);
    is ($got, $expected, "_calc_implicit_hydrogens(@$test) == $expected");
}

# functional test
my $mol = Chemistry::Mol->new;
my $a1 = $mol->new_atom(symbol => 'C');
my $a2 = $mol->new_atom(symbol => 'O', formal_charge => -1);
my $a3 = $mol->new_atom(symbol => 'N', formal_charge => 1);
$mol->new_bond(atoms => [$a1, $a2]);
$mol->new_bond(atoms => [$a1, $a3], order => 2);
$mol->add_implicit_hydrogens;
is ( $a1->implicit_hydrogens,  1,   'C==1');
is ( $a2->implicit_hydrogens,  0,   'O==0');
is ( $a3->implicit_hydrogens,  2,   'N==2');
