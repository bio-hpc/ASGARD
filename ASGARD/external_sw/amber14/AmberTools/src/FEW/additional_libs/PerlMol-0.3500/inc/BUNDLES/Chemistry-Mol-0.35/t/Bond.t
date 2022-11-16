use strict;
use warnings;

use Chemistry::Mol;
use Chemistry::File::Formula;
#use Test::More "no_plan";
use Test::More tests => 17;

my $mol  = Chemistry::Mol->parse('CC', format => 'formula');
my $bond = $mol->new_bond(atoms => [ $mol->atoms ] );
my ($a1, $a2) = $mol->atoms;

is ( scalar $mol->atoms,    2,      "mol atom count" );
is ( scalar $mol->bonds,    1,      "mol bond count" );
is ( scalar $bond->atoms,   2,      "bond atom count" );
is ( scalar $a1->bonds,     1,      "atom bond count" );
is ( scalar $a2->bonds,     1,      "atom bond count" );

$bond->delete;
ok ( 1, "deleted the bond" );

is ( scalar $mol->atoms,    2,      "mol atom count" );
is ( scalar $mol->bonds,    0,      "mol bond count" );
is ( scalar $bond->atoms,   2,      "bond atom count" );
is ( scalar $a1->bonds,     0,      "atom bond count" );
is ( scalar $a2->bonds,     0,      "atom bond count" );

$mol->add_bond($bond);
ok ( 1, "readded the bond" );

is ( scalar $mol->atoms,    2,      "mol atom count" );
is ( scalar $mol->bonds,    1,      "mol bond count" );
is ( scalar $bond->atoms,   2,      "bond atom count" );
is ( scalar $a1->bonds,     1,      "atom bond count" );
is ( scalar $a2->bonds,     1,      "atom bond count" );


