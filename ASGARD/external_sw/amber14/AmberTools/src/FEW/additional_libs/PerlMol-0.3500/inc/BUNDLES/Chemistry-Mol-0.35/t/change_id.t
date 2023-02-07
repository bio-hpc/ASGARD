use Test::More;

# These tests make sure that if the id of e.g. an atom changes, the containing
# object is notified

#plan 'no_plan';
plan tests => 5;

use Chemistry::File::Dumper;

my $mol = Chemistry::Mol->read("t/mol.pl");
isa_ok( $mol, 'Chemistry::Mol' );
is( $mol->atoms(1)->id,     'a1',           'id before' );
ok( $mol->atoms(1) == $mol->by_id('a1'),    'id matches before' );

$mol->atoms(1)->id('xyz123');
is( $mol->atoms(1)->id,     'xyz123',           'id after' );
ok( $mol->atoms(1) == $mol->by_id('xyz123'),    'id matches after' )
    or diag sprintf "got %s, expected %s", $mol->atoms(1), $mol->by_id('xyz123');
