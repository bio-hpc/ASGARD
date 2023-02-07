use strict;
use warnings;

use Test::More;
use Chemistry::Mol;
use Chemistry::File::Dumper;

BEGIN {
    if (eval 'use Test::Exception; 1') {
        plan tests => 9;
        #plan 'no_plan';
    } else {
        plan skip_all => "You don't have Test::Exception installed";
    }
}

# make sure that we throw some expected exceptions

###### MOL ######

throws_ok { Chemistry::Mol->read('t/empty.mol') } 
    qr/guess format/, "unknown format (read)";

throws_ok { Chemistry::Mol->read('no_file.mol', format => 'mdl') } 
    qr/No class installed/, "no class installed";

throws_ok { Chemistry::Mol->read('no_file.mol', format => 'dumper') } 
    qr/open file/, "can't open";

throws_ok { Chemistry::Mol->write('no_file.mol') } 
    qr/guess format/, "unknown format (write)";

###### ATOM ######

my $mol = Chemistry::Mol->read('t/mol.pl');
my ($a1, $a2, $a3) = $mol->atoms;

throws_ok { Chemistry::Atom::angle($a1, $a2) } 
    qr/three/, "three atoms needed for angle";

throws_ok { Chemistry::Atom::dihedral($a1, $a2, $a3) } 
    qr/four/, "four atoms needed for dihedral";

throws_ok { $a1->distance('xyz') } 
    qr/undefined for objects of type/, "distance to non-object";

throws_ok { $a1->angle(1, 2) } 
    qr/not an object/, "angle to non-object";

throws_ok { $a1->angle(bless({}), bless({})) } 
    qr/neither an atom/, "angle to funny object";


