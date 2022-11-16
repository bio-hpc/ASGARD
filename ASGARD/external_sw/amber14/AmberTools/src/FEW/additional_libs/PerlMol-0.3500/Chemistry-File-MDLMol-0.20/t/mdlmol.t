use Test::More tests => 13;
BEGIN { use_ok('Chemistry::File::MDLMol') };

# read a molecule and see if some of the data is ok
my $mol = Chemistry::Mol->read("t/1.mol");

isa_ok($mol, 'Chemistry::Mol', '$mol');
is($mol->name, "trans-Difluorodiazene", "name");
is($mol->attr("mdlmol/line2"), "  -ISIS-            3D", "line2");
is($mol->attr("mdlmol/comment"), "r23 N2F2 FN=NF", "comment");
is($mol->atoms(2)->symbol, "N", "symbol");
is($mol->atoms(3)->y3, 1.2409, "coords");
is($mol->bonds(1)->type, 2, "bond type");
is($mol->bonds(1)->order, 2, "bond order");

# try one with charges and radicals
$mol = Chemistry::Mol->read('t/chg_rad.mol');
is($mol->atoms(1)->formal_radical,      2, 'radical');
is($mol->atoms(1)->implicit_hydrogens,  2, 'implicit_hydrogens');
is($mol->atoms(4)->formal_charge,       -1, 'charge');

# see if we can change the change, write it and parse it back
$mol->atoms(4)->formal_charge(1);
my $s = $mol->print(format => 'mdl');
$mol = Chemistry::Mol->parse($s, format => 'mdl');
is($mol->atoms(4)->formal_charge, 1, 'charge rw');


