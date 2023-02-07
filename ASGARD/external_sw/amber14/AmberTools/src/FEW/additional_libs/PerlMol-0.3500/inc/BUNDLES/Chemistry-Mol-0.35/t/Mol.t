#use Test::More "no_plan";
use Test::More tests => 20;
BEGIN { 
    use_ok('Chemistry::Mol');
};

# Constructors
my $mol = Chemistry::Mol->new;
isa_ok($mol, 'Chemistry::Mol', '$mol');
isa_ok($mol, 'Chemistry::Obj', '$mol');
my $atom = Chemistry::Atom->new(Z => 6, coords => [0, 0, 3], name => "carbon");
isa_ok($atom, 'Chemistry::Atom', '$atom');
isa_ok($atom, 'Chemistry::Obj', '$atom');
my $atom2 = Chemistry::Atom->new(Z => 8, coords => [4, 0, 0], id => 'xyz');
my $bond = Chemistry::Bond->new(atoms => [$atom, $atom2], type => '=');
isa_ok($bond, 'Chemistry::Bond', '$bond');
isa_ok($bond, 'Chemistry::Obj', '$bond');

# Mol methods
$mol->add_atom($atom, $atom2);
is(scalar $mol->atoms, 2, '$mol->atoms');
ok($mol->atoms(1) eq $atom, '$mol->atoms(1) eq $atom');
ok($mol->by_id('xyz') eq $atom2, '$mol->by_id');
ok($mol->atoms_by_name('carbon') eq $atom, '$mol->atoms_by_name');
$mol->add_bond($bond);
is(scalar $mol->bonds, 1, '$mol->bonds');
ok($mol->bonds(1) eq $bond, '$mol->bonds(1) eq $bond');
my $atom3;
ok($atom3 = $mol->new_atom(symbol => "N"), '$mol->new_atom');


# mass
$mol = Chemistry::Mol->new;
$mol->new_atom(symbol => 'O');
$mol->new_atom(symbol => 'H');
$mol->new_atom(symbol => 'H');
ok(abs($mol->mass - 18.01528) < 0.0001, '$mol->mass');
$mol->atoms(1)->mass(18); 
ok(abs($mol->mass - 20.015) < 0.01, '$mol->mass');

# sprout_hydrogens
$mol = Chemistry::Mol->new;
$mol->new_atom(symbol => 'O', implicit_hydrogens => 2);
is( 0+$mol->atoms,      1,      'before sprout_hydrogens' );
$mol->sprout_hydrogens;
is( 0+$mol->atoms,      3,      'after sprout_hydrogens' );
$mol->collapse_hydrogens;
is( 0+$mol->atoms,      1,      'before sprout_hydrogens' );


# Bond methods
is($bond->length, 5, '$bond->length');

