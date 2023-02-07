use Test::More;

# Tests for miscelaneous methods that play around with the molecular graph
# as a whole, such as clone, combine, and separate

use strict;
use warnings;

#plan 'no_plan';
plan tests => 11;

use Chemistry::File::Dumper;

my $mol = Chemistry::Mol->read("t/mol.pl");
$mol->bonds(7)->delete;

# clone test
my $mol2 = $mol->clone;
is_deeply( $mol, $mol2, "clone" );

# separate test
my @mols = $mol->separate;
is ( scalar @mols, 2, 'got 2 things' );
is ( scalar (grep $_->isa('Chemistry::Mol'), @mols), 2, 'separate: two mols' );
is ( $mols[0]->formula, 'CClH2',    'mol 1 is CClH2' );
is ( $mols[1]->formula, 'CHO2',     'mol 2 is CHO2' );
my $a1 = $mol->atoms(2);
my $a2 = $mols[0]->atoms(2);
my $nb_before = $a1->neighbors;
my $nb_after  = $a2->neighbors;
is ( $nb_after, $nb_before, "bond count for $a2 equal to $a1 ($nb_before)?" );

# combine - new
my $comb_new = Chemistry::Mol->combine(@mols);
isa_ok($comb_new, 'Chemistry::Mol');
for my $method (qw(atoms bonds formula)) {
    is ( scalar $comb_new->$method, scalar $mol->$method, "combine-new; $method" );
}

# combine - in place
my $comb_inplace = $mols[0]->combine($mols[1]);
is_deeply ( $comb_inplace, $mol, "combine-in place" );

#use Chemistry::File::SMILES; $mol->printf("%S\n");
#$_->printf("%f\n") for @mols;

