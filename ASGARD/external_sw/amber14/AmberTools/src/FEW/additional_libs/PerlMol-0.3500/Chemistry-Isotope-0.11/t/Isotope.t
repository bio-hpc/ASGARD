use Test::More;
use Data::Dumper;

#plan 'no_plan';
plan tests => 6;

my $TOL = 0.00001;

use Chemistry::Isotope qw(:all);

# Mass tests

my $m = isotope_mass(235, 92);
is_float($m, 235.043923094753, "235U");

$m = isotope_mass(235, 6);
ok(! defined $m, "235C");

$m = isotope_mass(400, 180);
ok(! defined $m, "400, 180");


# Abundance tests

my $ab = isotope_abundance('C');
is_float($ab->{12}, 98.93, "12C abundance");
is_float($ab->{13}, 1.07, "13C abundance");

$ab = isotope_abundance('Xyz');
is ( scalar keys %$ab, 0,   'Xyz abundance');

#################

sub is_float {
    my ($got, $expected, $name) = @_;
    ok( abs($got-$expected) < $TOL, $name);
}
