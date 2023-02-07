#!/home/ivan/bin/perl

use strict;
use warnings;
#use Test::More 'no_plan';
#use Test::More tests => 114;
use Test::More;

if ($ENV{TEST_BUILDER}) {
    plan tests => 114;
} else {
    plan skip_all => "to run builder tests, please define \$ENV{TEST_BUILDER}";
}

use Chemistry::File::Dumper;
use Chemistry::InternalCoords::Builder ':all';

my $TOL = 0.002;

my $mol = Chemistry::Mol->read('t/builder.pl');
build_zmat($mol);
check($mol, "t/builder.out");

$mol = Chemistry::Mol->read('t/builder.pl');
build_zmat($mol, bfs => 0);
check($mol, "t/builder_nobfs.out");

sub check {
    my ($mol, $fname) = @_;
    my @result;
    push @result, [
        $_->symbol, 
        idx_val($_->internal_coords->distance),
        idx_val($_->internal_coords->angle),
        idx_val($_->internal_coords->dihedral),
    ] for $mol->atoms;


    open F, "<$fname" or die "couldn't open $fname: $!\n";
    my @expected = map { [ split(" ") ] } <F>;

    is(scalar @result, scalar @expected, "same size");

    for (my $i = 0; $i < @expected; $i++) {
        my $got = $result[$i];
        my $exp = $expected[$i];
        is($exp->[0], $got->[0], "symbol $i");
        is($exp->[1], $got->[1], "distance ref $i");
        is($exp->[3], $got->[3], "angle ref $i");
        is($exp->[5], $got->[5], "dihedral ref $i");
        is_float($exp->[2], $got->[2], "distance val $i");
        is_float($exp->[4], $got->[4], "angle val $i");
        is_float($exp->[6], $got->[6], "dihedral val $i");
    }
}


sub idx_val {
    my ($atom, $val) = @_;
    my $idx = $atom ? $atom->attr("zmat/index") || 0 : 0;
    ($idx, $val || 0);
}

sub is_float {
    my ($got, $exp, $name) = @_;
    ok(abs($got-$exp) < $TOL, $name);
}
