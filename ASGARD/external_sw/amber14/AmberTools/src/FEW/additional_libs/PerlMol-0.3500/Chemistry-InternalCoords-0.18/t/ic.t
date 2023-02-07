use Test::More 'no_plan';

use Chemistry::Mol;
use strict;
use warnings;

my @files;
my $TOL = 0.0001;

BEGIN { 
    @files = glob "*.ic";
    use_ok('Chemistry::InternalCoords');
}

for my $in_file (@files) {
    my $out_file = $in_file;
    $out_file =~ s/ic$/out/;
    open IN, "<", $in_file or die "couldn't open $out_file: $!\n";
    open OUT, "<", $out_file or die "couldn't open $out_file: $!\n";
    my @rows_in = map { [split] } <IN>;
    my @rows_out = map { [split] } <OUT>;

    my $mol = Chemistry::Mol->new;
    my $n = 0;
    for my $row_in (@rows_in) {
        my $row_out = shift @rows_out;
        my $atom = $mol->new_atom;
        my $ic = Chemistry::InternalCoords->new($atom, @$row_in);
        my @calc_coords = $ic->add_cartesians->array;
        $n++;
        for my $axis (qw(x y z)) {
            my $got = shift @calc_coords;
            my $expected = shift @$row_out;
            ok(abs($got-$expected) < $TOL, 
                "$in_file: $axis($n); expected $expected, got $got");
        }
    }
}
