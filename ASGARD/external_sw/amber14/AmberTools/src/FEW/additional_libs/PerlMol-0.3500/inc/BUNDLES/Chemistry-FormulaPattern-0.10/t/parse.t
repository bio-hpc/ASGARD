#!/home/ivan/bin/perl

use strict;
use warnings;
use Test::More;
use Chemistry::FormulaPattern;

my @files = glob "t/parse/*.pl";

plan tests => 2 * @files;

# have to be global to get into the "do" below
our ($formula_patt, $expected_ranges, $expected_options);

for my $file (@files) {
    do $file;
    my $patt = Chemistry::FormulaPattern->new($formula_patt);
    my $got_ranges = $patt->{formula_pattern};
    my $got_options = $patt->{options};
    is_deeply($got_ranges, $expected_ranges, "$formula_patt ranges");
    is_deeply($got_options, $expected_options, "$formula_patt options");
}

