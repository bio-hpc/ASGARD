#!/home/ivan/bin/perl

use strict;
use warnings;
use Test::More;
use Chemistry::File::FormulaPattern;
use Chemistry::File::Formula;

my @lines = <DATA>;

plan tests => 2*@lines;

for my $line (@lines) {
    chomp $line;
    my ($patt_str, $formula_str, $expected_match) = split " ", $line;
    my $patt = Chemistry::Pattern->parse($patt_str, 
        format => 'formula_pattern');
    my $mol =  Chemistry::Mol->parse($formula_str, format=>'formula');

    my $got_match = $patt->match($mol);
    is($got_match, $expected_match, "$line");
    $got_match = $patt->match($mol);
    is($got_match, 0, "$line");
}

__DATA__
CH4         CH4     1
CH4         CH3     0
C*          CH4     1
CH3-4       CH3     1
CH3-4       CH4     1
CH3-4       CH2     0
CH3-        CH4     1
CH3-        CH2     0
