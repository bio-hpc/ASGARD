use Test::More tests => 2;

use strict;
use warnings;
use Chemistry::File::XYZ;

my $mol = Chemistry::Mol->read("cyclobutene.xyz");

# default format
my $got = $mol->print(format => "xyz");
open F, "cyclobutene2.xyz" or die "Couldn't open cyclobutene2.xyz: $!\n";
my $expected = do { local $/; <F> };

is($got, $expected, "cyclobutene.xyz (default)");

# custom format
$got = $mol->print(format => "xyz", symbol => 0, count => 0, name => 0);
open F, "cyclobutene3.xyz" or die "Couldn't open cyclobutene3.xyz: $!\n";
$expected = do { local $/; <F> };

is($got, $expected, "cyclobutene.xyz (custom)");
