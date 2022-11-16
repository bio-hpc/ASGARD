#!/home/ivan/bin/perl

use Chemistry::File::SMILES;
use Test::More 'no_plan';

my $file = shift || "test_names.txt";
my @mols = Chemistry::Mol->read($file, format => 'smiles');
open F, "<$file" or die "couln't open $file: $!\n";
my $expected = do { local $/; <F> };

$out = Chemistry::Mol->print(format => 'smiles', mols => \@mols, name => 1,
    aromatic => 1);
is($out, $expected, "multiline $file");
