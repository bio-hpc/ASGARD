#!/home/ivan/bin/perl

use Chemistry::File::SMILES;
use Test::More;

my $file = shift || "unique_smiles.txt";
open F, "$file" or die "couldn't open $file\n";
my @lines = <F>;
plan tests => @lines*1;

for my $line (@lines) {
    my ($smiles, $expected_usmiles) = split " ", $line;
    my $mol = Chemistry::Mol->parse($smiles, format => 'smiles');
    my $usmiles = $mol->print(format => 'smiles', unique => 1);
    is($usmiles, $expected_usmiles, "$smiles => $expected_usmiles");
}
