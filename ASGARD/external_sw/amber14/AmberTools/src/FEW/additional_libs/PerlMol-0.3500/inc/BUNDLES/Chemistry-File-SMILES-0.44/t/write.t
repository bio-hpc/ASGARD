#!/home/ivan/bin/perl

use Chemistry::File::SMILES;
use Test::More;

my $file = shift || "test_smiles.txt";
open F, "$file" or die "couldn't open $file\n";
my @smiles = <F>;
plan tests => scalar @smiles;

for my $smiles (@smiles) {
    chomp $smiles;
    my $mol = Chemistry::Mol->parse($smiles, format => 'smiles');
    my $out = $mol->print(format => 'smiles', aromatic => 1);
    is($out, $smiles, "$smiles");
}
