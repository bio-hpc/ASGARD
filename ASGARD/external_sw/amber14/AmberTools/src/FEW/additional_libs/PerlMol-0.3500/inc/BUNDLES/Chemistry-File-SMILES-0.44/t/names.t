#!/home/ivan/bin/perl

use Chemistry::File::SMILES;
use Test::More;

my $file = shift || "test_names.txt";
open F, "$file" or die "couldn't open $file\n";
my @smiles = <F>;
plan tests => @smiles * 2;

for my $smiles (@smiles) {
    chomp $smiles;
    my $mol = Chemistry::Mol->parse($smiles, format => 'smiles');
    my $out = $mol->print(format => 'smiles', aromatic => 1) . "\t" . $mol->name;
    is($out, $smiles, "$smiles (name => 0)");
    $out = $mol->print(format => 'smiles', name => 1, aromatic => 1);
    is($out, $smiles, "$smiles (name => 1)");
}
