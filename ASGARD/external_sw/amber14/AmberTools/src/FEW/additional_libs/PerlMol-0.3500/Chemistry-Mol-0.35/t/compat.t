use Test::More;

plan tests => 1;
ok(1);

eval {
    require Chemistry::File::SMILES;
    my $v = Chemistry::File::SMILES->VERSION;
    if ($v < 0.43) {
        diag "You have Chemistry::File::SMILES version $v installed. It is not compatible with this version of Chemistry::Mol. Please upgrade to Chemistry::File::SMILES 0.43 or higher.\n";
    }
};
