use Test::More;
use Chemistry::Mol 0.24;
use Chemistry::Pattern 0.21;
use strict;

my @files;

BEGIN { 
    eval qq{ 
        use Chemistry::Ring 0.11 'aromatize_mol';
        use Chemistry::File::SMILES 0.40;
    };
    @files = glob "t/pats/*.pat" unless $@;
    plan tests => 1 + @files * 2;
    use_ok('Chemistry::File::SMARTS');
};


for my $file (@files) {
    open F, $file or die "couldn't open $file\n";   
    my ($patt_str, $options, $mol_str, @expected_matches) = map { /: ([^\n\r]*)/g } <F>;
    
    my ($mol, $patt);
    Chemistry::Atom->reset_id;
    $patt = Chemistry::Pattern->parse($patt_str, format => "smarts");
    is ($patt->name, $patt_str, "\$patt->name($patt_str)");
    $mol = Chemistry::Mol->parse($mol_str, format => 'smiles');
    aromatize_mol($mol);
    $patt->options(split " ", $options);

    my @matches;
    while ($patt->match($mol) ) {
        my @ret = $patt->atom_map;
        push @matches, "(@ret)";
    }
    push @matches, "()";

    is_deeply(\@matches, \@expected_matches, "$file: $mol_str =~ /$patt_str/");
}
