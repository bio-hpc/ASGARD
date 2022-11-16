use Test::More;
use Chemistry::Mol;
use Chemistry::Ring::Find ':all';
use strict;

my @files = glob "t/rings2/*.ring";

eval "use Chemistry::File::SMILES";
if ($@) {
    plan skip_all => "You don't have Chemistry::File::SMILES installed";
} else {
    plan tests => 0 + @files;
}

for my $file (@files) {
    open F, $file or die "couldn't open $file: $!\n";   
    my ($smiles, $options, @expected_rings) = 
        map { /: ([^\n\r]*)/g } <F>;
    
    my ($mol, $patt);
    Chemistry::Atom->reset_id;
    Chemistry::Bond->reset_id;
    $mol = Chemistry::Mol->parse($smiles, format => 'smiles');
    my %opts = split " ", $options;

    my @rings    = find_rings($mol, %opts);
    #use List::Util 'shuffle'; @rings = shuffle @rings;
    @rings = map { $_->[1] } 
        sort { $a->[0] cmp $b->[0] } 
        map { [ join(" ", $_->atoms) => $_ ] } @rings;

    my @got_rings;
    for my $ring (@rings) {
        my @atoms    = $ring->atoms;
        my @bonds    = $ring->bonds;
        my $aromatic = $ring->is_aromatic;

        push @got_rings, "atoms(@atoms); bonds(@bonds); aromatic($aromatic)";
    }

    if ($opts{sssr}) {
        is(scalar @got_rings, scalar @expected_rings, "$file: $smiles");
    } else {
        is_deeply(\@got_rings, \@expected_rings, "$file: $smiles");
    }
}
