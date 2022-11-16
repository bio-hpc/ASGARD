use strict;
use warnings;
no warnings qw(qw);

use Test::More;
use Chemistry::3DBuilder qw(build_3d);
use Chemistry::File::SMILES;
use Chemistry::File::SMARTS;
use Chemistry::Pattern;
use Chemistry::Atom qw(angle_deg dihedral_deg distance);


my @mols = qw( 
    C CC CCC CCCC CCCCC C(CC)CC C(C)(C)(C)C 
    C=C CC=CC C#C C=C=C [H][H]
    O N
);

plan tests => scalar @mols;

my @bond_types = (
    {
        name => 'single',
        patt => Chemistry::Pattern->parse('CC', format => 'smarts'),
        val => 1.5,
    },
    {
        name => 'double',
        patt => Chemistry::Pattern->parse('*=*', format => 'smarts'),
        val => 1.35,
    },
    {
        name => 'triple',
        patt => Chemistry::Pattern->parse('*#*', format => 'smarts'),
        val => 1.2,
    },
    {
        name => 'hydrogen',
        patt => Chemistry::Pattern->parse('*[H]', format => 'smarts'),
        val => 1.1,
    },
);

my @ang_types = (
    {
        name => 'sp3',
        patt => Chemistry::Pattern->parse('*[*X4,NX3,OX2,OX3]*', 
            format => 'smarts'),
        val => 109.47,
    },
    {
        name => 'sp2',
        patt => Chemistry::Pattern->parse('*~[CX3,NX2]~*', format => 'smarts'),
        val => 120,
    },
    {
        name => 'sp',
        patt => Chemistry::Pattern->parse('*~[CX2]~*', format => 'smarts'),
        val => 180,
    },
);

my @dih_types = (
    {
        name => 'proper sp3',
        patt => Chemistry::Pattern->parse('*[CX4][CX4]*', 
            format => 'smarts'),
        mod => 120,
        min => 60,
    },
    {
        name => 'proper sp2',
        patt => Chemistry::Pattern->parse('*C=C*', 
            format => 'smarts'),
        mod => 180,
        min => 0,
    },
    {
        name => 'improper sp3',
        patt => Chemistry::Pattern->parse('*[*X4,NX3,OX3,OX2](*)*', 
            format => 'smarts'),
        mod => 120,
        min => 0,
    },
    {
        name => 'improper sp2',
        patt => Chemistry::Pattern->parse('*C(=C)*', 
            format => 'smarts'),
        mod => 360,
        min => 180,
    },
);

MOL: for my $s (@mols) {
    my $mol = Chemistry::Mol->parse($s, format => 'smiles');
    build_3d($mol);

    my $n_bond = 0;
    for my $bt (@bond_types) {
        while ($bt->{patt}->match($mol)) {
            $n_bond++;
            my $len = distance($bt->{patt}->atom_map);
            if (abs($len - $bt->{val}) > 0.01) {
                ok(0, "'$s', len = '$len'");
                next MOL;
            }
        }
    }

    my $n_ang = 0;
    for my $at (@ang_types) {
        while ($at->{patt}->match($mol)) {
            $n_ang++;
            my $ang = angle_deg($at->{patt}->atom_map);
            if (abs($ang - $at->{val}) > 0.1) {
                ok(0, "'$s', ang='$ang', pat=$at->{name}: '$at->{patt}{name}'");
                next MOL;
            }
        }
    }

    my $n_dih = 0;
    for my $dt (@dih_types) {
        while ($dt->{patt}->match($mol)) {
            $n_dih++;
            my $dih = dihedral_deg($dt->{patt}->atom_map);
            $dih = Math::Complex::Re($dih) 
                if UNIVERSAL::isa($dih, 'Math::Complex');
            if (abs(($dih - $dt->{min} + 1) % $dt->{mod}) > 2) {
                my $path = join("", map { $_->symbol } $dt->{patt}->atom_map)
                        . join("", map { $_->order } $dt->{patt}->bond_map);
                ok(0, "'$s', dih='$dih', pat=$dt->{name}: '$dt->{patt}{name}', path='$path'");
                next MOL;
            }
        }
    }
    ok(1, "$s; $n_bond bonds, $n_ang angles, $n_dih dihedrals");
}


