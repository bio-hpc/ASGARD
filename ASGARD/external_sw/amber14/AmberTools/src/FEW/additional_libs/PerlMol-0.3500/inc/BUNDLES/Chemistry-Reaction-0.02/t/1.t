#!/usr/bin/perl

use strict;
use warnings;
use Chemistry::File::SMILES;
use Chemistry::Pattern;
use Chemistry::Reaction;
use Test::More tests => 16;

open F, "t/reactions.txt" or die;
my @reacts;
while (<F>) {
    chomp;
    my ($subst, $prod, $name, $map) = split " ", $_;
    my $s = Chemistry::Pattern->parse($subst, format=>'smiles');
    my $p  = Chemistry::Pattern->parse($prod, format=>'smiles');
    my %m;
    if ($map) {
        my @smap = split ",", $map;
        for (my $i = 1; $i <= $s->atoms; $i++) {
            $m{$s->atoms($smap[$i-1])} = $p->atoms($i);
        }
    } else {
        for (my $i = 1; $i <= $s->atoms; $i++) {
            $m{$s->atoms($i)} = $p->atoms($i);
        }
    }
    my $react = Chemistry::Reaction->new($s, $p, \%m, name => $name);
    push @reacts, $react;
}

open F, "t/1.out" or die;

my @mols = Chemistry::Mol->read('t/mols.smi');
for my $react (@reacts) {
    for my $mol (@mols) {
        my @products;
        my $got = sprintf "%s\t%s\t", $react->name, $mol->name;
        my $subst = $react->substrate;
        while ($subst->match($mol)) {
            my $new_mol = $mol->clone;
            my @map = map { $new_mol->by_id($_) } $subst->atom_map;
            $react->forward($new_mol, @map);
            push @products, $new_mol;
        }
        $got .= join(", ", map {$_->print(format=>'smiles')} @products);
        my $expected = <F>;
        chomp $expected;
        is($got, $expected, "$expected");
    }
}
