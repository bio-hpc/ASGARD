#!/usr/bin/perl

# Combinatorial library enumeration example.
# Sent by Konrad Koehler, Karo Bio AB

# Read in two SMILES files, one containing amines and the second acyl chlorides
# and print out SMILES strings of the product amides.
# Input files: <amines.smi> and <acid_chlorides.smi>
# Output: STDOUT

use strict;
use warnings;
use Chemistry::File::SMILES;
use Chemistry::File::SMARTS;
use Chemistry::Ring 'aromatize_mol';

#  read in two files of reagents
my @react1 = Chemistry::Mol->read('amines.smi');
my @react2 = Chemistry::Mol->read('acid_chlorides.smi');

# define SMARTS patterns for reactive functional groups

# 1. amine (but not tertiary) or aromatic nitrogen with attached hydrogen atom
my $fg1     = '[N&!D3,n&h1,n&H1]';
my $fg1_pat = Chemistry::Pattern->parse("$fg1", format => 'smarts');

# 2. acid chloride
my $fg2     = 'C(Cl)=O';
my $fg2_pat = Chemistry::Pattern->parse("$fg2", format => 'smarts');

# loop over both sets of reactants

my @prod;
for (my $i = 0 ; $i <= $#react1 ; $i++) {
    for (my $j = 0 ; $j <= $#react2 ; $j++) {

        # the next two lines are needed whenever SMARTS patterns depend on
        # aromatic/nonaromatic properties
        aromatize_mol($react1[$i]);
        aromatize_mol($react2[$j]);

        # combine two reagents to form a new disconected molecule
        my $name = $react1[$i]->name . "+" . $react2[$j]->name;
        $prod[$i][$j] = Chemistry::Mol->new(name => $name);
        $prod[$i][$j]->combine($react1[$i], $react2[$j]);

        # find reactive functional groups
        $fg1_pat->match($prod[$i][$j]);
        my @nfg1 = $fg1_pat->atom_map;
        $fg2_pat->match($prod[$i][$j]);
        my @nfg2 = $fg2_pat->atom_map;

        # test to make sure both reactants have required functional groups
        # if not, skip
        if ($nfg1[0] && $nfg2[0] && $nfg2[1]) {

            # delete displaced atoms (i.e., H and Cl) to create two radicals
            my $h_count = $nfg1[0]->implicit_hydrogens();
            $h_count = $h_count - 1;
            $nfg1[0]->implicit_hydrogens($h_count);
            $prod[$i][$j]->delete_atom($nfg2[1]);

            # add a bond between the two radicals to form the product
            $prod[$i][$j]
              ->new_bond(atoms => [$nfg1[0], $nfg2[0]], order => '1');
            my $smi =
              $prod[$i][$j]->print(format => 'smiles', unique => 1, name => 1);

            #  print out the product!
            print "$smi\n";
        }
    }
}
