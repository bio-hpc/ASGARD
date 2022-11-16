#!/home/ivan/bin/perl

# usage: ./molgrep.pl <SMARTS pattern> <file(s)...>

# prints to STDOUT a table with the filenames, names, formulas, and SMILES
# strings for the matching molecules in the files given in the command line

use strict;
use warnings;

# use all available file types
use Chemistry::File ':auto';

# for molecules with no explicit bonds
use Chemistry::Bond::Find 'find_bonds';

# so that ring and aromatic SMARTS work
use Chemistry::Ring 'aromatize_mol';

# make sure we have enough command-line arguments
unless (@ARGV >= 2) {
    print "$0 <SMARTS pattern> <file(s)...>\n";
    exit;
}

# create the pattern object
my $patt_str = shift @ARGV;
my $patt = Chemistry::Pattern->parse($patt_str, format => 'smarts');

# print the header
print "Filename\tName\tFormula\tSMILES\n";

# loop over the filenames
for my $filename (@ARGV) {

    # read the file. Note that there can be multiple molecules in one file
    my @mols = Chemistry::Mol->read($filename);

    # loop over the molecules from this file
    for my $mol (@mols) {

        # if the molecule doesn't have any bonds, because it came from an XYZ
        # file or something like that, find them out from the 3d coordinates.
        find_bonds($mol, orders => 1) unless $mol->bonds;

        # so that ring and aromatic SMARTS work
        aromatize_mol($mol);

        # check if molecule matches
        if ($patt->match($mol)) {
            no warnings qw(uninitialized); # in case $mol has no name

            # print the molecular information. In $mol->sprintf,
            # %n is the name, %f the formula, and %S the canonical SMILES
            print "$filename\t" . $mol->sprintf("%n\t%f\t%S\n");
        }
    }
}

