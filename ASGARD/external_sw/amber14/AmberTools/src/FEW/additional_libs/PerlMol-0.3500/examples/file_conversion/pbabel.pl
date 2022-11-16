#!/home/ivan/bin/perl

# This sample program is called 'pbabel' because it is inspired by the Babel
# program and it is written in Perl using PerlMol.

# USAGE
# ./pbabel.pl [options] <input file> <output file>

use strict;
use warnings;

# load all File I/O modules
use Chemistry::File ':auto';

# load helper modules for finding bonds, aromatizing, and generating coords
use Chemistry::3DBuilder 'build_3d';
use Chemistry::Ring 'aromatize_mol';
use Chemistry::Bond::Find 'find_bonds';

# parse options
use Getopt::Long;
Getopt::Long::Configure ("bundling");
my ($input_format, $output_format, $build_3d, $find_bonds, $aromatize,
    $no_kekulize, $unique, $name); 
my $result = GetOptions(
    "i=s"    => \$input_format,
    "o=s"    => \$output_format,
    "3"      => \$build_3d,
    "b"      => \$find_bonds,
    "a"      => \$aromatize,
    "K"      => \$no_kekulize,
    "u"      => \$unique,
    "n"      => \$name,
);

# check if everything is OK, or else print usage info
if (@ARGV != 2 or !$result) {
    my @formats = Chemistry::Mol->formats;
    print <<OPTIONS;
$0 [options] <input file> <output file>

Options:
    -i <input format>
    -o <output format>
    -3      build 3d coordinates 
    -b      find bonds
    -a      aromatize
    -K      don't kekulize (only for reading SMILES)
    -u      unique (only for writing SMILES)
    -n      include name (only for writing SMILES)

Available file formats (if omitted they will be guessed):
OPTIONS
print "    $_\n" for @formats;
    exit;
}

my ($input_file, $output_file) = @ARGV;

# read the input file
my @mols = Chemistry::Mol->read($input_file, 
    format   => $input_format,
    kekulize => ! $no_kekulize,     # only used by SMILES
);

# do optional procession on the molecules
for my $mol (@mols) {
    build_3d($mol)                  if $build_3d;
    find_bonds($mol, orders => 1)   if $find_bonds;
    aromatize_mol($mol)             if $aromatize;
}

# write the output file
$mols[0]->write($output_file, 
    format => $output_format, 
    mols   => \@mols,           # only used by multi-molecule files 
                                # such as SMILES and SDF
    name   => $name,            # only used by SMILES
    unique => $unique,          # only used by SMILES
);

