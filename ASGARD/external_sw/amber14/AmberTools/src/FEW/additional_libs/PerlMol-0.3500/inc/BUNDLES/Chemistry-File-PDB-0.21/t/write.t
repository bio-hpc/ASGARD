#!/home/ivan/bin/perl

use strict;
use warnings;
use Test::More tests => 1;

use Chemistry::File::PDB;

my $mol = Chemistry::MacroMol->read("test.pdb");
$mol->write("out.pdb");

my ($got, $expected);
{
    local $/;
    open F, "<test.pdb" or die "couldn't open test.pdb: $!\n";
    $expected = <F>;
    open F, "<out.pdb" or die "couldn't open test.pdb: $!\n";
    $got = <F>;
}

is($got, $expected, "pdb write test");
