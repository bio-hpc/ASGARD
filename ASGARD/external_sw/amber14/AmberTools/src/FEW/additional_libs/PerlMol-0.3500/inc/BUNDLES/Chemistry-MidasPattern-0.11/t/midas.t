#!/home/ivan/bin/perl

use strict;
use warnings;
use Test::More;

use Chemistry::File::MidasPattern;

my @files = glob("t/pats/*.pat");
my $n = @files;

eval "use Chemistry::File::PDB";
if ($@) {
    plan skip_all => "You don't have Chemistry::File::PDB installed";
} else {
    plan tests => $n;
}

my $mol = Chemistry::MacroMol->read(shift || "test.pdb");

for my $fname (@files) {
    open F, "<$fname" or die "couldn't open $fname: $!\n";
    my ($str, @expected) = map { /: (.*)/g } <F>;

    my $patt = Chemistry::MidasPattern->parse($str, format => 'midas');
    $patt->match($mol);
    my @got = map { 
        sprintf "%s\t%s",  $_->attr("pdb/residue_name"), $_->name
    } $patt->atom_map;
    push @got, sprintf("%d atoms", 0+$patt->atom_map);

    is_deeply(\@got, \@expected, "$str");
    #use Data::Dumper; print Dumper \@got, \@expected;
}
