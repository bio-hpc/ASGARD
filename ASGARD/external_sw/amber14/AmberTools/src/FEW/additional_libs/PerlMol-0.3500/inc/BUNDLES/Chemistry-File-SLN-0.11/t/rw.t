#!/home/ivan/bin/perl

use Test::More;
use Chemistry::File::SLN;

my $fname = "t/rw.sln";
open F, "<$fname" or die "couln't open $fname: $!\n";

my @lines = <F>;

plan tests => 0+@lines;

for my $line (@lines) {
    chomp $line;
    my $mol = Chemistry::Mol->parse($line, format => 'sln');

    my $got = $mol ? $mol->print(format => 'sln', attr => 1, name=>1) 
        : "ERROR: $line";
    is($got, $line, $line);
}
