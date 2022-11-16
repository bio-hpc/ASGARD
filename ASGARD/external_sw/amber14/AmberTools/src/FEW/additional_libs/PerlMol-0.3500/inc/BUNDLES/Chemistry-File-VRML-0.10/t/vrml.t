use strict;
use warnings;
use Test::More;
use Chemistry::File::VRML;
use Chemistry::File::Dumper;

plan tests => 1;

my $mol = Chemistry::Mol->read('test.dump', format => 'dumper');
my $got = $mol->print(
    format => 'vrml',
    center => 1,
    style  => 'ballAndWire',
    color  => 'byAtom',
);

open F, "< test.wrl" or die;
my $expected = do { local $/; <F> };

ok ($got eq $expected, "test.dump -> test.wrl");
