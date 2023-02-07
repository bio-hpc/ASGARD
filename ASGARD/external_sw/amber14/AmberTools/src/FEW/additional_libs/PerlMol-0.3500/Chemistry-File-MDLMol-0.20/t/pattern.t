use strict;
use warnings;
use Test::More;
use Chemistry::File::MDLMol;

if (eval "use Chemistry::Pattern; use Chemistry::Ring; 1") {
    plan tests => 5;
    #plan 'no_plan';
    ok(1, "loaded Chemistry::Pattern");
} else {
    plan skip_all => 'Chemistry::Pattern or Chemistry::Ring not installed';
}

#$Chemistry::File::MDLMol::DEBUG = 1;
my $patt = Chemistry::Pattern->read('t/query.mol');
isa_ok ($patt, 'Chemistry::Pattern');

my @files = (
    { name => 't/mol1.mol', matches => 2, },
    { name => 't/mol2.mol', matches => 2, },
    { name => 't/mol3.mol', matches => 0, },
);

#$Chemistry::Pattern::DEBUG=1;

for my $file (@files) {
    my $mol = Chemistry::Mol->read($file->{name});
    Chemistry::Ring::aromatize_mol($mol);
    my $n = 0;
    #print "\n\n\n\n******************\nTRYING $file->{name}\n\n";
    $n++ while $patt->match($mol);
    is ($n, $file->{matches}, "$file->{name} matches $file->{matches} times?");
}

