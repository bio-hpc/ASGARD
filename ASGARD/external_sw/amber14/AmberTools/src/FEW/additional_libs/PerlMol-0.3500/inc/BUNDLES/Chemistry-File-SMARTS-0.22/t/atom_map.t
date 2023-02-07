use Test::More;
use Chemistry::File::SMARTS;

plan tests => 3;

my $s = Chemistry::Pattern->parse('C(=[O:1])[O:2]C', format => 'smarts');
is ($s->atoms(1)->name, undef);
is ($s->atoms(2)->name, 1);
is ($s->atoms(3)->name, 2);
