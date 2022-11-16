use Test::More;
use Chemistry::File::SLN;

# this test is to make sure that the bug reported at 
# https://sourceforge.net/tracker/?func=detail&atid=646299&aid=1172562&group_id=106995
# doesn't happen again

plan tests => 2;

my $mol = Chemistry::Mol->parse('OCC', format => 'sln');
$mol->atoms(3)->formal_charge(-1);

is ($mol->print(format => 'sln'), 'OCC[-1]',  'once');
is ($mol->print(format => 'sln'), 'OCC[-1]',  'twice');
