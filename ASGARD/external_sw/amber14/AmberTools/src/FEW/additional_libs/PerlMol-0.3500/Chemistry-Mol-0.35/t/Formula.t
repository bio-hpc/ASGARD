use Test::More;
my @lines;

BEGIN { 
    open F, "<", "formula_tests.txt"
        or die "couldn't open formula_tests.txt; $!";
    @lines = <F>;
    close F;
    plan tests => 5 + @lines;
    use_ok('Chemistry::File::Formula');
};

# Constructors
my $mol = Chemistry::Mol->parse("CH4O", format => "formula");
isa_ok($mol, 'Chemistry::Mol', 'parse isa mol');
ok($mol->atoms == 6, "enough atoms");
my $formula = $mol->formula("%s%d{<sub>%d</sub>}");
is($formula, "CH<sub>4</sub>O", "formula format");

$mol = Chemistry::Mol->parse("1[Ph(Me)3]2", format => "formula");
my $fh = $mol->formula_hash;
is_deeply($fh, {C => 18, H => 28}, "formula hash 1[Ph(Me)3]2");

for my $line (@lines) {
    chomp $line;
    my ($test_formula, $expected) = split /\t/, $line;
    my $got = Chemistry::Mol->parse($test_formula, format => "formula")
        ->print(format=>'formula');
    is($got, $expected, "$test_formula = $expected");
}
