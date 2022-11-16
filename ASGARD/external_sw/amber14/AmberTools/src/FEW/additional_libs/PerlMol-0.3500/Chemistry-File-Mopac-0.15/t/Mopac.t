use Test::More 'no_plan';

my @files;
my $TOL = 0.0001;

BEGIN { 
    @files = glob "*.mop";
    #plan tests => 1 + @files;
    use_ok('Chemistry::File::Mopac');
}

for my $mop_file (@files) {
    my $out_file = $mop_file;
    $out_file =~ s/mop$/out/;
    open F, "<", $out_file or die "couldn't open $out_file: $!\n";
    my $mol = Chemistry::Mol->read($mop_file, format => 'mop');
    my @rows = map { [split] } <F>;
    for my $row (@rows) {
        my $n = shift @$row;
        my $atom = $mol->atoms($n);
        my $symbol = shift @$row;
        my @calc_coords = $atom->coords->array;
        is($atom->symbol, $symbol, "$mop_file: symbol($n)");
        for my $axis (qw(x y z)) {
            my $got = shift @calc_coords;
            my $expected = shift @$row;
            ok(abs($got-$expected) < $TOL, 
                "$mop_file: $axis($n); got $got, expected $expected");
        }
    }
}
