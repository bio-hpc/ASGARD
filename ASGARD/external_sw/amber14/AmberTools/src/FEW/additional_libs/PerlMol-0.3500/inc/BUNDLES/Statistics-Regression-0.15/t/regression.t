use strict;
#use warnings; # not available before 5.006

use lib qw(blib/lib inc);
use Test::More tests => 5;

my $TOLERANCE = 1e-6;

BEGIN {
    # make sure module loads correctly
    use_ok('Statistics::Regression'); 
}

# do the regression example
my $reg= Statistics::Regression->new( 3, "sample regression", [ "const", "someX", "someY" ] );
$reg->include( 2.0, [ 1.0, 3.0, -1.0 ] );
$reg->include( 1.0, [ 1.0, 5.0, 2.0 ] );
$reg->include( 20.0, [ 1.0, 31.0, 0.0 ] );
$reg->include( 15.0, [ 1.0, 11.0, 2.0 ] );

my @theta = $reg->theta;
my $rsq   = $reg->rsq;


# expected output
my @expected_theta = (0.295033929673043, 0.672270203578038, 1.06878470080197);
my $expected_rsq   = 0.808185547954473;

# Compare expected with actual output
for my $i (0 .. 2) {
    ok (abs($expected_theta[$i] - $theta[$i]) < $TOLERANCE, "theta[$i]")
        or diag("expected $expected_theta[$i]; got $theta[$i]");
}

ok (abs($expected_rsq - $rsq) < $TOLERANCE, "rsq")
    or diag("expected $expected_rsq; got $rsq");

