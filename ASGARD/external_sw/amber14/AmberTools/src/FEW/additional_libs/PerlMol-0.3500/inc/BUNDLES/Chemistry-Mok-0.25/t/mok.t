use strict;
use Test::More;

use Chemistry::Mok;

eval "use Chemistry::File::MDLMol";
if ($@) {
    plan skip_all => "Chemistry::File::MDLMol is needed for testing";
} else {
    my @fnames = glob('*.mok');
    plan tests => @fnames * 2;
    
    for my $fname (@fnames) {
        $fname =~ s/\.mok$//;
        open CODE, "<", "$fname.mok" or die "couldn't open $fname.mok: $!\n";
        my $code = join "", <CODE>;

        open EXPECTED, "<", "$fname.out" 
            or die "couldn't open $fname.out: $!\n";
        my @expected = <EXPECTED>;

        open GOT, '-|', "$^X -Mblib mok -f $fname.mok *.mol" 
            or die "couln't open mok pipe: $!\n";
        my @got = <GOT>;
        is_deeply(\@got, \@expected, "mok -f $fname.mok *.mol");

        SKIP: {
            if ($^O eq 'MSWin32') {
                $code =~ s/\n//g;
                skip "MSWin32", 1 if $code =~ /#/;
                open GOT, '-|', "$^X -Mblib mok \"$code\" *.mol"
                    or die "couln't open mok pipe: $!\n";
            } else {
                open GOT, '-|', "$^X -Mblib mok '$code' *.mol" 
                    or die "couln't open mok pipe: $!\n";
            }
            @got = <GOT>;
            is_deeply(\@got, \@expected, "mok '<$fname.mok>' *.mol");
        }
    }
}

