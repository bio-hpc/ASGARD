use strict;
use warnings;

use Test::More;
use Chemistry::Mol;
use Chemistry::File::Dumper;

if (eval 'use Compress::Zlib; 1') {
    plan tests => 6;
    #plan 'no_plan';
} else {
    plan skip_all => "You don't have Compress::Zlib installed";
}

my $mol;
my $in  = 't/mol.pl.gz';
my $out = 't/tmp/mol.pl.gz';
my $na  = 8; # expected number of atoms

####### read tests

$mol = Chemistry::Mol->read($in, gzip => 1, format => 'dumper');
isa_ok($mol, "Chemistry::Mol", "explicit decompressed read");
ok($mol->atoms == $na, "has $na atoms");

$mol = Chemistry::Mol->read($in, format => 'dumper');
isa_ok($mol, "Chemistry::Mol", "implicit decompressed read");
ok($mol->atoms == $na, "has $na atoms");


####### write tests

mkdir 't/tmp';

$mol->write($out, format => 'dumper', gzip => 1);
is_gzipped($out);
unlink $out;

$mol->write($out, format => 'dumper');
is_gzipped($out, "implicit compression on output");
unlink $out;

rmdir 't/tmp';

sub is_gzipped {
    my ($fname, $comment) = @_;
    my $header;

    open F, $fname or die;
    read F, $header, 2;
    ok($header eq "\x1f\x8b", $comment || "compressed ok");
    close F;
}


