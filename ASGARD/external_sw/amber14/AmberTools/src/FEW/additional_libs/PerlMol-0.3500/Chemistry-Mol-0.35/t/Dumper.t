use Test::More tests => 6;

BEGIN { use_ok('Chemistry::File::Dumper') };

my $mol = Chemistry::Mol->read("t/mol.pl", format => 'dumper');
isa_ok($mol, "Chemistry::Mol", 'read');

my $mol_auto = Chemistry::Mol->read("t/mol.pl");
isa_ok($mol, "Chemistry::Mol", 'read (autodetect)');

my $s = $mol->print(format=>'dumper');
my $s2 = $mol_auto->print(format=>'dumper');
is($s, $s2, 'dump and compare');

$mol->write('t/test1.pl');
my $mol3 = Chemistry::Mol->read("t/test1.pl");
my $s3 = $mol->print(format=>'dumper');

is($s, $s3, 'write and read and compare');
is_deeply($mol, $mol3, 'deep compare');
unlink 't/test1.pl';

