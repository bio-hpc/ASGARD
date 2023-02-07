use Test::More;
my @files = (glob("lib/Chemistry/*.pm lib/Chemistry/Ring/*.pm"));
my $n = @files;

eval 'use Test::Pod';
if ($@) {
    plan skip_all => "You don't have Test::Pod installed";
} else {
    plan tests => $n;

}

for my $file (@files) {
    pod_file_ok($file, "POD for '$file'");
}
