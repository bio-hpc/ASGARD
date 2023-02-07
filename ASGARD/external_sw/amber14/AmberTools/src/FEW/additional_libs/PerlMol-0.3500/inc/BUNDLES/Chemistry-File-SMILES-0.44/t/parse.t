use Test;

our @files;
BEGIN {
    @files = glob "t/*.sm";
    plan tests => 1 + @files; 
};

use Chemistry::File::SMILES;
ok(1); # If we made it this far, we're ok.

no warnings 'uninitialized';
my $parser = Chemistry::File::SMILES->new_parser(
    add_atom => sub {
        my $c=shift; 
        pop; # to get rid of the new name property (this is a very old test)
        local $"=',';
        $c->{out} .= "ATOM$c->{i}(@_)\n"; 
        $c->{i}++;
    },
    add_bond => sub {
        my $c=shift; 
        local $"=',';
        $c->{out} .= "BOND(@_)\n";
    }
);

for $fname (@files) {
    open F, $fname or die;
    my $content;
    { local undef $/; $content = <F>; }
    my ($s) = $content =~ /(.*)$/m;
    my $c = {i=>1, out=>"$s\n"};
    eval {$parser->parse($s, $c);};
    ok($c->{out}, $content);
}

