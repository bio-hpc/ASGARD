use strict;
use warnings;

use Test::More;

BEGIN { 
    #plan 'no_plan';
    plan tests => 15;
    use_ok('Chemistry::File');
}


# simple constructor test
my $f = Chemistry::File->new;
isa_ok($f, "Chemistry::File");

require Chemistry::File::Dumper;

# file reader test
my $fname = 't/mol.pl';
my $file = Chemistry::File::Dumper->new(file => $fname);
isa_ok($file, "Chemistry::File::Dumper");
my $mol = $file->read(format => 'dumper');

isa_ok($mol, "Chemistry::Mol", 'read file');
is(scalar $mol->atoms, 8, "atoms");


# string reader test
open F, "<$fname" or die;
my $s = do {local $/; <F>};
$file = Chemistry::File::Dumper->new(file => \$s);
$mol  = $file->read;

isa_ok($mol, "Chemistry::Mol", 'read string');
is(scalar $mol->atoms, 8, "atoms");


# subclass test

package MolList;
use base qw(Chemistry::File);

sub read_header {
    my ($self) = @_;
    my $fh = $self->fh;
    my $name = <$fh>;
    chomp $name;
    $self->name($name);
}

sub read_footer {
    my ($self) = @_;
    my $fh = $self->fh;
    my $footer = join '', <$fh>;
    chomp $footer;
    $self->attr('list/footer', $footer);
}

sub slurp_mol {
    my ($self) = @_;
    my $fh = $self->fh;
    my $s = <$fh>;
    chomp $s;
    return if $s eq '--END--';
    $s;
}

sub parse_string {
    my ($self, $s) = @_;
    my $fh = $self->fh;
    my ($name, $price) = split "\t", $s;
    my $mol = Chemistry::Mol->new(name => $name);
    $mol->attr('list/price', $price);
    push @{$self->{mol_list}}, $mol;
    $mol;
}

package main;

my $list = MolList->new(file => 't/list.txt');
my @mols = $list->read;

is(scalar @mols, 3, "read mollist");
is($list->name, "This is a list of molecules", "list name");
is($list->attr('list/footer'), 
    "As you can see, \neverything you need \nis on this list.", "list footer");
is($mols[1]->name, "ethane", "mol name");
is($mols[1]->attr('list/price'), '$200', "mol price");


# Chemistry::Mol->file 
$file = Chemistry::Mol->file($fname);
isa_ok($file, "Chemistry::File::Dumper");
$mol = $file->read(format => 'dumper');

isa_ok($mol, "Chemistry::Mol", 'Chemistry::Mol->read');
is(scalar $mol->atoms, 8, "atoms");

