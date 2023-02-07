use strict;
use warnings;

use Chemistry::File::SMILES;
use Test::More;
use List::Util qw(sum);

#plan tests => 8;
plan 'no_plan';

my $mol;
my $total_bo;
my $s_out;
my $s_in = 'c1ccccc1C';
my $bo_total;

# kekulize
$mol = Chemistry::Mol->parse($s_in, format => 'smiles', kekulize => 0);
$bo_total = sum( map {$_->order} $mol->bonds );
is( $bo_total,      7,      'kekulize => 0' );

$mol = Chemistry::Mol->parse($s_in, format => 'smiles', kekulize => 1);
$bo_total = sum( map {$_->order} $mol->bonds );
is( $bo_total,      10,      'kekulize => 1' );

# aromatic
$s_out = $mol->print(format => 'smiles', aromatic => 0);
is( $s_out,     'C1=CC=CC=C1C',    'aromatic => 0');

$s_out = $mol->print(format => 'smiles', aromatic => 1);
is( $s_out,     $s_in,    'aromatic => 1');

# unique
$s_out = $mol->print(format => 'smiles', unique => 1);
is( $s_out,     'Cc1ccccc1',    'unique => 1');

# auto_number
$s_out = $mol->print(format => 'smiles', aromatic => 1, auto_number => 1);
is( $s_out, '[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[CH3:7]',  'auto_number => 1');

# number
$mol->atoms(7)->name(345);
$s_out = $mol->print(format => 'smiles', aromatic => 1, number => 1);
is( $s_out, 'c1ccccc1[CH3:345]',  'number => 1');

# read numbers
$s_in = 'CC(=[O:1])[O:2]CC';
$mol = Chemistry::Mol->parse($s_in, format => 'smiles');
is( $mol->atoms(3)->name, 1, 'atom name(1)' );
is( $mol->atoms(4)->name, 2, 'atom name(2)' );

