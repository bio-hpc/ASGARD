use strict;
use warnings;

#use Test::More "no_plan";
use Test::More tests => 19;

BEGIN { 
    use_ok('Chemistry::Obj');
};

my ($obj, $obj2, $obj3);

# constructor
$obj = Chemistry::Obj->new;
isa_ok( $obj, 'Chemistry::Obj',  'blank obj' );

# obj using standard attributes
$obj = Chemistry::Obj->new(name => 'joe', id => 'joe01', type => 'funny');
is( $obj->name, 'joe',      'name' );
is( $obj->id,   'joe01',    'joe01' );
is( $obj->type, 'funny',    'type' );
is( "$obj",     'joe01',    'stringify' );

# relational operators
$obj2 = Chemistry::Obj->new(name => 'joe', id => 'joe01', type => 'funny');
$obj3 = $obj;
is( $obj,   $obj2,      'eq' );
ok( $obj != $obj2,      '!=' );
ok( $obj == $obj3,      '==' );

$obj2->id('joe02');
ok( $obj ne $obj2,      'ne' );

# user attributes
$obj->attr(color => 123);
is ($obj->attr('color'),    123,    'attr');

# attr(list)
$obj->attr(1,2,3,4);
is ($obj->attr(3),          4,      'attr list');
is ($obj->attr('color'),    123,    'attr list');

# attr(hashref)
$obj->attr({ a => 1, b => 2 });
is ($obj->attr('a'),        1,      'attr hashref');
is ($obj->attr('color'),    undef,  'attr hashref');

# attr()
my $attr = $obj->attr;
is ($attr->{b},        2,      'attr get hashref');

# del_attr
$obj->del_attr('b');
$attr = $obj->attr;
ok( ! exists $attr->{b},    'del_attr');

# accessor
package _test;
our @ISA = ('Chemistry::Obj');
Chemistry::Obj::accessor(qw(a b c d));

package main;
$obj = _test->new(a => 1, b => 2);
is ($obj->a, 1,   'accessor');

$obj->a(3)->c(4);
is ($obj->c, 4,   'chained accessor');


