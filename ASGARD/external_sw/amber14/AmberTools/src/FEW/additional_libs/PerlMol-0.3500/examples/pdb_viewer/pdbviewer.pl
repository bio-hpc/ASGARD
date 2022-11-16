#!/usr/bin/env perl

# Adapted from an original script by Simon Cozens.

use strict;
use warnings;
use OpenGL::Simple qw(:all);
use OpenGL::Simple::GLUT qw(:all);
use OpenGL::Simple::Viewer;
use Chemistry::File::PDB;
use Chemistry::Bond::Find ':all';

my ($mass_scale,$sphericity) = (10,20);


my %colour = (
    red     => [ 1,   0,   0,   1 ],
    yellow  => [ 1,   1,   0,   1 ],
    orange  => [ 1,   0.5, 0,   1 ],
    green   => [ 0,   1,   0,   1 ],
    cyan    => [ 0,   1,   1,   1 ],
    blue    => [ 0,   0,   1,   1 ],
    magenta => [ 1,   0,   1,   1 ],
    grey    => [ 0.5, 0.5, 0.5, 1 ],
    white   => [ 1,   1,   1,   1 ],
);

my %ccache;
my @colours = values %colour;
my $iter = 0;

my %element_colours = (
    C => $colour{grey},
    O => $colour{red},
    N => $colour{blue},
    H => $colour{white},
    S => $colour{yellow},
    P => $colour{orange},

);


my @ballpoints = ();
my @ballsticks = ();

$|=1;
my $filename = shift or die("Need PDB filename");
print "Reading molecule..";
my $mol = Chemistry::MacroMol->read($filename);
print " finding bonds..";
find_bonds($mol);
print " done.\n";

my $displaylist;

glutInit;

my $v  = new OpenGL::Simple::Viewer(
        title => 'PerlMol',
        draw_geometry => sub { 
                if (!$displaylist) { build_displaylist(); }
                glCallList($displaylist);
        },
        screenx => 512, screeny => 512,
);

glClearColor(1,1,1,1);

glutMainLoop;

exit 0;



sub recenter {
    my @center;
    my @atoms = $mol->atoms;
    for (@atoms) {
        my @coords = $_->coords->array;
        $center[$_] += $coords[$_] for 0..2;
    }
    $center[$_] /= @atoms for 0..2;
    my $center_v = Math::VectorReal->new(@center);
    for (@atoms) { $_->coords($_->coords - $center_v) }
}

sub make_model {
    my $i = 0;
    recenter();
    my @atoms = $mol->atoms;
    for my $atom ( @atoms ) {
        my $mass   = log( 1 + $atom->mass ) / $mass_scale;
        my $color  = $element_colours{ $atom->symbol } || $colour{cyan};

        my @coords = $atom->coords->array;
        push @ballpoints, [ $color, $mass, @coords ];
    }
    for my $bond ( $mol->bonds ) {
        my ( $from, $to ) = $bond->atoms;
        my @from = $from->coords->array;
        my @to   = $to->coords->array;
        push @ballsticks, [ \@from, \@to ];
    }

}

sub visualize {
    if ( !@ballpoints ) { make_model() }
    for (@ballpoints) {
        my ( $color, $mass, @coords ) = @$_;
        glColor(@$color);
        glPushMatrix;
        glTranslate(@coords);
        glutSolidSphere( $mass , $sphericity, $sphericity);
        glPopMatrix;
    }

    glColor( @{ $colour{'grey'} } );
    glLineWidth(4);
    for (@ballsticks) { glBegin(GL_LINES); glVertex(@$_) for @$_; glEnd; }
}

sub build_displaylist {
        $displaylist = glGenLists(1);
        glNewList($displaylist,GL_COMPILE);
        visualize();
        glEndList();
}
