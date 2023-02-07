#!/home/ivan/bin/perl
# usage: ./pep.pl <SEQUENCE>
# prints an mdl molfile to stdout

# curent version only supports the residues A,G,L,V
# the residues are assumed to be in the files listed below in %files,
# which must have the atoms N, CA, and C as the first three atoms in the 
# file

use strict;
use warnings;
#use diagnostics;

use Chemistry::File::MDLMol;
use Chemistry::File::SMARTS;
use Chemistry::InternalCoords::Builder qw(build_zmat);

{
    # the files with the template aminoacids (3D coords)
    my %files = (
        A => 'ala.mol',
        V => 'val.mol',
        G => 'gly.mol',
        L => 'leu.mol',
    );

    # This defines the essential part of an aminoacid
    my $patt = Chemistry::Pattern->parse('C((=O)CN([C,H])[H])O[H]', 
        format => 'smarts');

    my %aa;

    my ($c, $o, $ca, $n, $hn);

    # return a molecule object from an aminoacid given its one-letter code
    sub amino_acid {
        my ($code) = @_;
        return $aa{$code} if $aa{$code};
        my $file = $files{$code} or die "unknown code '$code', stopping ";
        $aa{$code} = Chemistry::Mol->read($file);
    }

    # start a new chain from a residue (a Chemistry::Mol object).
    # returns the new chain (a Chemistry::Mol object)
    sub start_chain {
        my ($res) = @_;

        # clone because we don't want to destroy the template
        my $chain = safe_clone($res);

        # find out which atom is which
        $patt->match($chain) or die "didn't match, stopping ";
        my ($oh, $ho, $hnb);
        ($c, $o, $ca, $n, $hn, $hnb, $oh, $ho) = $patt->atom_map;

        # chop off the carboxyl OH
        $oh->delete; $ho->delete;

        # generate internal coordinates. 
        build_zmat($chain, bfs => 0);
        $chain;
    }

    # add a new residue (a Chemistry::Mol object) to a chain
    sub add_residue {
        my ($chain, $res) = @_;

        # clone because we don't want to destroy the template
        $res = safe_clone($res);

        # find out which atom is which
        $patt->match($res) or die "didn't match, stopping ";
        my @map = $patt->atom_map;

        # delete the hydroxyl and an HN atom
        $_->delete for splice @map, 5;

        # generate internal coordinates. 
        # The bfs => 0 is needed to make sure
        # that the first three atoms are N, CA, and C, so that we can use 
        # them reliably for positioning the new residue with respect to
        # the previous one
        build_zmat($res, bfs => 0);

        # add the new residue
        $chain->combine($res);

        # since combine clones $res, we want to find the atoms of interest
        # in the new part of the chain, by "translating" the atom_map
        my ($c2, $o2, $ca2, $n2, $hn2) = 
            map { $chain->by_id($_->id) } @map;

        # create the peptide bond
        $chain->new_bond(atoms => [$n2, $c], order => 1);

        # add some roughly reasonable internal coordinates to position the
        # new residue
        $n2 ->internal_coords($c, 1.5, $ca, 120, $o, 180);
        $ca2->internal_coords($n2, 1.5, $c, 120, $ca, 180);
        $c2 ->internal_coords($ca2, 1.5, $n2, 120, $c, 180);

        # make sure that the N is planar
        $hn2->internal_coords($n2, 1.1, $c, 120, $o, 180)
            if $hn2->symbol eq 'H'; # make sure it's not Proline

        # save the new terminal residue's atoms of interest
        ($c, $o, $ca, $n, $hn) = ($c2, $o2, $ca2, $n2, $hn2);
    }

    sub end_chain {
        my ($chain) = @_;

        # add the terminal hydroxyl group
        my $ox = $chain->new_atom(
            symbol => 'O', 
            internal_coords => [ $c, 1.5, $ca, 120, $o, 180 ]
        );
        my $h = $chain->new_atom(
            symbol => 'H', 
            internal_coords => [ $ox, 1.1, $c, 105, $o, 0 ]
        );
        $chain->new_bond(atoms => [$ox, $c], order => 1);
        $chain->new_bond(atoms => [$ox, $h], order => 1);

        # finally, calculate the cartesian coordinates for all the atoms in
        # the chain
        $_->internal_coords->add_cartesians for $chain->atoms;
    }

}

# this is like Chemistry::Mol::clone, but assigns a fresh ID to the cloned
# molecule, atoms, and bonds. This should be added to Chemistry::Mol
sub safe_clone {
    my ($mol) = @_;
    my $clone = $mol->clone;
    for ($clone, $clone->atoms, $clone->bonds) {
        # this uses the undocumented nextID method in Chemistry::Mol,
        # Chemistry::Atom, and Chemistry::Bond
        $_->id($_->nextID);
    }
    $clone;
}

# takes a sequence string such as "LAVA" and returns a Chemistry::Mol object
sub build_sequence {
    my ($s) = @_;
    $s or die "no sequence, stopping ";
    my @seq = split //, $s;

    my $code = shift @seq;
    my $res = amino_acid($code);

    my $chain = start_chain($res);
    $chain->name($s);

    for my $code (@seq) {
        $res = amino_acid($code);
        add_residue($chain, $res);
    }
    end_chain($chain);

    return $chain;
}


##############################

# generate the structure and print it out as an MDL molfile
my $chain = build_sequence(shift);
print $chain->print(format => 'mdl');

