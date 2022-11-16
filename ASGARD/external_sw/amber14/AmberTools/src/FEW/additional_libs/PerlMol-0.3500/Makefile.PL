use strict;
use warnings;
#use diagnostics;
use inc::Module::Install;

name('PerlMol');
author('Ivan Tubert-Brohman <itub@cpan.org>');
abstract('Perl modules for chemoinformatics and computational chemistry');
license('perl');
version_from('PerlMol.pm');

requires(qw(
    perl                5.006
    Test::Simple        0.45
    Scalar::Util        1.07
    Storable            2.04
    Text::Balanced      1.89
));

bundle(qw(
    Chemistry::File::VRML           0.10
    Chemistry::3DBuilder            0.10
    Chemistry::Mol                  0.35
    Chemistry::Mok                  0.25
    Chemistry::File::SMILES         0.44
    Chemistry::File::SMARTS         0.22
    Chemistry::File::PDB            0.21
    Chemistry::Pattern              0.26
    Chemistry::MidasPattern         0.11
    Chemistry::File::MDLMol         0.20
    Chemistry::InternalCoords       0.18
    Chemistry::MacroMol             0.06
    Chemistry::Ring                 0.18
    Chemistry::Canonicalize         0.10
    Chemistry::Bond::Find           0.21
    Chemistry::File::SLN            0.10
    Chemistry::FormulaPattern       0.10
    Chemistry::File::Mopac          0.15
    Chemistry::File::XYZ            0.11
    Chemistry::Isotope              0.11
    Chemistry::Reaction             0.02
    Math::VectorReal                1.00
    Statistics::Regression          0.15
    Parse::Yapp::Driver             1.05

));

check_nmake();

&WriteAll();


