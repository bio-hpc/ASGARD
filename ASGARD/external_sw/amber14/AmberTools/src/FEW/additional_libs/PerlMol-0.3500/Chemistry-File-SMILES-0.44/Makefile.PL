use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SMILES',
    'VERSION_FROM'	=> 'SMILES.pm', 
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol              0.30
        Chemistry::Canonicalize     0.10
        Chemistry::Bond::Find       0.21
        Chemistry::Ring             0.15
        List::Util                  0
    )},
);
