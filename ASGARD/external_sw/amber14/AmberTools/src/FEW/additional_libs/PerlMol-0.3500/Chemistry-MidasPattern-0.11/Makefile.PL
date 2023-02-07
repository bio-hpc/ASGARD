use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::MidasPattern',
    'VERSION_FROM'	=> 'lib/Chemistry/MidasPattern.pm', 
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol         0.24
        Chemistry::MacroMol    0.05
        Chemistry::Pattern     0.20
        Chemistry::File::PDB   0.21
        Test::Simple           0
    )},
);
