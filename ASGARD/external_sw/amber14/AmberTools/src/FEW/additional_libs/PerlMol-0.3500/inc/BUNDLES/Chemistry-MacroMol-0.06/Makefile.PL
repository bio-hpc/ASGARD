use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::MacroMol',
    'VERSION_FROM'	=> 'MacroMol.pm', 
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol      0.24 
        Scalar::Util        0
    )},
);
