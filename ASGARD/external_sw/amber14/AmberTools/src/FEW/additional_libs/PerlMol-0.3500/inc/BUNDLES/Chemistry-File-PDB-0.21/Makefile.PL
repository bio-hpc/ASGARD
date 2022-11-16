use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::PDB',
    'VERSION_FROM'	=> 'PDB.pm', 
    'PREREQ_PM'		=> {
        'Chemistry::MacroMol'   => 0.05,
        'Chemistry::Mol'        => '0.30',
    },
);
