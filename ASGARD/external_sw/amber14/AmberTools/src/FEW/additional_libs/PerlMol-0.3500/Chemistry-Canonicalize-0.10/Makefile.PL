use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::Canonicalize',
    'VERSION_FROM'	=> 'Canonicalize.pm', 
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => '0.24',
    },
);
