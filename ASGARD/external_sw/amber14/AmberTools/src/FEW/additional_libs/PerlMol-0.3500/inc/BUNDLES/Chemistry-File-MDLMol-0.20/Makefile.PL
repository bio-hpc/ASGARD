use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::MDLMol',
    'VERSION_FROM'	=> 'MDLMol.pm',
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => "0.35" 
    },
);
