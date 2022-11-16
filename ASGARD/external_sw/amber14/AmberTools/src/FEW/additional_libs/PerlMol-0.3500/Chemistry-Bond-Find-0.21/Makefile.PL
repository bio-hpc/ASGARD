use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::Bond::Find',
    'VERSION_FROM'	=> 'Find.pm', # finds $VERSION
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => '0.24',
    },
);
