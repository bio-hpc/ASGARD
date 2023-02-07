use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::Isotope',
    'VERSION_FROM'	=> 'lib/Chemistry/Isotope.pm',
    'PREREQ_PM'		=> {
        'Test::Simple'      => 0,
    },
);
