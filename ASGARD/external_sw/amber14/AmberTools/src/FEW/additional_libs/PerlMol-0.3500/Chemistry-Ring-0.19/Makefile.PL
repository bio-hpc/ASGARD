use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::Ring',
    'VERSION_FROM'	=> 'lib/Chemistry/Ring.pm', # finds $VERSION
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol              0.24
        Statistics::Regression      0.15
    )},
);
