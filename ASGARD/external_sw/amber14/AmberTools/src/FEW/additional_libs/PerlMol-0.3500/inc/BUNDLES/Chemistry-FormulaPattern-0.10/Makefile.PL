use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::FormulaPattern',
    'VERSION_FROM'	=> 'lib/Chemistry/FormulaPattern.pm', 
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol      0.24
        Chemistry::Pattern  0.20
        Text::Balanced      0
        Test::More          0
    ),
    },
);
