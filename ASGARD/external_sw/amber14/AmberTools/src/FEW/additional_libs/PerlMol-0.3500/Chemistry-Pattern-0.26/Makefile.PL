use 5.006;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Chemistry::Pattern',
    'VERSION_FROM'	=> 'Pattern.pm', # finds $VERSION
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => '0.24',
        'Test::Simple' => 0,
    },
);
