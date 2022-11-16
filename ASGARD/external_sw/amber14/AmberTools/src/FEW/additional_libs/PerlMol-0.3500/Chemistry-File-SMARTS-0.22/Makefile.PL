use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SMARTS',
    'VERSION_FROM'	=> 'SMARTS.pm', 
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol          0.24
        Chemistry::Pattern      0.25
        Chemistry::Ring         0.16
        Chemistry::File::SMILES 0.41
        List::Util              0
        Test::More              0
    )}, 
);
