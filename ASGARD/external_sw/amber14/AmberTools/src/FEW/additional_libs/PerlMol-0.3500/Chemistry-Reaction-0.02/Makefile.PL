use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::Reaction',
    'VERSION_FROM'	=> 'lib/Chemistry/Reaction.pm', 
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol          0.26
        Chemistry::Pattern      0.25
        Chemistry::File::SMILES 0.42
    )}, 
);

