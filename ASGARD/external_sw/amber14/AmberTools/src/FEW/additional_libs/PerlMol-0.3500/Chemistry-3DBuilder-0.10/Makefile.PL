use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Chemistry::3DBuilder',
    VERSION_FROM        => 'lib/Chemistry/3DBuilder.pm', 
    PREREQ_PM           => {qw(
        Math::VectorReal                1.00
        Chemistry::Mol                  0.33
        Chemistry::InternalCoords       0.18
        Chemistry::Ring                 0.18
        Chemistry::Pattern              0.25
        Chemistry::File::SMILES         0.43
        Chemistry::File::SMARTS         0.21
    )},
);
