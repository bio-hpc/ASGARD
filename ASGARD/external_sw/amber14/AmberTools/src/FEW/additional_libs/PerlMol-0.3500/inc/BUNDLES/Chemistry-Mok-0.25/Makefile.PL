use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Chemistry::Mok',
    VERSION_FROM        => 'lib/Chemistry/Mok.pm', 
    PREREQ_PM           => {
        Chemistry::Mol          => 0.30,
        Chemistry::Pattern      => 0.21,
        Chemistry::File::SMARTS => 0.11,
        Chemistry::Bond::Find   => 0.21,
        Chemistry::Isotope      => 0.10,
        Chemistry::3DBuilder    => 0.10,
        Test::Simple            => 0,
        Scalar::Util            => 1.01,
        Text::Balanced          => 1.87,
        Math::VectorReal        => 1.0,
    },
    EXE_FILES           => ['mok'],
);
