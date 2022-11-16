use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'              => 'Chemistry::InternalCoords',
    'VERSION_FROM'      => 'InternalCoords.pm',
    'PREREQ_PM'         => {
        "Chemistry::Mol" => "0.25",
        "Scalar::Util" => 0,
        "Chemistry::Canonicalize" => 0.10,
        "Chemistry::Bond::Find" => 0.20,
    },
    'PMLIBDIRS'         => ['InternalCoords'],
);
