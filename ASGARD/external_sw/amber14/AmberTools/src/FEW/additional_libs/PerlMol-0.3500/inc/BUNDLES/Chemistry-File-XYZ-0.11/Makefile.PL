use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::XYZ',
    'VERSION_FROM'	=> 'XYZ.pm', 
    'PREREQ_PM'		=> {
        'Chemistry::Mol'        => '0.24',
    },
);
