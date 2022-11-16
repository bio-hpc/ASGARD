use 5.006;

use ExtUtils::MakeMaker;
WriteMakefile(
    'NAME'	=> 'Chemistry::File::VRML',
    'VERSION_FROM' => 'lib/Chemistry/File/VRML.pm',
    'PREREQ_PM'		=> {
        'Chemistry::Mol'        => '0.33',
    },
);
