use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SLN',
    'VERSION_FROM'	=> 'lib/Chemistry/File/SLN.pm', # finds $VERSION
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol              0.26
        Chemistry::Canonicalize     0.10
        Chemistry::Bond::Find       0.21
        Chemistry::Ring             0.15
        Parse::Yapp                 0
        List::Util                  0
        Test::More                  0
    )},
);

sub MY::postamble {
<<'EOT';
PARSER = lib/Chemistry/File/SLN/Parser.pm

#lib/Chemistry/File/SLN/Parse.pm: Parser.yp
$(PARSER) :: Parser.yp
	$(PERL) -S yapp -m 'Chemistry::File::SLN::Parser' -o 'lib/Chemistry/File/SLN/Parser.pm' Parser.yp
EOT
}
