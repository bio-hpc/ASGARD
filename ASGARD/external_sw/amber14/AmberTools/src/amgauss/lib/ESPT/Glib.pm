package ESPT::Glib;
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(rparser);

=head1 NAME

ESPT::Glib - Gaussian library module

=head1 SYNOPSIS

        use ESPT::Glib;

	rparser($object);

=head1 DESCRIPTION

This module contains subroutines for analzing Gaussian files.

=cut

our $VERSION = '0.01';

### Version History ###
# 0.01	rparser(jobtype)

### To Do ###

=head1 SUBROUTINES

=over 10

=item B<rparser($object)>

Gaussian route line parser.

=back

=cut

# parse the route line
sub rparser {

# grab the object to store data in
my $gauss = shift;

# grab keywords
my @keywords = split /(?<!,)\s+/, $gauss->{ROUTE};

# keyword regexpressions
my @bases = ("gen", "[ceopst346]+-[1-3]+[\+]*g", "d95v*", "shc", "lanl","tz", "(?:aug-)*cc-pv[dqt56]z"); 
push @bases, ("sv", "sdd", "midix", "epr",  "ugbs", "mtsmall", "dg[dtz]+vp", "6-31g");

my @exchange = ("h*f*[sb](?:handh)*", "xa(?:lpha)*", "pw91", "mpw", "g96", "m*pbe", "o", "vsxc");
push @exchange, ( "hcth", "tpss", "lsda");

my @jobtypes = ("sp", "opt","ts", "freq", "irc(?:max)*", "scan", "polar", "admp", "bomd", "force");
push @jobtypes, ("stable", "volume", "density=check", "guess=only", "rearchive", "mixed", "saddle");

my @theories = ("amber", "dreiding", "uff","[cimz]+ndo", "am1", "pm3m*", "hf","mp[2-5]");
push @theories, ("ci", "cc[ds]{1,2}", "qci", "g[1-3]", "cbs", "w1", "cas", "gvb", "sac-ci");
                
# parsing
KEY: foreach (@keywords) {
        
        # save to KEYWORDS
        push @{$gauss->{KEYWORDS}}, $_;
                
        # print options
        next KEY if (/^#[npt]*\Z/ );

	# Job Type
	# SP runs using theory/basis notation
	$gauss->{JOBTYPE} = "SP" if ( /.+\/.+/ && $gauss->get("JOBTYPE") eq "undef" );

	# OPT-SP runs using theory/basis//theory/basis notation
	$gauss->{JOBTYPE} = "OPT SP" if ( /.+\/.+\/\/.+\/.+/ && $gauss->get("JOBTYPE") eq "undef" );

	J: foreach my $jt (@jobtypes) {
		if ( /^([fp]*$jt)/ ) {
		  my $tmp = $1;
		  # account for Opt Freq runs
		  if ( $gauss->{COMPLETE} == 0 && ($gauss->{JOBTYPE} =~ /OPT/ ) ) {
			$gauss->{JOBTYPE} = join " ", $gauss->{JOBTYPE}, uc($tmp);
			next KEY;
		  } else {	
			$gauss->{JOBTYPE} = uc($tmp);
			next KEY;
		  }
		}
	}
         
        # theory
        T: foreach my $theory (@theories) {
                if ( /^((?:[ur]*)$theory[a-b0-9\(\)]*)\/*/ ) {
                        $gauss->{THEORY} = uc($1);
                        next KEY unless ( /\// );
                }
        }
         
        # keywords with options
        next KEY if ( /[=\(](?![1-3dpf,]{1,7}\))/ );
                
        # functional
        F: foreach my $functional (@exchange) {
                if ( /^((?:[ur])*$functional(?:[13]*)t*[belmnpsvwy125-9]*)\/*/ ) {
                        $gauss->{THEORY} = "DFT";
                        $gauss->{FUNCTIONAL} = uc($1);
                        $gauss->{FUNCTIONAL} =~ s/AND/and/;
                        next KEY unless ( /\// );
                }
        }
        # basis set
        B: foreach my $basis (@bases) {
                if ( /\/*($basis(?:.*))/ ) {
                        $gauss->{BASIS} = uc($1);
                        # enumerate the * & ** notation
                        $gauss->{BASIS} =~ s/\*\*/(d,p)/;
                        $gauss->{BASIS} =~ s/\*/(d)/;
                        next KEY unless ( /\// );
 
                }
        }
}
print "Gaussian job type = ", $gauss->{JOBTYPE}, "\n" if $gauss->{DEBUG} >= 1;
}

1;
__END__

=head1 VERSION

0.01

=head1 AUTHOR

Jason L. Sonnenberg, E<lt>sonnenberg.11@osu.edu<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2006 by Jason L. Sonnenberg

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. I would like to hear of any
suggestions for improvement.

=cut
