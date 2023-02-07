#!/usr/bin/perl

use lib "/home/kfwong/Work/VOTH/ONR/amber10/src/amgauss/lib";
use File::Basename;
use ESPT::Aprmtop 0.02;
use ESPT::Aprmcrd 0.01;
use ESPT::Gfchk 0.02;
use strict;

=head1 NAME

AmGauss - The AMBER <-> Gaussian file generator
   
=head1 SYNOPSIS

B<amgauss> [ B<-b> basis ] [ B<-c> charge ] [ B<--cartesian> ] [ B<-d> ] [ B<-ext> extension ] 
[ B<-h> ] [ B<--ignore-irc> ] [ B<-j> type ]
[ B<-k> keywords ] [ B<-m> multiplicity ] [ B<-mem> memory ] [ B<-o> name ] [ B<-prmcrd> file ]
[ B<-q> ] [ B<-qmmm> layers ] [ B<-t> theory ] [ B<-title> title ] F<filename>

=head1 DESCRIPTION

This program generates Gaussian input files from AMBER prmtop and prmcrd files.  It
also generates AMBER EVB data from Gaussian formatted checkpoint files (fchk) files. 
 
=cut

### Version History ###
# 1.0   Create .inp files from prmtop & prmcrd files
#	Create .EVB files from .fchk files
# 1.1	Options to pass molecular charge. Gaussian calculation
#	title; extra Gaussian route line keywords; Added ONIOM
#	input; quiet option; atomic symbols from integer
#	atomic masses
# 1.2	parsing of L123 IRC files, ignoreIRC option, cartesian option
#	unified EVB writing routine, redundant internal coordinates,
#	
#
# 2007 Jason L. Sonnenberg

### To Do List ###
# Erorr checking 
# PBS commands in .inp files
# handle isotopes using element(iso=X)
# charge multiplicity sanity check
# table of differences in red. int. coords, gradients, Hessian
# ability to run formchk if passed a .chk file

### Main Program ###
our $version = "1.2";

# check for arguments
usage() if ( $#ARGV < 0 );

help() if $ARGV[0] eq "-h";
help() if $ARGV[0] eq "--help";

our ($basis, $cart, $charge, $debug, $extension, $file, $ignoreIRC, $jobtype, @keywords, $logfile);
our ($memory, $multiplicity, $name, $prmcrdfile, $qmmm, $title, $theory);

# parse arguments
for (my $i=0; $i<=$#ARGV; $i++) {
	$basis = $ARGV[$i + 1] if $ARGV[$i] eq "-b";
        $charge = $ARGV[$i + 1] if $ARGV[$i] eq "-c";
	$cart = 1 if $ARGV[$i] eq "--cartesian";
	$debug = 1 if $ARGV[$i] eq "-d";
	$extension = $ARGV[$i + 1] if $ARGV[$i] eq "-ext";
	$ignoreIRC = 1 if $ARGV[$i] eq "--ignore-irc";
	$jobtype = $ARGV[$i + 1] if $ARGV[$i] eq "-j";
	@keywords = split( /,/, uc($ARGV[$i + 1]) ) if $ARGV[$i] eq "-k";
	if ( $ARGV[$i] eq "-o") {
		$name = $ARGV[$i + 1];
		$name =~ s/\.[A-Za-z\d]+$//;
	}
	$multiplicity = $ARGV[$i + 1] if $ARGV[$i] eq "-m";
	$memory = $ARGV[$i + 1] if $ARGV[$i] eq "-mem";
	$prmcrdfile = $ARGV[$i + 1] if $ARGV[$i] eq "-prmcrd";
	$theory = uc($ARGV[$i + 1]) if $ARGV[$i] eq "-t";
	$debug = -1 if $ARGV[$i] eq "-q";
	if ( $ARGV[$i] eq "-qmmm" ) {
		$qmmm = uc($ARGV[$i + 1]) if $ARGV[$i + 1] =~ /.+:.+:*/ ;
		$theory = join "", "ONIOM(", $qmmm, ")";
	}
	$title = $ARGV[$i + 1] if $ARGV[$i] eq "-title";
        $file = $ARGV[$i] if $i == $#ARGV;
} 

# set defaults
$cart ||= 0;
$debug ||= 0;
$qmmm ||= "B3LYP/6-31G(d):AMBER";
$ignoreIRC ||= 0;

=head1 OPTIONS

Command  line  option  specifications are processed from left to right and may 
be specified more than once. If  conflicting  options  are specified, later  
specifications override earlier ones. The filename passed to amgauss must be
either a prmtop or fchk file.  Amgauss automatically determines what type of file 
was passed and locates F<filename>.crd files if present. 

=over 16

=item B<-b> basis

Selects the basis set used by Gaussian. The default is 6-31G(d). Please 
see the Gaussian manual for available basis sets.

=item B<-c> charge

Sets the molecular charge.  Charge must be an integer.  Default molecular charge
is computed from the molecular mechanic charges.

=item B<--cartesian>

Include Cartesian gradient and Hessian data in the EVB file.

=item B<-d>

Turn on debug printing

=item B<-ext> extension

Sets the extension of the resulting file. The defaults are inp for Gaussian input files and EVB 
for AMBER data files.

=item B<-h>

=item B<--help>

Print full AmGauss documentation via perldoc. Can not be used with other options.

=item B<--ignore-irc>

Ignore any IRC data contained in the fchk file.  Mainly used for debugging.

=item B<-j> type

Sets the Gaussian job type. An optimization followed by calculation of the frequencies (Opt 
Freq) is default. Please see the Gaussian manual for available job types.

=item B<-k> keywords

Extra keywords for the Gaussian route line.  Multiple keywords should be separated by commas 
without separating spaces. ONIOM theory levels should be passed via B<-qmmm>. Please see the 
Gaussian manual for available keywords.

=item B<-m> multiplicity

Molecular spin multiplicity, 2S+1. Defaults to the number of unpaired electrons + 1.

=item B<-mem> memory

Amount of dynamic memory, passed via the %mem command, to be used by 
Gaussian. The amount of memory is followed by B<KB>,
B<KW>, B<MB>, B<MW>, B<GB>, B<GW> with no intervening spaces. Defaults to 500MB.

=item B<-o> name

Output file name. The default is to use F<filename> with all extensions and paths removed.

=item B<-prmcrd> file

AMBER prmcrd filename.  If the prmcrd filename is not F<filename>.crd then this option is 
mandatory for creating EVB files.

=item B<-q>

Run in quite mode and do not print progress messages.

=item B<-qmmm> layers

Request a QM/MM calculation in Gaussian. Two or three colon seperated theory levels can be 
passed for the high, medium and low layers in ONIOM. To specify which atoms belong to each
layer, please open the resulting Gaussian input file in GaussView. Defaults to B3LYP/6-31G(d):AMBER.

=item B<-t> theory

Sets the theory level used by Gaussian.  The default is Density Functional 
Theory using the B3LYP functional.  For available theory levels, please 
consult the Gaussian manual.	

=item B<-title> title

Sets the title used in the Gaussian calculation. Titles containing spaces must be passed in double quotes.
Defaults to the title contained in the AMBER prmtop file.

=back

=cut

# check for prmtop or fchk file
open(FILEIN,$file) || die "Could not read $file\n$!\n";

# get path and base filename
(my $base, my $dir, my $ext) = fileparse($file, qr/\.[topfchk]*/);
$name ||= $base;

our ($input, $prmcrd);

# determine file type and set extension
while (<FILEIN>){
        # skip blank lines
        next if /^$/;
	
	# AMBER prmtop
	if ( /^%FLAG\sATOM_NAME/ ) {
		$input = ESPT::Aprmtop->new();
		$extension ||= "inp";
		last;
	} 

        # Gaussian fchk
        if ( /^Number\s+of\s+atoms\s+I\s+\d+/ ) {
                $input = ESPT::Gfchk->new();   
      		$extension ||= "EVB";
         	last;
        }
}
close(FILEIN);

print "Processing a ", $input->get("PROGRAM"), " ", $input->get("TYPE"), " file.\n" if $debug >= 0;

# read file contents
$input->debug($debug);
$input->analyze($file, "Alpha");

# Build output file name
our $out = join "", ">", $name, ".", $extension;

# Grab other files, extract data and write output.
if ( $input->get("PROGRAM") eq "AMBER") {

	# open the prmcrd file
	$prmcrdfile ||= join "", $dir, $name, ".crd";
	open (PRMCRDFILE, $prmcrdfile) || die "Could not open $prmcrdfile\n$!\n";

	$prmcrd = ESPT::Aprmcrd->new();
	$prmcrd->debug($debug);
	$prmcrd->analyze($prmcrdfile, "Alpha");
	print "Read file ", $prmcrdfile, "\n" if $debug >= 0;

	# Set Gaussian job defaults
	# Precedence: command line > built-in defaults
	$basis ||= "6-31G(d)";
	$memory ||= "500mb";
	$jobtype ||= "OPT FREQ";
	$theory ||= "B3LYP";
	$multiplicity ||= $input->get("MULTIPLICITY");

	## write Gaussian input file ##
	writeGaussian($out);

}else {
	## write AMBER EVB file ##

	# Data is printed in 5 columns of scientific notation 
	# with an eight decimal mantissa. 

	# write out a seperate .EVB file for each IRC point
	# Note that Hessian data is not computed for standard
	# IRC runs.
	if ( $input->get("IRCPOINTS") > 0 && $ignoreIRC == 0) {

	  print $input->get("IRCPOINTS"), " IRC data points found.\n";
		
	  for (my $i=0; $i<$input->get("IRCPOINTS"); $i++) {
			
		# Adjust output file name. If the step sizes are
		# identical then files will be overwritten. 
		# use a mantissa of five to match table in Gaussian logfile
		$out = join "", ">", $name, "_", sprintf("%.5f", $input->get("IRCSTEP", $i)), ".", $extension;
		writeEVB($out, 1, $i);			
	  }
	} else {
          # write out .EVB file for non-IRC runs
	  writeEVB($out);
	}
}

print "Successfully wrote ", $extension, " file(s).\n" if $debug >= 0;

## Subroutines ##

# display help on usage
sub help {
	system("perldoc amgauss");
	exit;
}

sub usage {
        print "\nAmGauss $version AMBER <-> Gaussian file generator\n";
        print "\nUsage: amgauss [options] filename \n";
	print "\t-b basis\tBasis set, default is 6-31G(d)\n";
	print "\t--cartesian\tInclude Cartesian gradient & Hessian in EVB file\n";
	print "\t-c charge\tMolecular charge\n";
        print "\t-d \t\tDebug print \n";
        print "\t-ext extension\tOutput file extension\n";
        print "\t\t\tdefault is inp for Gaussian, EVB for AMBER\n";
	print "\t-h\t\tprint full documentation\n";
	print "\t--ignore-irc\tIgnore IRC data in fchk file\n";
	print "\t-j type\t\tGaussian job type, default is Opt Freq\n";
	print "\t-m multiplicity\tMolecular spin state\n";
	print "\t\t\tdefaults to number of unpaired electrons + 1\n";
	print "\t-mem memory\tGaussian memory allocation, defaults to 500MB\n";
	print "\t-o name\t\tOutput file name, defaults to filename\n";
	print "\t-prmcrd file\tAMBER prmcrd filename\n";
	print "\t-q\t\tQuiet mode\n";
	print "\t-qmmm layers\tQM/MM theory levels, defaults to B3LYP/6-31G(d):AMBER\n";
	print "\t-t theory\tTheory level, default is B3LYP\n";
	print "\t-title title\tGaussian calculation title, defaults to prmtop title\n";
	print "\tfilename\tA valid prmtop or fchk file.\n";
	print "\n";
        exit;
}

# write out an EVB data file
# eventually this should be made completely general and have everything
# passed as parameters; note the use of prototypes
sub writeEVB($;$$) {
	# process parameters
	my ($filename, $ircflag, $ircpoint) = @_;
	
	# establish defaults
	my $handle = \*FILEOUT;
	$handle = \*STDOUT if $debug >= 1;
	$ircflag ||=0;

       	# Open output file
	print "EVB output file name is $filename\n" if $debug >= 1;
        unless ( $debug == 1) {
          	open($handle, $filename) || die "Could not create $filename\n$!\n";
                print "Opened file ", $filename, "\n" if $debug >= 0;
        }

	# print data dimension(s)
	print $handle "[external evb redundant internal data dimension]\n";
	print $handle " ", $input->get("NREDINT"),"\n\n";

	if ($cart == 1) {
	        print $handle "[external evb cartesian data dimension]\n";
        	print $handle " ", 3 * $input->get("NATOMS"), "\n\n";
	}

	# print coordinates
	print $handle "[redundant internal coordinates]\n";
	for (my $i=0; $i<$input->get("NREDINT"); $i++) {
	  for (my $j=0; $j<4; $j++) {
		print $handle sprintf " % 7u ", $input->get("REDINTCOORD",$i, $j);
		print $handle "\n" if (4*$i + $j + 1)/4 - int((4*$i + $j + 1 )/4) == 0;
	  }
	}
	print $handle "\n";
	
        print $handle "[cartesian coordinates]\n";
        for (my $i=0; $i<$input->get("NATOMS"); $i++) {
                # x, y, z coords.
                for (my $j=0; $j<3; $j++) {
                  if ( $ircflag == 0 ) {      
			print $handle sprintf " % 9.8E", $input->get("CARTCOORD", $i, $j);
		  } else {
			print $handle sprintf " % 9.8E", $input->get("IRCCOORD", $ircpoint, $i, $j);
		  }
                        print $handle "\n" if (3*$i + $j + 1)/5 - int((3*$i + $j + 1)/5) == 0;
                }
        }
	print $handle "\n" if ( 3*$input->get("NATOMS") % 5 );
	print $handle "\n";

	# electronic energy
        print $handle "[electronic energy]\n";
	if ( $ircflag == 0 ) {
	        print $handle " ", $input->get("EELEC"), "\n\n";
	} else {
		print $handle " ", $input->get("IRCENERGY", $ircpoint), "\n\n";
	}
	
	# gradient(s)
	print $handle "[redundant internal gradient]\n";
	for (my $i=0; $i<$input->get("NREDINT"); $i++) {
		print $handle sprintf " % 9.8E", $input->get("REDINTGRADIENT", $i);
             	print $handle "\n" if ($i+1)/5 - int(($i+1)/5) == 0;
        }
	print $handle "\n" if ( $input->get("NREDINT") % 5 );
	print $handle "\n";

	if ( $cart == 1 ) {
	        print $handle "[cartesian gradient]\n";
        	for (my $i=0; $i<$input->get("NATOMS")*3; $i++) {
 		  if ( $ircflag == 0) {
                	print $handle sprintf " % 9.8E", $input->get("GRADIENT", $i);
	  	  } else {
			print $handle sprintf " % 9.8E", $input->get("IRCGRADIENT", $ircpoint, $i);
		  }
                  print $handle "\n" if ($i+1)/5 - int(($i+1)/5) == 0;
        	}
	        print $handle "\n" if ( 3*$input->get("NATOMS") % 5 );
        	print $handle "\n";
	}

	# Hessian data
	print $handle "[redundant internal hessian]\n";
	for (my $i=0; $i<$input->get("NREDINT"); $i++) {
	  for (my $j=0; $j<$input->get("NREDINT"); $j++) {
              print $handle sprintf " % 9.8E", $input->get("REDINTHESSIAN", $i, $j);
              print $handle "\n" if ($input->get("NREDINT")*$i + $j + 1)/5 - int(($input->get("NREDINT")*$i + $j + 1)/5) == 0;
	  }
	}
        print $handle "\n" if ( $input->get("NREDINT")*($input->get("NREDINT")+1)*0.5 % 5 );
        print $handle "\n"; 

	unless ( $cart == 0 || $ircflag == 1 ) {
          print $handle "[cartesian hessian]\n";
          for (my $i=0; $i<$input->get("NATOMS")*3; $i++) {
            for (my $j=0; $j<$input->get("NATOMS")*3; $j++) {
              print $handle sprintf " % 9.8E", $input->get("HESSIAN", $i, $j);
              print $handle "\n" if ($input->get("NATOMS")*3*$i + $j + 1)/5 - int(($input->get("NATOMS")*3*$i + $j + 1)/5) == 0;
            }
      	  }
          print $handle "\n" if ( (3*$input->get("NATOMS"))*(3*$input->get("NATOMS")+1)*0.5 % 5 );
          print $handle "\n"; 
	}

	# atomic masses
        print $handle "[mass]\n";
        for (my $i=0; $i<$input->get("NATOMS"); $i++) {
     	       print $handle  sprintf " % 9.8E", $input->get("MASS", $i);
               print $handle "\n" if ($i+1)/5 - int(($i+1)/5) == 0;
        }
	print $handle "\n" if ( $input->get("NATOMS") % 5 );
        print $handle "\n";

	# Redundant internal TS vector
#	print $handle "[redundant internal TS vector]\n";	
        close ($handle) if $debug <= 0;
}

# write out a Gaussian input file
# eventually this should be made completely general and have everything
# passed as parameters; note the use of prototypes
sub writeGaussian($) {

	# process parameters
	my ($filename) = shift;

        my $handle = \*FILEOUT;
        $handle = \*STDOUT if $debug >= 1;

	# Open output file
	print "Gaussian output file name is $filename\n" if $debug >= 1;
	unless ( $debug == 1) {
		open($handle, $filename) || die "Could not create $filename\n$!\n";
		print "Opened file ", $filename, "\n" if $debug >= 0;
	}

	# % commands
	print $handle "%chk=$name\n%mem=$memory\n\n";

	# Route
	print $handle "#P ", uc($theory), " ";
	print $handle uc($basis), " " unless $theory =~ m/ONIOM/;
 	print $handle uc($jobtype);
	for (my $i=0; $i<$#keywords + 1; $i++) {
		print $handle " ", $keywords[$i];
	}
	print $handle "\n\n";

	# title
	$title ||= $input->get("TITLE");
	print $handle $title, "\n\n";

	# Charge & Multiplicity
	$charge ||= $input->get("CHARGE");
	print $handle $charge, " ", $multiplicity, "\n";

	# Atoms & cartesian coordinates
	for (my $i=0; $i<$input->get("NATOMS"); $i++) {
		#print atomic symbols for clarity, although Gaussian accepts atomic numbers
		print $handle $input->get("ATOMS", $i);
		for (my $j=0; $j<3; $j++) {
			print $handle sprintf " % 10.7f", $prmcrd->get("CARTCOORD", $i, $j);
		}
		print $handle "\n";
	}
	print $handle "\n";

	print "Charge = ", $input->get("CHARGE"), " Multiplicity = ", $input->get("MULTIPLICITY"), "\n" if $debug >= 0;
        close ($handle) if $debug <= 0;
}

1;

=head1 EXAMPLES

=over

=item amgauss

Called with no parameters at all, Amgauss will display usage information. 
If B<-h> or B<--help> is passed then the full Amgauss documentation is displayed via perldoc.

=item amgauss foo.top

Amgauss reads foo.top and foo.crd files from AMBER then creates a Gaussian input file named foo.inp.
If foo.crd is not present in the same place as foo.crd then Amgauss will fail.

=item amgauss -prmcrd foo.xyz -ext com -m 2 foo.top

Amgauss reads foo.top and uses foo.xyz as the AMBER prmcrd file to create the Gaussian input file foo.com.
The multiplicity in the .com file is set to 2.

=item amgauss -title "TS optimization for EVB" -c 2 -k nosymm,pop=nbo foo.top

Amgauss reads foo.top and foo.crd files from AMBER and then creates a Gaussian input file named foo.inp. 
In foo.inp the molecular charge is set to 2, NOSYMM and POP=NBO are added to the Gaussian route, line 
and the calculation title is TS optimization for EVB.

=item amgauss -o test foo.fchk

Amgauss reads foo.fchk file from Gaussian and creates test.EVB necessary for AMBER EVB jobs.

=back

=head1 NOTES

Amgauss uses the atomic masses found in the prmtop file to determine atomic symbols.  For elements 
possessing the same integer atomic mass, such as Ar and Ca, Amgauss relies on the atom name section of
the prmtop file.  If those atom names are not atomic symbols followed by an optional integer and separated
by a space(s), then the resulting Gaussian input file will require editing.

=head1 VERSION

1.2

=head1 AUTHOR

Jason L. Sonnenberg, E<lt>sonnenberg.11@osu.eduE<gt>

=cut


