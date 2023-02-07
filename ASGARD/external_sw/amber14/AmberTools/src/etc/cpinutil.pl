#!/usr/bin/perl -w
BEGIN {                         # Add executable dir to INC to find CPin.pm
  $0 =~ /^(.*?)\/([^\/]+)$/;	# $1 Path $2 filename
  push @INC, $1;
}
use strict;
use IO::File;
use Getopt::Long;
use CPin;

main();

sub main {
  my @origArgs = @ARGV;
  my ($pH, $stepPeriod, $dummy, $stateStr,$stateFile, $sysName, 
      $minpKa, $maxpKa, $resNumStr, $resNameStr, $notResNumStr, $notResNameStr);
  GetOptions("states=s"=>\$dummy, # Clear flags to find filename
	     "outfile|statefile=s" => \$dummy, "system=s" => \$dummy,
	     "minpKa=f" => \$dummy, "maxpKa=f" => \$dummy, 
	     "resnum=s" =>\$dummy, "resname=s" =>\$dummy,
             "notresnum=s" => \$dummy, "notresname=s" => \$dummy,
             "help|usage|?" => \&usageMessage
            );
  
  my $file = defined($ARGV[0]) ? $ARGV[0] :  '-';
  @ARGV = @origArgs;
  my $fh = new IO::File($file);
  $_ = $fh->getc();
  $fh->ungetc(ord $_);
  my $cpin;
  if (/^\d/) {                  # Old style charges file
    $cpin = new CPinFile($fh);
  } elsif ($_ eq '&' or $_ eq ' ') {         # New (namelist) style charges file
    $cpin = new CPinNamelist($fh);
  } else {
    $cpin = new CPinFromPDBFile($fh);
    $cpin->SysName("Unknown");
  }
  GetOptions("states=s" => \ $stateStr, # Have input stream, now really process flags
	     "outfile|statefile=s" => \ $stateFile,
	     "system=s" => \$sysName,
	     "minpKa=f" => \$minpKa,
	     "maxpKa=f" => \$maxpKa,
	     "resnum=s" => \$resNumStr,
             "resname=s" =>\$resNameStr,
             "notresnum=s" => \$notResNumStr,
             "notresname=s" => \$notResNameStr
	    );
  if (defined($sysName)) {
    $cpin->SysName($sysName);
    $cpin->setpKas() if $cpin->isa('CPinFromPDBFile');
  }
  if (defined($stateStr)) {
    $cpin->setStates(split(/,/,$stateStr));
  }
  if (defined($stateFile)) {
    $cpin->setStatesFromOutputFile($stateFile);
  }

  # Filter residues
  my $i = 0;
  my (%SelResNums,%SelResNames,%UnselResNums, %UnselResNames);
  %SelResNums = map {$_ => $_} split(/,/, $resNumStr) if defined $resNumStr;
  %SelResNames = map {uc($_) => $_} split(/,/, $resNameStr) if defined $resNameStr;
  %UnselResNums = map {$_ => $_} split(/,/, $notResNumStr) if defined $notResNumStr;
  %UnselResNames = map {uc($_) => $_} split(/,/, $notResNameStr) if defined $notResNameStr;
  
  while ($i < $cpin->ResCount()) {
    my $res = $cpin->getRes($i);
    my $deleteRes = 0;
    if (defined $minpKa and defined $res->pKa() and $res->pKa() < $minpKa) {
      $deleteRes = 1;
    }
    if (defined $maxpKa and defined $res->pKa() and $res->pKa() > $maxpKa) {
      $deleteRes = 1;
    }
    if (defined $resNumStr and not exists $SelResNums{$res->ResNum()}) {
      $deleteRes = 1;
    }
    if (defined $resNameStr and not exists $SelResNames{uc($res->ResName())}) {
      $deleteRes = 1;
    }
    if (defined $notResNumStr and exists $UnselResNums{$res->ResNum()}) {
      $deleteRes = 1;
    }
    if (defined $notResNameStr and exists $UnselResNames{uc($res->ResName())}) {
      $deleteRes = 1;
    }
    if ($deleteRes) {
      $cpin->DeleteResidueNum($i);
    } else {
      $i++;
    }
  }
  print $cpin->outputNml();
}

sub usageMessage() {
  print <<USAGE;
    
CPinUtil -- utility for manipulating constant pH MD charges files
  
cpinutil.pl <cpinfile or PDB file> [-states=<val,val...>] [-system=<string>]
               [-minpKa=<val>] [-maxpKa=<val>] [-resnum=<val,val...>]
               [-resname=<val,val...>] [-notresnum=<val,val...>]
               [-notresname=<val,val...>]

   cpinfile or PDB file  existing cpin file, cprestrt file or PDB file for system
                         (can also read these from STDIN)
   states                list of initial states for residues
   system                System name (for looking up experimental pKa values)

 The following options delete residues from the cpin file so they don't titrate
   minpKa                Eliminate residues with experimental pKa less than val
   maxpKa                Eliminate residues with experimental pKa higher than val
   resnum                Eliminate residues numbers *not* in the list
   notresnum             Eliminate residues with number that *are* in the list
   resname               Eliminate residue types *not* in the list
   notresname            Eliminate residue types that *are* in the list

 Example:

 ambpdb <prmcrd | cpinutil.pl -system HEWL -maxpKa 7 -notresname HIP >cpin

 This example would create a cpin file from prmtop and prmcrd, loading 
 experimental pKa values for lysozyme. Only non-histidine residues with 
 experimental pkA less than 7 would titrate.

USAGE
exit;
}

