#/usr/bin/perl

# Script for automatic extraction of energy data from *statistics.out files of MM-PBSA / MM-GBSA
# calculations for a set of ligands.
# 
# Usage: perl extract_WAMMenergies.pl <structure file> <path> pb<no.>_gb<no.> <Start>_<Stop>_<Offset>
#
# structure file : Text-file containing list of names of ligand structures that shall be considered;
#                  one ligand name per line.
# path : Directory containing the results of the MM-PBSA and/or MM-GBSA calculations,
#        e.g.: /home/work/FEW/calc_r_1t
# pb<no.>_gb<no.> : FEW internal IDs for types of MM-PBSA / MM-GBSA calculation used.
# <Start>_<Stop>_<Offset> : Numbers of snapshots taken into consideration in MM-PBSA / MM-GBSA
#                           calculations. Start, Stop, and Offset need to be identical to the
#                           values specified for the respective flags in the WAMM command file.
###################################################################################################################

my $ref_file = $ARGV[0];
my $calc_folder = $ARGV[1];
my $pb_gb_type = $ARGV[2];
my $snaps_inter = $ARGV[3];

# Checking input
if(@ARGV < 4){
	print "Usage: perl extract_WAMMenergies.pl <structure file> <path> pb<no.>_gb<no.> <Start>_<Stop>_<Offset>\n";
	exit;
}
if(!($pb_gb_type =~ m/pb(\d)\_gb(\d)/)){
	print "The input flag for specification of the type of MM-PBSA / MM-GBSA calculation\n";
	print "that shall be considered does not have the expected pattern. It should have\n";
	print "the format 'pb<no.>_gb<no.>'\n";
	exit;
}
if(!($snaps_inter =~ m/(\d+)\_(\d+)\_(\d+)/)){
	print "The input flag for specification of the snapshots that were considered in the\n";
	print "MM-PBSA / MM-GBSA analysis from which energy values shall be extracted does not\n";
	print "have the expected pattern. The snapshot information need to be provided in\n";
	print "the format '<Start>_<Stop>_<Offset>\n";
	exit;
}

# Read information about structures to regard
my @structs;

open(REF, $ref_file) || die "Cannot open file $ref_file.\n";
while(my $l = <REF>){
	chomp($l);
	push(@structs, $l);
}
close REF;

# Read information from *statistics.out files
my $start = 0;
my $ele = 0;
my $vdw = 0;
my $pbsur = 0;
my $pbcal = 0;
my $pbtot = 0;
my $gb = 0;
my $gbsur = 0;
my $gbtot = 0;

# Prepare output file
open(EXTRACT, ">".$pb_gb_type.".txt") || die "Cannot open file $pb_gb_type.txt.\n";
if(($pb_gb_type eq "pb0_gb1")||($pb_gb_type eq "pb0_gb2")||($pb_gb_type eq "pb0_gb5")
 ||($pb_gb_type eq "pb1_gb0")||($pb_gb_type eq "pb3_gb0")||($pb_gb_type eq "pb4_gb0")){
	print EXTRACT "Ligand\tELE\tVDW\tNP_SOLV\tP_SOLV\tE_TOT\n";
}
elsif($pb_gb_type eq "pb4_gb1"){
	print EXTRACT "Ligand\tELE\tVDW\tPB_NP_SOLV\tPB_P_SOLV\tPBTOT\tGB_NP_SOLV\tGB_P_SOLV\tGBTOT\n";
}
else{
	print "The type of MM-PBSA / MM-GBSA analysis option you specified in the program call\n";
	print "is not supported. Please refer to the AmberTools Manual for available input options.\n";
	unlink("$pb_gb_type.txt");
	exit;
}

foreach my $struct (@structs){

	my $out_file = $calc_folder."/".$struct."/s".$snaps_inter."/".$pb_gb_type."/".$struct."_statistics.out";
	$start = 0;
		
	open(OUT, $out_file) || die "Cannot open $out_file for reading.\n";
		
	while(my $o = <OUT>){
		chomp($o);
			
		if($o =~ m/DELTA/){
			$start = 1;
		}
			
		if($start == 1){
			my @o = split(/\s+/, $o);
				
			if($o =~ m/^ELE /){
				$ele = $o[1];
			}
			if($o =~ m/^VDW /){
				$vdw = $o[1];
			}
			if($o =~ m/^PBSUR /){
				$pbsur = $o[1];
			}
			if($o =~ m/^PBCAL /){
				$pbcal = $o[1];
			}
			if($o =~ m/^PBTOT /){
				$pbtot = $o[1];
			}
			if($o =~ m/^GB /){
				$gb = $o[1];
			}
			if($o =~ m/^GBSUR/){
				$gbsur = $o[1];
			}
			if($o =~ m/^GBTOT/){
				$gbtot = $o[1];
			}
		}
	}
	close OUT;
	
	# Write extracted data
	if(($pb_gb_type eq "pb0_gb1")||($pb_gb_type eq "pb0_gb2")||($pb_gb_type eq "pb0_gb5")){
		print EXTRACT "$struct\t$ele\t$vdw\t$gbsur\t$gb\t$gbtot\n";
	}
	elsif(($pb_gb_type eq "pb1_gb0")||($pb_gb_type eq "pb3_gb0")||($pb_gb_type eq "pb4_gb0")){
		print EXTRACT "$struct\t$ele\t$vdw\t$pbsur\t$pbcal\t$pbtot\n";
	}
	else{
		print EXTRACT "$struct\t$ele\t$vdw\t$pbsur\t$pbcal\t$pbtot\t$gbsur\t$gb\t$gbtot\n";
	}
}
close(EXTRACT);
