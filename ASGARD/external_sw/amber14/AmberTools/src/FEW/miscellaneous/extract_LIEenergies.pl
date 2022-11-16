#/usr/bin/perl

# Perl-Script for automated extraction of the calculated differences in electrostatic
# and vdW interaction energies of the ligand in the bound and free state from the output
# files of the LIE calculation for a set of ligands. It is assumed that LIE analyses were
# prepared with the LIEW module of FEW.
#
# Usage: perl extract_LIEenergies.pl <structure reference file> <path> <name of LIE output file>
#
# structure reference file: Text file containing the names of the ligands that shall be
#                           considered (one ligand name per line) and experimentally 
#                           measured IC50 or binding free energies in tab-separated format.
#                           Example:
#                           #Ligands	dG
#                           Lig_5	-0.5394
#                           Lig_17	-1.3409
#                           ....
# path: Directory in which LIE calculations were conducted, e.g. /home/work/FEW_test/lie_am1
# name of LIE output file: Name of final result file of LIE calculations, e.g. LIE_s301_500_1.txt
#################################################################################################### 

my $struct_file = $ARGV[0];
my $path = $ARGV[1];
my $lie_out_file = $ARGV[2];

if(($struct_file eq "")||($path eq "")||($lie_out_file eq "")){
	print "Usage: perl extract_energies.pl <structure reference file> <path> <name of LIE output file>\n";
	exit;
}

# Read structure info
my @structs;
my %pIC;

open(S_FILE, "$path/$struct_file") || die "Cannot open file ".$path."/".$struct_file." for reading.\n";
my $header = <S_FILE>;
my $e_tpye;
my @h_line;
if($header =~ m/^#/){
	@h_line = split(/\s+/, $h_line);
	if($h_line[1] ne ""){
		$e_type = $h_line[1];
	}
}
	
while(my $s_line = <S_FILE>){
	chomp($s_line);
	my @s_line = split(/\t/, $s_line);
	push(@structs, $s_line[0]);
	$pIC{$s_line[0]} = $s_line[1];
}
close S_FILE;

# Read data
my %LIE;

foreach my $s (@structs){
	my $lie_file = $path."/".$s."/".$lie_out_file;
	open(LIE_FILE, $lie_file) || die "Cannot open file $lie_file for reading.\n";
	while(my $l = <LIE_FILE>){
		chomp($l);
		my @l = split(/\:/, $l);
		if($l =~ m/^Delta electrostatic energy/){
			$l[1] =~ s/ +//;
			my @ELE = split(/\s+/, $l[1]);
			$LIE{$s}->{"ELE"} = $ELE[0];
		}
		if($l =~ m/^Delta van der Waals energy/){
			$l[1] =~ s/ +//;
			my @VDW = split(/\s+/, $l[1]);
			$LIE{$s}->{"VDW"} = $VDW[0];
		}
		if($l =~ m/^SASA difference/){
			$l[1] =~ s/ +//;
			my @SASA = split(/\s+/, $l[1]);
			$LIE{$s}->{"SASA"} = $SASA[0];
		}
	}	 
}

# Print data
my $out_file = $path."/LIE_results.txt";
open(OUT, ">$out_file") || die "Cannot open file $out_file for writing.\n";
if($e_type ne ""){
	print OUT "Structure\t".$e_type."\tELE\tVDW";
}
else{
	print OUT "Structure\tdG\tELE\tVDW";
}

if(exists $LIE{$structs[0]}->{"SASA"}){
	print OUT "\tSASA\n";
}
else{
	print OUT "\n";
}

foreach my $s (@structs){
	if(exists $LIE{$s}->{"SASA"}){
		print OUT $s."\t".$pIC{$s}."\t".$LIE{$s}->{"ELE"}."\t".$LIE{$s}->{"VDW"}."\t".$LIE{$s}->{"SASA"}."\n";
	}
	else{
		print OUT $s."\t".$pIC{$s}."\t".$LIE{$s}->{"ELE"}."\t".$LIE{$s}->{"VDW"}."\n";
	}
}
close OUT;
