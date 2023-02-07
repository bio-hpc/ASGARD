#/usr/bin/perl

# perl-script for automated calculation of free energies of binding for ligands
# of membrane proteins using APBS and an implicit membrane slab model.

use FindBin;
use File::Copy;
use strict;
my $current_dir = $FindBin::Bin;

# Read input information from file mmpbsa_FEWmem.in
my %c_in;
open(IN, "mmpbsa_FEWmem.in") || die "Cannot open file mmpbsa_FEWmem.in for reading.\n";
while(my $in_line = <IN>){
	chomp($in_line);
	my @in_line = split(/\t/, $in_line);
	$c_in{$in_line[0]} = $in_line[1]; 
}
close IN;

# Run MM and SA calculation using mm_pbsa.pl
my $mmpbsa_c = "perl ".$c_in{'MMPBSAexe'}." mmpbsa.in > mmpbsa.log 2> mmpbsa.log";
call_prog($mmpbsa_c, "mmpbsa.log", $c_in{'CALCdir'});


# Prepare APBS input files for the generation of maps for receptor and complex and
# for the calculation of the solvation free energy of complex, receptor, and ligand.
copy($c_in{'DUMMY'}, "apbs_dummy.in") || die "Cannot copy ".$c_in{'DUMMY'}." to apbs_dummy.in\n";
copy($c_in{'SOLV'}, "apbs_solv.in") || die "Cannot copy ".$c_in{'SOLV'}." to apbs_solv.in\n";


foreach my $spec ("dummy", "solv"){
	my $out_str = "";
	open(IN_FILE, "apbs_".$spec.".in") || die "Cannot open file apbs_".$spec.".in for reading.\n";
	while(my $f_line = <IN_FILE>){
		chomp($f_line);
		if($f_line =~ m/DIME/){
			$f_line =~ s/DIME/$c_in{'DIME'}/g;
		}
		if($f_line =~ m/GLEN_L/){
			$f_line =~ s/GLEN_L/$c_in{'GLEN_L'}/g;
		}
		if($f_line =~ m/GLEN_M/){
			$f_line =~ s/GLEN_M/$c_in{'GLEN_M'}/g;
		}
		if($f_line =~ m/GLEN_S/){
			$f_line =~ s/GLEN_S/$c_in{'GLEN_S'}/g;
		}
		if($f_line =~ m/PDIE/){
			$f_line =~ s/PDIE/$c_in{'PDIE'}/g;
		}
		if($f_line =~ m/SDIE/){
			$f_line =~ s/SDIE/$c_in{'SDIE'}/g;
		}
		if($f_line =~ m/IONCONC/){
			$f_line =~ s/IONCONC/$c_in{'IONCONC'}/g;
		}
		if(($c_in{'DEC'} == 1)&&($f_line =~ m/calcenergy total/)){
			$f_line =~ s/total/comps/;
		}
		$out_str .= $f_line ."\n";
	}
	close(IN_FILE);

	unlink("apbs_".$spec.".in");

	open(NEW_IN_FILE, ">apbs_".$spec.".in") || die "Cannot open file apbs_".$spec.".in for writing.\n";
	print NEW_IN_FILE $out_str;
	close NEW_IN_FILE;
}


# APBS calculation
my %tdc_solv_ener;
my %bdc_solv_ener;
my %sdc_solv_ener;

my %tdc_solv_part;
my %tdc_ref_part;
my %bdc_solv_part;
my %bdc_ref_part;
my %sdc_solv_part;
my %sdc_ref_part;
my $kj2kt = 1000.0 / (300.0 * 8.314);
my $kj2kcal = 1/4.184;
my $read_pqr_file;
my %dec_assoc;
my @backbone_atoms = ("CA", "C", "O", "HA", "N", "H");

# Perform APBS calculation for complex and receptor
foreach my $tag ("com", "rec", "lig"){

	$read_pqr_file = 0; # Flag for reading first PQR file in case DEC=1
	my $read_header = 1; # Flag for reading header of *.all.out files
	my $all_out_file = $c_in{'STRUCT'}."_".$tag.".all.out";	
	open(NEW_ALL_OUT, ">".$all_out_file."_new") || die "Cannot open file ".$all_out_file."_new for writing.\n";
	
	for(my $i=$c_in{'START'}; $i<=$c_in{'STOP'}; $i+=$c_in{'OFFSET'}){

		# Prepare next round
		unlink "Current.pqr" if(-l "Current.pqr");
		unlink "Current.pqr" if(-e "Current.pqr");
		unlink "Receptor.pqr" if(-l "Receptor.pqr");
		unlink "Receptor.pqr" if(-e "Receptor.pqr");
		
		my $new_all_out_str = "";
		my $new_BDC_str = "";
		my $new_SDC_str = "";

		# Define current input PQR structures
		my $current = "pqr_snaps/".$c_in{'STRUCT'}."_".$tag.".pqr.".$i;
		symlink($current, "Current.pqr");
		
		# Define reference receptor structure for grid centering
		my $rec_reference = "pqr_snaps/".$c_in{'STRUCT'}."_rec.pqr.".$i;
		symlink($rec_reference, "Receptor.pqr");

		# Generate maps
		my $maps_com_rec_c = $c_in{'APBSexe'}." apbs_dummy.in > apbs_dummy.log 2> apbs_dummy.log";
		call_prog($maps_com_rec_c, "apbs_dummy.log", $c_in{'CALCdir'});
		if($c_in{'VERBOSE'} != 1){
			unlink("apbs_dummy.log");
		}
		else{
			rename("apbs_dummy.log", "apbs_dummy_".$tag.".log.".$i);
		}

		# Define map_tags
		my @map_tags;
		if($c_in{'FOCUS'} == 0){
			@map_tags = ("dielx_S", "diely_S", "dielz_S", "kappa_S", "charge_S");
		}
		else{
			@map_tags = ("dielx_L", "dielx_M", "dielx_S", "diely_L", "diely_M", "diely_S", "dielz_L", "dielz_M",
			"dielz_S", "kappa_L", "kappa_M", "kappa_S", "charge_L", "charge_M", "charge_S");
		}

		# For receptor and complex incorporate membrane into maps
		if($tag ne "lig"){
			if($c_in{'FOCUS'} == 0){
				foreach my $size ("S"){
					# Insert membrane
					my $incorporate_c = $c_in{'DRAWmem'}." dielx_".$size.".dx ".$c_in{'ZMEM'}." ".$c_in{'TMEM'}." ".$c_in{'PDIE'}." 0.0 ".$c_in{'IONCONC'}." ".$c_in{'R_TOP'}." ".$c_in{'R_BOTTOM'};
					$incorporate_c .= " >draw_membrane.log 2>draw_membrane.log";
					call_prog($incorporate_c, "draw_membrane.log", $c_in{'CALCdir'});
					unlink("draw_membrane.log");
				}
			}
			else{
				foreach my $size ("L", "M", "S"){
					# Insert membrane
					my $incorporate_c = $c_in{'DRAWmem'}." dielx_".$size.".dx ".$c_in{'ZMEM'}." ".$c_in{'TMEM'}." ".$c_in{'PDIE'}." 0.0 ".$c_in{'IONCONC'}." ".$c_in{'R_TOP'}." ".$c_in{'R_BOTTOM'};
					$incorporate_c .= " >draw_membrane.log 2>draw_membrane.log";
					call_prog($incorporate_c, "draw_membrane.log", $c_in{'CALCdir'});
					unlink("draw_membrane.log");
				}
			}
				
			# Insert second dielectric region, if dielectric constant of the membrane differs from 2 
			# or a second slab region was defined. Third slab region might be defined in addition.
			if(($c_in{'DIELC_SEC_SLAB'} != 0)&&($c_in{'T_SEC_SLAB'} != 0)){
				modify_membrane($c_in{'DIELC_MEM'}, $c_in{'T_SEC_SLAB'}, $c_in{'DIELC_SEC_SLAB'}, $c_in{'T_THIRD_SLAB'}, $c_in{'DIELC_THIRD_SLAB'},  $c_in{'ZMEM'}, $c_in{'TMEM'}, \@map_tags);
			}
		}
		
		# For ligand use created dx-files without implicit membrane representation
		if($tag eq "lig"){
			foreach my $map_tag (@map_tags){
				rename($map_tag.".dx", $map_tag."m.dx");
			}
		}

		# Compute solvation free energy
		my $solv_c = $c_in{'APBSexe'}." apbs_solv.in > apbs_solv.log 2> apbs_solv.log";
		my $solv_c = $c_in{'APBSexe'}." --output-file=apbs_solv.out apbs_solv.in > apbs_solv.log 2> apbs_solv.log";
		call_prog($solv_c, "apbs_solv.log", $c_in{'CALCdir'});


		# Analyze results
		
		# 1.) of decomposition calculation
		if($c_in{'DEC'} == 1){
			
			# Read first PQR file
			my %total_resno;
			if($read_pqr_file == 0){
				my $pqr_file_IDs = "pqr_snaps/".$c_in{'STRUCT'}."_".$tag.".pqr.".$i;
				my $atom_id = 0;
				open(PQR_ID, $pqr_file_IDs) || die "Cannot open file $pqr_file_IDs for reading.\n";
				while(my $pqrID_l = <PQR_ID>){
					if(($pqrID_l =~ m/ATOM/)||($pqrID_l =~ m/HETATM/)){
						my($card, $atnum, $atom, $alt, $resname, $resno,   $x, $y, $z) =
						unpack("a6    a5  x   a4     a     a3    x2  a4   x4   a8  a8  a8", $pqrID_l);
		
						$resno =~ s/^\s+//g;
						$resno =~ s/\s+$//g;
						$resname =~ s/^\s+//g;
						$resname =~ s/\s+$//g;
						$atnum =~ s/^\s+//g;
						$atnum =~ s/\s+$//g;
						$atom =~ s/^\s+//g;
						$atom =~ s/\s+$//g;
						
						if(($tag eq "com")||($tag eq "rec")){
							$atom_id = $atnum - 1; # Since atom numbering starts from 0 in apbs-output file
							$dec_assoc{$tag}->{$atom_id}->{'resno'} = $resno;
						}
						if($tag eq "lig"){
							$dec_assoc{$tag}->{$atom_id}->{'resno'} = 1;
							$atom_id++;
						}
						$dec_assoc{$tag}->{$atom_id}->{'atomname'} = $atom;
						$dec_assoc{$tag}->{$atom_id}->{'resname'} = $resname;
						$total_resno{$tag} = $resno;
					}
				}
			}
			
			# Read APBS output
			my $id;
			open(APBS_SOLV_OUT, "apbs_solv.out") || die "Cannot open file apbs_solv.out for reading.\n";
			while(my $apbs_solv_out = <APBS_SOLV_OUT>){
				chomp($apbs_solv_out);
				
				# Detect ID and determine if section shall be read or not.
				if($apbs_solv_out =~ m/^\s+id/){
					$apbs_solv_out =~ s/^\s+//g;
					$apbs_solv_out =~ s/\s+$//g;
					my @id = split(/\s+/, $apbs_solv_out);
					$id = $id[1];
				}
				
				# Read atom entries of solvation energy calculation part
				if((($c_in{'FOCUS'} == 0)&&($id == 1))||(($c_in{'FOCUS'} == 1)&&($id == 3))){
					if($apbs_solv_out =~ m/^\s+atom/){
						my @atom_id = split(/\s+/, $apbs_solv_out);
						
						# TDC
						$tdc_solv_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} = $tdc_solv_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} + $atom_id[3];
						
						# BDC
						if((grep {$_ eq $dec_assoc{$tag}->{$atom_id[2]}->{'atomname'}} @backbone_atoms)
						&&(($dec_assoc{$tag}->{$atom_id[2]}->{'resname'} ne "ACE")
						 &&($dec_assoc{$tag}->{$atom_id[2]}->{'resname'} ne "NME"))){
							$bdc_solv_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} = $bdc_solv_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} + $atom_id[3];
						}
						# SDC
						else{
							$sdc_solv_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} = $sdc_solv_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} + $atom_id[3];
						}
					}
				}
					
				# Read atom entries of reference energy calculation part
				elsif((($c_in{'FOCUS'} == 0)&&($id == 2))||(($c_in{'FOCUS'} == 1)&&($id == 6))){
					if($apbs_solv_out =~ m/^\s+atom/){
						my @atom_id = split(/\s+/, $apbs_solv_out);
						
						# TDC
						$tdc_ref_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} = $tdc_ref_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} + $atom_id[3];
					
						# BDC
						if((grep {$_ eq $dec_assoc{$tag}->{$atom_id[2]}->{'atomname'}} @backbone_atoms)
						 &&(($dec_assoc{$tag}->{$atom_id[2]}->{'resname'} ne "ACE")
						 &&($dec_assoc{$tag}->{$atom_id[2]}->{'resname'} ne "NME"))){
							$bdc_ref_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} = $bdc_ref_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} + $atom_id[3];
						}
						# SDC
						else{
							$sdc_ref_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} = $sdc_ref_part{$tag}->{$i}->{$dec_assoc{$tag}->{$atom_id[2]}->{'resno'}} + $atom_id[3];
						}
					}
				}
				
				# For all other IDs do nothing.
				else{
					next;
				}
			}
			close APBS_SOLV_OUT;
			
			# Calculate solvation free energy per residue
			for(my $resno = 0; $resno <= $total_resno{$tag}; $resno++){
				# TDC
				$tdc_solv_ener{$tag}->{$i}->{$resno} = $tdc_solv_part{$tag}->{$i}->{$resno} - $tdc_ref_part{$tag}->{$i}->{$resno};
				$tdc_solv_ener{$tag}->{$i}->{$resno} = $tdc_solv_ener{$tag}->{$i}->{$resno} * $kj2kcal;
				$tdc_solv_ener{$tag}->{$i}->{$resno} = sprintf("%.3f", $tdc_solv_ener{$tag}->{$i}->{$resno});
				# BDC
				$bdc_solv_ener{$tag}->{$i}->{$resno} = $bdc_solv_part{$tag}->{$i}->{$resno} - $bdc_ref_part{$tag}->{$i}->{$resno};
				$bdc_solv_ener{$tag}->{$i}->{$resno} = $bdc_solv_ener{$tag}->{$i}->{$resno} * $kj2kcal;
				$bdc_solv_ener{$tag}->{$i}->{$resno} = sprintf("%.3f", $bdc_solv_ener{$tag}->{$i}->{$resno});
				# SDC
				$sdc_solv_ener{$tag}->{$i}->{$resno} = $sdc_solv_part{$tag}->{$i}->{$resno} - $sdc_ref_part{$tag}->{$i}->{$resno};
				$sdc_solv_ener{$tag}->{$i}->{$resno} = $sdc_solv_ener{$tag}->{$i}->{$resno} * $kj2kcal;
				$sdc_solv_ener{$tag}->{$i}->{$resno} = sprintf("%.3f", $sdc_solv_ener{$tag}->{$i}->{$resno});
			}
			
			# Read current entry in *<tag>.all.out file and generated new entry
			my $start_reading_all_out = 0;
			my $snap_no = 0;
			
			open(ALL_OUT, $all_out_file) || die "Cannot open file ".$all_out_file."for reading.\n";
			while(my $all_out = <ALL_OUT>){
				chomp($all_out);
				if(($all_out =~ m/MM/)&&($read_header == 1)){
					$new_all_out_str .= $all_out."\n";
					$new_all_out_str .= "GB\n";
					$new_all_out_str .= "MS\n";
					$new_all_out_str .= "GB_SURFTEN 0.072\n";
					$new_all_out_str .= "GB_SURFOFF 0.00\n";
					$read_header = 0;
					next;
				}
				if($all_out =~ m/^(\d+)/){
					$snap_no = $1;
					if($snap_no == $i){
						$start_reading_all_out = 1;
						$new_all_out_str .= $snap_no."\n";
					}
					elsif($snap_no < $i){
						$start_reading_all_out = 0;
					}
					else{
						last;
					}
				}
				if(($start_reading_all_out == 1)&&($all_out =~ m/^TDC/)){
					my @TDC = split(/\s+/, $all_out);
					my $PB_TDC = sprintf("%3s%7s%10s%10s%10s%10s%10s\n", "TDC", $TDC[1], $TDC[2], $TDC[3], $TDC[4],
					                     $tdc_solv_ener{$tag}->{$i}->{$TDC[1]}, $TDC[6]);
					$new_all_out_str .= $PB_TDC;
				}
				if(($start_reading_all_out == 1)&&($all_out =~ m/^BDC/)){
					my @BDC = split(/\s+/, $all_out);
					my $PB_BDC = sprintf("%3s%7s%10s%10s%10s%10s%10s\n", "BDC", $BDC[1], $BDC[2], $BDC[3], $BDC[4],
					                     $bdc_solv_ener{$tag}->{$i}->{$BDC[1]}, $BDC[6]);
					$new_BDC_str .= $PB_BDC;
				}
				if(($start_reading_all_out == 1)&&($all_out =~ m/^SDC/)){
					my @SDC = split(/\s+/, $all_out);
					my $PB_SDC = sprintf("%3s%7s%10s%10s%10s%10s%10s\n", "SDC", $SDC[1], $SDC[2], $SDC[3], $SDC[4],
					                     $sdc_solv_ener{$tag}->{$i}->{$SDC[1]}, $SDC[6]);
					$new_SDC_str .= $PB_SDC;
				}
			}
			close ALL_OUT;
			$new_all_out_str .= $new_SDC_str;
			$new_all_out_str .= $new_BDC_str;
		}
		
		# 2.) Analyze results of normal calculation (no decomposition)
		else{
			# Read output data
			open(APBS_SOLV_OUT, "apbs_solv.log") || die "Cannot open file apbs_solv.log for reading.\n";
			while(my $apbs_solv_out = <APBS_SOLV_OUT>){
				chomp $apbs_solv_out;
				if($apbs_solv_out =~ m/Global net ELEC energy = (-?\d+\.\d+E(\+|-)\d+) kJ\/mol/){
					my $ener = $1;
					$ener = $ener * $kj2kt;
					$ener = sprintf("%.6f", $ener);
					$tdc_solv_ener{$tag}->{$i} = $ener;
					last;
				}
			}
			close APBS_SOLV_OUT;
			
			my $line_association = 0;
			my $snap_no = 0;
			
			open(ALL_OUT, $all_out_file) || die "Cannot open file ".$all_out_file."for reading.\n";	
			while(my $all_out = <ALL_OUT>){
				chomp($all_out);
				if(($all_out =~ m/MM/)&&($read_header == 1)){
					$new_all_out_str .= $all_out."\n";
					$new_all_out_str .= "PB\n";
					$new_all_out_str .= "MS\n";
					$new_all_out_str .= "PB_SURFTEN 0.00542\n";
					$new_all_out_str .= "PB_SURFOFF 0.92\n";
					$read_header = 0;
					next;
				}		
				if($all_out =~ m/^(\d+)/){
					$snap_no = $1;
					if($snap_no == $i){
						$line_association = 1;
						$new_all_out_str .= $snap_no."\n";
					}
					elsif($snap_no < $i){
						$line_association = 0;
					}
					else{
						last;
					}
				}
				if(($line_association == 1)&&($all_out =~ m/^ BOND/)){
					$line_association = 2;
					$new_all_out_str .= $all_out."\n";
				}
				if(($line_association == 2)&&($all_out =~ m/^ VDWAALS/)){
					$line_association = 3;
					$new_all_out_str .= $all_out."\n";
				}
				if(($line_association == 3)&&($all_out =~ m/^ 1-4 VDW/)){
					$line_association = 4;
					$new_all_out_str .= $all_out."\n";
					$new_all_out_str .= "corrected reaction field energy:";
					$new_all_out_str .= sprintf("%16s\n", $tdc_solv_ener{$tag}->{$snap_no});
					next;
				}
				if(($line_association == 4)&&($all_out =~ m/^surface area =\s+(\d+\.\d+)/)){
					my $surf_area = $1;
					$new_all_out_str .= $all_out."\n";
					$new_all_out_str .= "ECAVITY =";
					$new_all_out_str .= sprintf("%11s\n", $surf_area);
					$new_all_out_str .= "EDISPER = 0.0000\n";
					$line_association = 5;
					next;
				}
			}
			close ALL_OUT;
		}

		print NEW_ALL_OUT $new_all_out_str;

		# Save output of solvation free energy calculation if verbose output was requested
		if($c_in{'VERBOSE'} != 1){
			# Clean up
			unlink("apbs_solv.log");
 			unlink glob "*.dx";
			unlink "io.mc";
		}
		else{
			rename("apbs_solv.log", "apbs_solv_".$tag.".log.".$i);
			rename("apbs_solv.out", "apbs_solv_".$tag.".out.".$i);
			foreach my $map_tag (@map_tags){
				rename($map_tag."m.dx", $map_tag."_".$tag."_".$i.".dx");
			}
		}
	}
	close NEW_ALL_OUT;
	unlink($all_out_file);
	rename($all_out_file."_new", $all_out_file);
}

# Clean up
unlink "Current.pqr" if(-l "Current.pqr");
unlink "Current.pqr" if(-e "Current.pqr");
unlink "Receptor.pqr" if(-l "Receptor.pqr");
unlink "Receptor.pqr" if(-e "Receptor.pqr");
unlink($c_in{'STRUCT'}."_statistics.out");
unlink($c_in{'STRUCT'}."_statistics.out.snap");

# Perform MM-PBSA statistics analysis
my $statistics_c = "";
if($c_in{'DEC'} == 0){
	$statistics_c = "perl $c_in{'MMPBSAstat'} 1 0 ".$c_in{'STRUCT'}."_statistics.in ".$c_in{'STRUCT'}."_statistics.out";
}
else{
	$statistics_c = "perl $c_in{'MMPBSAstat'} 1 1 ".$c_in{'STRUCT'}."_statistics.in ".$c_in{'STRUCT'}."_statistics.out";
}
$statistics_c .= " > ".$c_in{'STRUCT'}."_statistics.log 2> ".$c_in{'STRUCT'}."_statistics.log";
call_prog($statistics_c, $c_in{'STRUCT'}."_statistics.log", $c_in{'CALCdir'});


######################################################################################
# Subroutines

# Calling programs
sub call_prog{
	my $command = shift;
    	my $log_file = shift;
	my $calc_dir = shift;

	if($command =~ /draw_membrane/){
		system($command);
	}
	else{
		system($command) == 0
		or do { copy($log_file, $calc_dir."/".$log_file); die "Execution of external program failed. Please make sure
that all required programs are correctly installed on your system and
consult the file $calc_dir/$log_file for further information.\n"; };
    		if($? != 0){
			exit;
    		}
	}
}


# Modify dielectric constants of membrane region
sub modify_membrane{
	my $dielc_mem = shift;
	my $second_slab_thickness = shift;
	my $second_slab_dielc = shift;
	my $third_slab_thickness = shift;
	my $third_slab_dielc = shift;
	my $membrane_boundary = shift;
	my $membrane_thickness = shift;
	my $map_tags_ref = shift;	
	my @map_tags = @{$map_tags_ref};
	my %upper_sec_slab_region;
	my %lower_sec_slab_region;
	my %upper_third_slab_region;
	my %lower_third_slab_region;
	my %membrane_region;
	
	# Definition of upper and lower headgroup regions
	$upper_sec_slab_region{'top'} = $membrane_boundary + $membrane_thickness;
	$upper_sec_slab_region{'bottom'} = $membrane_boundary + $membrane_thickness - $second_slab_thickness;
	$lower_sec_slab_region{'top'} = $membrane_boundary + $second_slab_thickness;
	$lower_sec_slab_region{'bottom'} = $membrane_boundary;
	$membrane_region{'top'} = $membrane_boundary + $membrane_thickness - $second_slab_thickness;
	$membrane_region{'bottom'} = $membrane_boundary + $second_slab_thickness;
	
	if(($third_slab_dielc != 0)&&($third_slab_thickness != 0)){
		$upper_third_slab_region{'top'} = $membrane_boundary + $membrane_thickness - $second_slab_thickness;
		$upper_third_slab_region{'bottom'} = $membrane_boundary + $membrane_thickness - $second_slab_thickness - $third_slab_thickness;
		$lower_third_slab_region{'top'} = $membrane_boundary + $second_slab_thickness + $third_slab_thickness;
		$lower_third_slab_region{'bottom'} = $membrane_boundary + $second_slab_thickness;
		$membrane_region{'top'} = $membrane_boundary + $membrane_thickness - $second_slab_thickness - $third_slab_thickness;
		$membrane_region{'bottom'} = $membrane_boundary + $second_slab_thickness + $third_slab_thickness;
	}
	
	$dielc_mem = sprintf("%E", $dielc_mem);
	$second_slab_dielc = sprintf("%E", $second_slab_dielc);
	$third_slab_dielc = sprintf("%E", $third_slab_dielc);
	
	
	# Read initial dx-file and dx-file with membrane 
	foreach my $map_tag (@map_tags){
	
		if($map_tag =~ m/diel/){

			# Read initial map
			my %epsilon_org;
			my $start_org = 0;
			my $dot_count_org = 0;
		
			my ($dime_x, $dime_y, $dime_z);
			my ($origin_x, $origin_y, $origin_z);
			my ($inter_x, $inter_y, $inter_z);
			my $start_mem = 0;	
			my %current_coords;
			my %coords_mem;
			my %epsilon_mem;
			my $dot_count_mem = 0;
			my $out_str = "";
		

			open(DX_ORG, $map_tag.".dx") || die "Cannot open file ".$map_tag.".dx for reading.\n";
			open(DX_MEM, $map_tag."m.dx") || die "Cannot open file ".$map_tag."m.dx for reading.\n";
			while(! eof(DX_ORG)){
		
				# Read line of initial (original) dx-file
				my $dx_org_line = <DX_ORG>;
			
				chomp($dx_org_line);
				$dx_org_line =~ s/^\s+//g;
				$dx_org_line =~ s/\s+$//g;
	
				if($dx_org_line =~ m/^attribute/){
					$start_org = 0;
				}
	
				if($start_org == 1){
					my @dot_org_data = split(/\s+/, $dx_org_line);
					my $temp_dot_count = 0;
		
					foreach my $d (@dot_org_data){
						$dot_count_org++;
						$temp_dot_count++;
					
						# Assign epsilon
						$epsilon_org{$temp_dot_count} = sprintf("%.6f", $d);
					}
				}
	
				if($dx_org_line =~ m/data follows/){
					$start_org = 1;
				}
			
			
				# Read line of dx-file with membrane
				my $dx_mem_line = <DX_MEM>;
			
				chomp($dx_mem_line);
				$dx_mem_line =~ s/^\s+//g;
				$dx_mem_line =~ s/\s+$//g;
	
				if($dx_mem_line =~ m/object 1 class gridpositions counts\s+(\d+)\s+(\d+)\s+(\d+)/){
					$dime_x = $1;
					$dime_y = $2;
					$dime_z = $3;
				}
	
				if($dx_mem_line =~ m/^origin/){
					my @origin = split(/\s+/, $dx_mem_line);
					$origin_x = sprintf("%.6f", $origin[1]);
					$origin_y = sprintf("%.6f", $origin[2]);
					$origin_z = sprintf("%.6f", $origin[3]);
		
					# Initialization of current coordinates
					$current_coords{'x'} = $origin_x;
					$current_coords{'y'} = $origin_y;
					$current_coords{'z'} = $origin_z;
				}
	
				if($dx_mem_line =~ m/^delta/){
					my @inter = split(/\s+/, $dx_mem_line);
					$inter[1] = sprintf("%.6f", $inter[1]);
					$inter[2] = sprintf("%.6f", $inter[2]);
					$inter[3] = sprintf("%.6f", $inter[3]);
					if($inter[1] != 0){
						$inter_x = $inter[1];
					}
					if($inter[2] != 0){
						$inter_y = $inter[2];
					}
					if($inter[3] != 0){
						$inter_z = $inter[3];
					}
				}
	
				if($dx_mem_line =~ m/^attribute/){
					$start_mem = 0;
				}
	
				if($start_mem == 1){
					my @dot_mem_data = split(/\s+/, $dx_mem_line);
					my $dot_str = "";
					my $temp_dot_count = 0;
				
					foreach my $d (@dot_mem_data){
						$dot_count_mem++;
						$temp_dot_count++;

						if(($dot_count_mem % $dime_z != 1)&&(($dot_count_mem % ($dime_z * $dime_y)) != 1)){
							$current_coords{'z'} = $current_coords{'z'} + $inter_z;
						}
						if(($dot_count_mem % $dime_z == 1)&&(($dot_count_mem % ($dime_z * $dime_y)) != 1)){
							$current_coords{'y'} = $current_coords{'y'} + $inter_y;
							$current_coords{'z'} = $origin_z;
						}
						if(($dot_count_mem % $dime_z == 1)&&(($dot_count_mem % ($dime_z * $dime_y)) == 1)){
							$current_coords{'x'} = $current_coords{'x'} + $inter_x;
							$current_coords{'y'} = $origin_y;
							$current_coords{'z'} = $origin_z;
						}
			
						# Assign coordinates
						$coords_mem{$temp_dot_count}->{'x'} = $current_coords{'x'};
						$coords_mem{$temp_dot_count}->{'y'} = $current_coords{'y'};
						$coords_mem{$temp_dot_count}->{'z'} = $current_coords{'z'};
			
						# Assign epsilon
						$epsilon_mem{$temp_dot_count} = sprintf("%.6f", $d);
					
						# Check whether dot lies in second slab region
						if(($second_slab_dielc != 0)&&($second_slab_thickness != 0)){
						
							if((($coords_mem{$temp_dot_count}->{'z'} <= $upper_sec_slab_region{'top'})&&($coords_mem{$temp_dot_count}->{'z'} >= $upper_sec_slab_region{'bottom'}))
							||(($coords_mem{$temp_dot_count}->{'z'} <= $lower_sec_slab_region{'top'})&&($coords_mem{$temp_dot_count}->{'z'} >= $lower_sec_slab_region{'bottom'}))){
					
								# Check whether dielc other than 80 was assigned in initial (original) maps
								#print $epsilon_org{$dot_count_mem}."\n";
								if(($epsilon_org{$temp_dot_count} == 80.000000)&&($epsilon_mem{$temp_dot_count} == 2.000000)){
									$d = $second_slab_dielc;
								}
							}
						}
						
						if(($third_slab_dielc != 0)&&($third_slab_thickness != 0)){
							if((($coords_mem{$temp_dot_count}->{'z'} < $upper_third_slab_region{'top'})&&($coords_mem{$temp_dot_count}->{'z'} >= $upper_third_slab_region{'bottom'}))
							||(($coords_mem{$temp_dot_count}->{'z'} <= $lower_third_slab_region{'top'})&&($coords_mem{$temp_dot_count}->{'z'} > $lower_third_slab_region{'bottom'}))){
								if(($epsilon_org{$temp_dot_count} == 80.000000)&&($epsilon_mem{$temp_dot_count} == 2.000000)){
									$d = $third_slab_dielc;
								}
							}
						}
						
						if(($dielc_mem != 0)&&($dielc_mem != 2)){
							if(($coords_mem{$temp_dot_count}->{'z'} < $membrane_region{'top'})&&($coords_mem{$temp_dot_count}->{'z'} > $membrane_region{'bottom'})){
								if(($epsilon_org{$temp_dot_count} == 80.000000)&&($epsilon_mem{$temp_dot_count} == 2.000000)){
									$d = $dielc_mem;
								}
							}
						}
						
						$dot_str .= $d." ";
					}
					$out_str .= $dot_str."\n";	
				}
				else{
					$out_str .= $dx_mem_line."\n";
				}
	
				if($dx_mem_line =~ m/data follows/){
					$start_mem = 1;
				}
			}
			close(DX_ORG);
			close(DX_MEM);
	
		
			if($dot_count_org != $dot_count_mem){
				print "ERROR: Different number of dots found in initial dx-file and in dx-file with membrane.\n";
				exit;
			}
		
			# Write new dx-file
			open(DX_HEAD, ">".$map_tag."withHead.dx") || die "Cannot open dx-file ".$map_tag."withHead.dx for writing.\n";
			print DX_HEAD $out_str;
			close DX_HEAD;
			
			unlink($map_tag."m.dx");
			rename($map_tag."withHead.dx", $map_tag."m.dx");
		}
	}
}
