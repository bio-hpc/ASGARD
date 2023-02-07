#/usr/bin/perl -w

# Program for conducting LIE calculations within the FEW workflow
# Input parameters: - Command-file
#                   - Ligand structure name
#

use File::Copy;
use File::Path;

my $command_file = $ARGV[0];
my $structure = $ARGV[1];
my $sander_exe;

# Read input parameters
if(@ARGV < 2){
  print "Wrong call of LIE.pl\n";
  print "Usage: LIE.pl <command file> <ligand structure name>\n";
  exit;
}

my %c_in;
my @traj_files;

open(C_IN, $command_file) || die "Cannot open file ".$command_file." for reading.\n";
while(my $c_in = <C_IN>){
  chomp($c_in);
  my $check_comment = substr($c_in, 0, 1);
  if(($c_in ne "")&&($check_comment ne "#")){
    my @c_in = split(/\s+/, $c_in);
    if($c_in[0] eq "parallel_call"){
	my $call_p = $c_in[1];
	for(my $p=2; $p<@c_in; $p++){
		$call_p .= " ".$c_in[$p];
	}
	$c_in[1] = $call_p;
    }
    if(($c_in[0] eq "trj_file")||($c_in[0] eq "trajectory_files")){
      push(@traj_files, $c_in[1]);
      next;
    }
    $c_in{$c_in[0]} = $c_in[1];
  }
}
close C_IN;

# Association between new and old keywords
my %new_to_old_keywords;

%new_to_old_keywords = (
		"output_path"  =>  "root_path",
		"no_of_rec_residues"  =>  "rec_res",
		"charge_method"  =>  "chrg_meth",
		"snaps_per_trajectory"  =>  "trj_snaps",
		"first_lie_snapshot"  =>  "start",
		"last_lie_snapshot"  =>  "stop",
		"offset_lie_snapshots"  =>  "offset",
		"sander_executable"  =>  "sander_exe",
		"parallel_lie_call"  =>  "parallel_call",
	);
	
	# Reset keywords
	foreach my $c_in_key (keys %c_in){
		if(exists $new_to_old_keywords{$c_in_key}){
			$c_in{$new_to_old_keywords{$c_in_key}} = $c_in{$c_in_key};
			delete($c_in{$c_in_key});
		}
	}

# Setting up paths and system definitions
my $analyze_dir = $c_in{'root_path'}."/lie_".$c_in{'chrg_meth'};
my %lie; # Hash for storing the the calculated energy contributions
my $tag_dir; 

# Define total number of residues in complex
my $total_res = $c_in{'rec_res'} + 1;
	
# Set default executable if no sander executable is specified
my $amberhome =  $ENV{"AMBERHOME"};
my $do_parallel;
my $sander_path;
if(($c_in{'parallel_call'})&&($c_in{'parallel_call'} ne "")){
	$do_parallel = $c_in{'parallel_call'};
}
if($c_in{'sander_exe'} =~ m/AMBERHOME/){
	$c_in{'sander_exe'} =~ s/AMBERHOME/$amberhome/;
}
if(($c_in{'sander_exe'})&&(-e $c_in{'sander_exe'})){
	$sander_exe = $c_in{'sander_exe'};
	my @path_components = split(/\//, $sander_exe);
	for(my $c=0; $c<$#path_components; $c++){
		$sander_path .= "/".$path_components[$c];
	}
}
if((! $c_in{'sander_exe'})||(($c_in{'sander_exe'})&&(! -e $c_in{'sander_exe'}))){
	if(-e "$amberhome/bin/sander"){
	  	print "\nWARNING: No sander executable specified. Running energy analysis with default settings.\n";
		if($do_parallel ne ""){
			$sander_exe = "$amberhome/bin/sander.MPI";
		}
		else{
			$sander_exe = "$amberhome/bin/sander";
		}
	}
}	


for my $tag ("com", "lig"){
  
	# Define complex/ligand directory
	$tag_dir = $analyze_dir."/".$structure."/".$tag;

	# Define topology-directory
	my $topo_dir = $analyze_dir."/".$structure."/".$tag."/topo";

	# Define system dependent components
	my @comp = ();
	my %comp_str;
	my @trajectories;
   
	if($tag eq "lig"){
		print "Energy calculation for ligand of $structure\n";
		@comp = ("tot", "lig", "wat");
		$comp_str{'tot'} = "total system";
		$comp_str{'lig'} = "ligand sub-system";
		$comp_str{'wat'} = "water sub-system";
	}
	if($tag eq "com"){
		print "Energy calculation for complex of $structure\n";
		@comp = ("tot", "res", "lig");
		$comp_str{'tot'} = "total system";
		$comp_str{'lig'} = "ligand sub-system";
		$comp_str{'res'} = "'not ligand' sub-system";
	}
	
	
	# Energy calculation
	#####################
		
	# Generate energy calculation directory
	my $ener_dir;
    
	print "Running energy calculation:\n";
	
	foreach my $comp (@comp){
	
		$ener_dir= $tag_dir."/ener_".$comp;	
		if(! -e $ener_dir){
			mkdir($ener_dir);
		}
 		chdir($ener_dir);
  
		# Prepare output
		my $ele_out_file = $ener_dir."/ele_".$c_in{'start'}."_".$c_in{'stop'}."_".$c_in{'offset'}.".dat";
		my $vdw_out_file = $ener_dir."/vdw_".$c_in{'start'}."_".$c_in{'stop'}."_".$c_in{'offset'}.".dat";
				
		open(ELE, ">$ele_out_file") || die "Cannot open file $ele_out_file for writing.\n";
		print ELE "Snapshot\tEEL\t1-4 EEL\n";
		open(VDW, ">$vdw_out_file") || die "Cannot open file $vdw_out_file for writing.\n";
		print VDW "Snapshot\tVDW\t1-4 VDW\n";    
			
		# Ensure that when ligand calculation is done only one processor is used to avoid crash
		# in the case where the number of residues is smaller or equal to the number of processors
		if($comp eq "lig"){
			$do_parallel = "";
			if(-e "$sander_path/sander"){
				$sander_exe = "$sander_path/sander";
			}
			else{
				$sander_exe = "$amberhome/bin/sander";
			}
		}
				
		# Prepare Sander input and output for calculation of vdW and electrostatic energies
		my $sander_in = $ener_dir."/sander.in";
		open(SANDER_IN, ">$sander_in") || die "Cannot open file $sander_in for writing.\n";
		print SANDER_IN "Single point calculation (sander)\n";
		print SANDER_IN " &cntrl\n";
		print SANDER_IN "  imin   = 5,\n";
		print SANDER_IN "  maxcyc = 0,\n";
		print SANDER_IN "  ncyc   = 0,\n";
		print SANDER_IN "  igb    = 6,\n";
		print SANDER_IN "\n";
		print SANDER_IN "  ntf    = 1,\n";
		print SANDER_IN "  ntb    = 0,\n";
		print SANDER_IN "  cut    = 999.0,\n";
		print SANDER_IN "  nsnb   = 9999,\n";
		print SANDER_IN " &end\n";
		close(SANDER_IN);

		my $topo_in = "../topo/".$comp.".top";
		my $crd_in = "../topo/".$comp.".crd";
			
		# Generate sander input and sasa output file for sasa calculation
		my $sasa_sander_in;
		my $sasa_sander_out;
		my $sasa_out_file;

		if(($c_in{'calc_sasa'} == 1)&&($comp eq "lig")){
				
			$sasa_out_file = $ener_dir."/sasa_".$c_in{'start'}."_".$c_in{'stop'}."_".$c_in{'offset'}.".dat";
			open(SASA, ">$sasa_out_file") || die "Cannot open file $sasa_out_file for writing.\n";
			print SASA "Snapshot\tSASA\n";
			
			$sasa_sander_in = $ener_dir."/sasa_sander.in";
			$sasa_sander_out = $ener_dir."/sasa_sander.out";
			
			open(SASA_SANDER, ">$sasa_sander_in") || die "Cannot write sander input file for ICOSA calculation $sasa_sander_in.\n";
			print SASA_SANDER "Sander input file for SASA calculations according to the ICOSA method\n";
			print SASA_SANDER " &cntrl\n";
			print SASA_SANDER "  ntf    = 1,       ntb    = 0,       dielc  = 1.0,\n";
			print SASA_SANDER "  imin   = 5,       maxcyc = 0,       ncyc   = 0,\n";
			print SASA_SANDER "  cut    = 999.0,   nsnb   = 99999,\n";
			if($tag eq "com"){
				print SASA_SANDER "  idecomp= 2,\n";
			}
			print SASA_SANDER "  igb    = 2,       saltcon= 0.00,     offset = 0.09,\n";
			print SASA_SANDER "  intdiel= 1.0,     extdiel= 80.0,\n";
			print SASA_SANDER "  gbsa   = 2,       surften= 1.0,\n";
			print SASA_SANDER " &end\n";

			if($tag eq "com"){
				print SASA_SANDER "Residues considered as REC\n";
				print SASA_SANDER "RRES 1 ".$c_in{'rec_res'}."\n";
				print SASA_SANDER "END\n";
				print SASA_SANDER "Residues considered as LIG\n";
				print SASA_SANDER "LRES ".$total_res." ".$total_res."\n";
				print SASA_SANDER "END\n";
				print SASA_SANDER "Residues to print\n";
				print SASA_SANDER "RES ".$total_res." ".$total_res."\n";
				print SASA_SANDER "END\n";
				print SASA_SANDER "END\n";
			}
			close SASA_SANDER;
		}
		
		# Run calculation
		my $traj_count = 0;
		foreach my $t (@traj_files){
			my $first_snap = $traj_count * $c_in{'trj_snaps'} + 1;
			$traj_count++;
			my $last_snap = $traj_count * $c_in{'trj_snaps'};
				
			# When start snapshot in trajectory
			if(($c_in{'start'} <= $last_snap)&&($c_in{'stop'} >= $first_snap)){
    
				# Execute sander
				print "Energy analysis for trajectory ".$t." of ".$comp_str{$comp}."\n";
				my @traj_name_components = split(/\./, $t);
				my $traj_base = $traj_name_components[0];
				my $mdcrd_in = "../s_".$comp."/".$traj_base."_nobox.mdcrd";
				my $c_sander = "$do_parallel $sander_exe -O -i sander.in -o sander.out -p $topo_in -c $crd_in -y $mdcrd_in 2> sander.log";
				call_prog($c_sander, $ener_dir."/sander.log");
         
				# Extract results
				my $start_read = 0;
				my %vdw = 0;
				my %eel = 0;
				my %vdw1_4 = 0;
				my %eel1_4 = 0;
				my $frame_count = 0;
         
				open(SANDER_OUT, "sander.out") || die "Cannot open file ".$sander_out." for reading\n";
				while(my $s_line = <SANDER_OUT>){
					chomp($s_line);
					if(($start_read == 1)&&($s_line =~ m/^ VDWAALS/)){
						my @s_line = split(/\s+/, $s_line);
						$eel{$frame_count} = $s_line[6];
						$vdw{$frame_count} = $s_line[3];
					}
					if(($start_read == 1)&&($s_line =~ m/^ 1-4 VDW/)){
						my @s_line = split(/\s+/, $s_line);
						$eel1_4{$frame_count} = $s_line[8];
						$vdw1_4{$frame_count} = $s_line[4];
						$start_read = 0; 
					}
					if($s_line =~ m/   NSTEP/){
						$start_read = 1;
						$frame_count++;
					}
					if($frame_count > $c_in{'trj_snaps'}){
						print "\nERROR: More snapshots present in trajectory than specified in 'trj_snaps'.\n";
						print "Please ensure that the total number of snapshots per trajectory, specified in\n";
						print "the command file, is correct.\n";
						exit;
					}
				}
				close SANDER_OUT;
				
				for(my $i=1; $i<=$frame_count; $i++){
					my $snap_no = $c_in{'start'} - $c_in{'offset'} + ($i * $c_in{'offset'});
					if($snap_no > $c_in{'stop'}){
						last;
					}
					if($snap_no >= $c_in{'start'}){
						print ELE $snap_no."\t".$eel{$i}."\t".$eel1_4{$i}."\n";
						print VDW $snap_no."\t".$vdw{$i}."\t".$vdw1_4{$i}."\n";
					}
				}
				
				# Calculate solvent accessible surface area of free and bound ligand				
				if(($c_in{'calc_sasa'} == 1)&&($comp eq "lig")){
					# Calculate surface area
					my ($crd_sasa, $topo_sasa, $mdcrd_sasa);
					if($tag eq "lig"){
						$crd_sasa = $crd_in;
						$topo_sasa = $topo_in;
						$mdcrd_sasa = $mdcrd_in;
					}
					else{
						$crd_sasa = "../topo/com.crd";
						$topo_sasa = "../topo/com.top";
						$mdcrd_sasa = "../s_com/".$traj_base."_nobox.mdcrd";
					}
					my $c_sasa = "$do_parallel $sander_exe -O -i sasa_sander.in -o sasa_sander.out -p $topo_sasa -c $crd_sasa -y $mdcrd_sasa";
					call_prog($c_sasa);
					
					# Read surface area from sander output file and save it to SASA output file
					my %sasa;
					my $frame_count_sasa = 0;
					
					open(MSLOG, $sasa_sander_out) || die "Cannot open molsurf log-file for reading.\n";
					while(my $m_line = <MSLOG>){
						chomp($m_line);
						if(($tag eq "com")&&($m_line =~ m/TDC\s+$total_res/)){
							$frame_count_sasa++;
							my @sasa = split(/\s+/, $m_line);
							$sasa{$frame_count_sasa} = $sasa[$#sasa];
						}
						if(($tag eq "lig")&&($m_line =~ m/^ ESURF\s+\=\s+(\d+)\.(\d+)/)){
							$frame_count_sasa++;
							$sasa{$frame_count_sasa} = $1.".".$2;
						}
						if($frame_count_sasa > $c_in{'trj_snaps'}){
							print "\nERROR: More snapshots present in trajectory than specified in 'trj_snaps'.\n";
							print "Please ensure that the total number of snapshots per trajectory, specified in\n";
							print "the command file, is correct.\n";
							exit;
						}
					}
					close MSLOG;
					
					for(my $i=1; $i<=$frame_count_sasa; $i++){
						my $snap_no = $c_in{'start'} -$c_in{'offset'} + ($i * $c_in{'offset'});
						if($snap_no > $c_in{'stop'}){
							last;
						}
						if($snap_no >= $c_in{'start'}){
							print SASA "$snap_no\t$sasa{$i}\n";
						}
					}	
				} # End SASA calculation
			} # End for snapshots >= start and <=end
		}  # End for each trajectory
		close ELE;
		close VDW;
		close SASA;
	} # End comp   
} # End tag


# Process results
my ($lig_ele_diff, $lig_vdw_diff, $lig_sasa, $com_ele_diff, $com_vdw_diff, $com_sasa);
my $results_file = $c_in{'root_path'}."/lie_".$c_in{'chrg_meth'}."/".$structure."/LIE_s".$c_in{'start'}."_".$c_in{'stop'}."_".$c_in{'offset'}.".txt";
open(RESULTS, ">$results_file") || die "Cannot open file $results_file for writing\n";
print RESULTS "Energy contributions:\n";
printf RESULTS "%-25s%13s%22s%18s\n", " ", "Mean energy", "Standard deviation", "Standard error";
  
for my $tag ("com", "lig"){
  
	# Define system dependent tags
	my @comp = ();
   
	if($tag eq "lig"){
		@comp = ("tot", "lig", "wat");
	}
	if($tag eq "com"){
		@comp = ("tot", "res", "lig");
	}
     	
	for my $ener ("vdw", "ele", "sasa"){
		foreach my $comp (@comp){
			
			# Ensure that SASA is only read if component is ligand
			if(($ener eq "sasa")&&($comp ne "lig")){
				next;
			}
				
			if(($ener eq "sasa")&&((! defined $c_in{'calc_sasa'})||($c_in{'calc_sasa'} == 0))){
				next;
			}
			
			# Read energy 
			###############	
			my $ener_dir = $c_in{'root_path'}."/lie_".$c_in{'chrg_meth'}."/".$structure."/".$tag."/ener_".$comp;
			my $ener_file = $ener_dir."/".$ener."_".$c_in{'start'}."_".$c_in{'stop'}."_".$c_in{'offset'}.".dat";
		  
			open(E_FILE, $ener_file) || die "Cannot open file $ener_file for reading.\n";
			my $header = <E_FILE>;
		  
			while(my $e_line = <E_FILE>){
				$snaps_count++;
				chomp $e_line;
				my @e_line = split(/\t/, $e_line);
				my $sum = $e_line[1]+$e_line[2];
				$lie{$tag}->{$comp}->{$ener}->{$e_line[0]} = $sum;
			}
			close E_FILE;
		}

		# Compute LIE contributions 
		############################
		my $ener_diff_file;
		my $sum_for_avg = 0;
		my $count = 0;
      	my %diff;
			
		# Calculate difference of energy contributions and determine average difference
		if($ener ne "sasa"){    
			$ener_diff_file = $c_in{'root_path'}."/lie_".$c_in{'chrg_meth'}."/".$structure."/".$tag."_".$ener."_s".$c_in{'start'}."_".$c_in{'stop'}."_".$c_in{'offset'}.".txt";
			
			open(ENER_DIFF, ">$ener_diff_file") || die "Cannot open file $ener_diff_file for writing\n";
		
			for(my $i=$c_in{'start'}; $i<=$c_in{'stop'}; $i+=$c_in{'offset'}){
				$count++;
				if($tag eq "lig"){  
					$diff{$i} = $lie{'lig'}->{'tot'}->{$ener}->{$i} - $lie{'lig'}->{'lig'}->{$ener}->{$i} - $lie{'lig'}->{'wat'}->{$ener}->{$i};
				}
				if($tag eq "com"){
					$diff{$i} = $lie{'com'}->{'tot'}->{$ener}->{$i} - $lie{'com'}->{'lig'}->{$ener}->{$i} - $lie{'com'}->{'res'}->{$ener}->{$i};
				}
				print ENER_DIFF $i."\t".$diff{$i}."\n";
				$sum_for_avg = $sum_for_avg + $diff{$i};
			}
	  		close ENER_DIFF;
		}
			
		# Determine average SASA
		if($ener eq "sasa"){
			for(my $i=$c_in{'start'}; $i<=$c_in{'stop'}; $i+=$c_in{'offset'}){
				$count++;
				$sum_for_avg = $sum_for_avg + $lie{$tag}->{'lig'}->{$ener}->{$i};
			}
		}
	  
		my $avg = $sum_for_avg / $count;
		$avg = sprintf("%.6f", $avg);
			
		# Determine standard error
		my $deviation = 0;
		for(my $i=$c_in{'start'}; $i<=$c_in{'stop'}; $i+=$c_in{'offset'}){
			if($ener ne "sasa"){
				$deviation = $deviation + (($diff{$i} - $avg)**2);
			}
			else{
				$deviation = $deviation + (($lie{$tag}->{'lig'}->{$ener}->{$i} - $avg)**2);
			}
		}
			
		my $total_snaps_no = $count;
			
		my $std = sqrt((1/($total_snaps_no-1))*$deviation);
		$std = sprintf("%.6f", $std);
		my $str_err = $std/sqrt($total_snaps_no);
		$str_err = sprintf("%.6f", $str_err);
	                                    
		if(($tag eq "lig")&&($ener eq "ele")){
			printf RESULTS "%-25s%13s%22s%18s\n", "Ligand non-bonded ELE: ", $avg, $std, $str_err;
			$lig_ele_diff = $avg;
		}
		if(($tag eq "lig")&&($ener eq "vdw")){
			printf RESULTS "%-25s%13s%22s%18s\n", "Ligand non-bonded VDW: ", $avg, $std, $str_err;
			$lig_vdw_diff = $avg;
		}
		if(($tag eq "lig")&&($ener eq "sasa")&&($c_in{'calc_sasa'} != 0)){
			printf RESULTS "%-25s%13s%22s%18s\n", "SASA of unbound ligand: ", $avg, $std, $str_err;
			$lig_sasa = $avg;
		}	 
		if(($tag eq "com")&&($ener eq "ele")){
			printf RESULTS "%-25s%13s%22s%18s\n", "Complex non-bonded ELE: ", $avg, $std, $str_err;
			$com_ele_diff = $avg;
		}
		if(($tag eq "com")&&($ener eq "vdw")){
			printf RESULTS "%-25s%13s%22s%18s\n", "Complex non-bonded VDW: ", $avg, $std, $str_err;
			$com_vdw_diff = $avg;
		}
		if(($tag eq "com")&&($ener eq "sasa")&&($c_in{'calc_sasa'} != 0)){
			printf RESULTS "%-25s%13s%22s%18s\n", "SASA of bound ligand: ", $avg, $std, $str_err;
			$com_sasa = $avg;
		}
	}
}
my $lie_ele = $com_ele_diff - $lig_ele_diff;
$lie_ele = sprintf("%.6f", $lie_ele);
my $lie_vdw = $com_vdw_diff - $lig_vdw_diff;
$lie_vdw = sprintf("%.6f", $lie_vdw);
my $lie_sasa = $com_sasa - $lig_sasa;
my $dG_bind = 0.16 * $lie_vdw + 0.5 * $lie_ele;
$dG_bind = sprintf("%.6f", $dG_bind);
	
print RESULTS "\n\n";
print RESULTS "================================================================================\n";
print RESULTS "Final LIE results:\n\n";
print RESULTS "Differences in interaction energies and SASA:\n\n";
printf RESULTS "%-38s%12s%9s\n", "Delta electrostatic energy (dE-ele): ", $lie_ele, "kcal/mol";
 printf RESULTS "%-38s%12s%9s\n", "Delta van der Waals energy (dE-vdW): ", $lie_vdw, "kcal/mol";
if($lie_sasa != 0){
	printf RESULTS "%-38s%12s%9s\n", "SASA difference (dSASA): ", $lie_sasa, "Å^^2";
}
print RESULTS "________________________________________________________________________________\n";
print RESULTS "\n";
print RESULTS "Estimated binding free energy:\n";
print RESULTS "\n";
print RESULTS "ATTENTION: Please note that the estimate of the binding free energy (dG-bind)\n";
print RESULTS "           provided below was calculated according to the equation\n";
print RESULTS "           dG-bind = alpha * dE-vdW + beta * dE-ele\n";
print RESULTS "           with fixed coefficients of alpha=0.16 and beta=0.5.\n";
print RESULTS "           We strongly recommend to thoroughly investigate whether using other\n";
print RESULTS "           coefficients could be beneficial in the specific case.\n";
print RESULTS "\n";
printf RESULTS "%-38s%12s%9s\n", "Binding energy estimate (dG-bind): ", $dG_bind, "kcal/mol";
close RESULTS;


######################################################################################
# Subroutines

# Calling programs
sub call_prog{
        my $command = shift;
		my $log_file = shift;

        system($command) == 0
        or die "Execution of external program failed. Please make sure that
all required programs are correctly installed on your system and consult the 
file $log_file for further information.\n";
        if($? != 0){
                exit;
        }
}
