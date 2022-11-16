#!/usr/bin/env bash
#
#	Genera fichero de equilibracion A
#	se le pasa la salida para crearlo, numSteeps de equilibrado y num de save step
#	$1 ruta deonde se guardara el fichero
#	$2 stepEquilibradoA: cuantos pasos va a tener la simulacion
#	$3 writeDataS: cada cuanto se van a escrtibir los datos de la simulacion 
#	$4 paddingGird: tamaño sobresaliente de la grid para que el cutoff (rcoulomb,rvdw)sea igual
#	$5 tc-grps
#	$6 tau_t
#	$7 ref_t
#	$8 energyGrops
#	$9 nstlist
#	$10 assign velocities from Maxwell distributio seed
#	
#
printf '%-20s %-2s %-20s %-40s\n' "title"				"=" "NVT `basename ${file_conf_nvt}  `"	""						            							> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "define"				"=" "-DPOSRES"					"; position restrain the protein"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											    "; Run parameters"											>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "integrator"			"=" "md"						"; leap-frog integrator"									>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "tinit" 				"=" "0"							"; inint step"  											>> ${file_conf_nvt} #comienza la equilibracion en 0
printf '%-20s %-2s %-20s %-40s\n' "nsteps"				"=" "${step_nvt}"				"; 2 * 50000 = 100 ps"										>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "dt"		    		"=" "0.002"						"; 2 fs"													>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""												"; Output control"											>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "nstxout"				"=" "${write_data}"				"; save coordinates every 1.0 ps"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "nstvout"				"=" "${write_data}"				"; save velocities every 1.0 ps"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "nstenergy"			"=" "${write_data}"				"; save energies every 1.0 ps"								>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "nstlog"				"=" "${write_data}"				"; update log file every 1.0 ps"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											    "; Bond parameters"											>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "continuation"	    "=" "no"						"; first dynamics run"										>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "constraint_algorithm" "=" "lincs"	    			"; holonomic constraints "									>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "constraints"	        "=" "all-bonds"					"; all bonds (even heavy atom-H bonds) constrained"			>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "lincs_iter"	        "=" "1"		    				"; accuracy of LINCS"										>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "lincs_order"	        "=" "4"		    				"; also related to accuracy"								>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											    "; Neighborsearching"										>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "cutoff-scheme"   	"=" "Verlet"																				>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "ns_type"		    	"=" "grid"						"; search neighboring grid cells"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "nstlist"		    	"=" "${nstlist}"				"; 20 fs, largely irrelevant with Verlet"					>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""												        "; Electrostatics"									>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "coulombtype"	   		"=" "PME"						"; Particle Mesh Ewald for long-range electrostatics"		>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "pme_order"	    	"=" "4"							"; cubic interpolation"										>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "fourierspacing"		"=" "0.16"						"; grid spacing for FFT"									>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""												"; Temperature coupling is on"								>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "tcoupl"				"=" "Nose-Hoover"	            "; modified Berendsen thermostat"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "tc-grps"				"=" "${tc_grps}"			    "; two coupling groups - more accurate"						>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "tau_t"				"=" "${tau_t}"      			"; time constant, in ps"									>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "ref_t"				"=" "${ref_t}"    				"; reference temperature, one for each group, in K"			>> ${file_conf_nvt} #proteina de humano 310
printf '%-20s %-2s %-20s %-40s\n' "" "" ""												"; Pressure coupling is off"								>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "pcoupl"				"=" "no" 						"; no pressure coupling in _NVT"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""												"; Periodic boundary conditions"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "pbc"					"=" "xyz"		   				"; 3-D PBC"													>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" "" 												"; Dispersion correction"									>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "DispCorr"			"=" "EnerPres"					"; account for cut-off vdW scheme"							>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "" "" "" 												"; Velocity generation"										>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "gen_vel"				"=" "yes"						"; assign velocities from Maxwell distribution"				>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "gen_temp"			"=" "300"						"; temperature for Maxwell distribution"					>> ${file_conf_nvt} #proteina de humano 310
printf '%-20s %-2s %-20s %-40s\n' "gen_seed"			"=" "${seedg}"					"; generate a random seed"									>> ${file_conf_nvt}
#______ADD
printf '%-20s %-2s %-20s %-40s\n' "rcoulomb"	    	"=" "${padding_grid}"			"; short-range electrostatic cutoff (in nm)" 				>> ${file_conf_nvt} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-20s %-40s\n' "rvdw"		    	"=" "${padding_grid}"			"; short-range van der Waals cutoff (in nm)" 				>> ${file_conf_nvt} ##igual que el tamaño de &paddingGrid
#printf '%-20s %-2s %-20s %-40s\n' "energygrps" 			"=" "${energy_grps}"			""															>> ${file_conf_nvt}
printf '%-20s %-2s %-20s %-40s\n' "comm_mode" 			"=" "linear"				    ";Remove center of mass translation "					    >> ${file_conf_nvt}
