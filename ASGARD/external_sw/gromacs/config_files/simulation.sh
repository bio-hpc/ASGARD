#!/usr/bin/env bash
printf '%-20s %-2s %-20s %-40s\n' "title"				"=" "MD `basename  ${file_conf_md}` " 										> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Run parameters"										>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "integrator"			"=" "md"					"; leap-frog integrator" 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "tinit" 				"=" "${int_md}"				""														>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "nsteps"				"=" "${step_md}"			"; 2 * 5000000 = 10000 ps (10 ns)" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "dt"		    		"=" "0.002"					"; 2 fs" 												>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""                                          "; Output control" 										>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "nstxout-compressed"  "=" "${write_data}"			"; save coordinates every 10.0 ps" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "nstvout"		        "=" "${write_data}"			"; save velocities every 10.0 ps" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "nstenergy"	        "=" "${write_data}"			"; save energies every 10.0 ps" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "nstlog"		        "=" "${write_data}"			"; update log file every 10.0 ps" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""                                			"; nstxout-compressed replaces nstxtcout" 				>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "compressed-x-grps"   "=" "System"    			"; replaces xtc-grps" 									>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Bond parameters"										>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "continuation"	     "=" "yes"					"; Restarting after NPT " 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "constraint_algorithm" "=" "lincs"	   			"; holonomic constraints "								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "constraints"	        "=" "all-bonds"				"; all bonds (even heavy atom-H bonds) constrained" 	>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "lincs_iter"	        "=" "1"		   				"; accuracy of LINCS"									>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "lincs_order"	        "=" "4"		    			"; also related to accuracy" 							>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Neighborsearching"									>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "cutoff-scheme"   	"=" "Verlet" 																		>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "ns_type"		    	"=" "grid"					"; search neighboring grid cells" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "nstlist"		    	"=" "${nstlist}"	    	"; 20 fs, largely irrelevant with Verlet scheme" 		>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "rcoulomb"	    	"=" "${padding_grid}"		"; short-range electrostatic cutoff (in nm)" 			>> ${file_conf_md} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-20s %-40s\n' "rvdw"		    	"=" "${padding_grid}"		"; short-range van der Waals cutoff (in nm)" 			>> ${file_conf_md} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-20s %-40s\n' "" "" "" 											"; Electrostatics"										>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "coulombtype"	    	"=" "PME"					"; Particle Mesh Ewald for long-range electrostatics" 	>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "pme_order"	    	"=" "4"		    			"; cubic interpolation" 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "fourierspacing"		"=" "0.16"					"; grid spacing for FFT" 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""										    "; Temperature coupling is on" 							>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "tcoupl"				"=" "Nose-Hoover"           "; modified Berendsen thermostat" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "tc-grps"				"=" "${tc_grps}"			"; two coupling groups - more accurate"				    >> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "tau_t"				"=" "${tau_t}"          	"; time constant, in ps"								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "ref_t"				"=" "${ref_t}"    			"; reference temperature, one for each group, in K"	    >> ${file_conf_md} #proteina de humano 310
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Pressure coupling is on" 							>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "pcoupl"		        "=" "Parrinello-Rahman"	    "; Pressure coupling on in NPT" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "pcoupltype"	        "=" "isotropic"	            "; uniform scaling of box vectors"						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "tau_p"		        "=" "2.0"		            "; time constant, in ps" 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "ref_p"		        "=" "1.0"		            "; reference pressure, in bar" 							>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "compressibility"     "=" "4.5e-5"	            "; isothermal compressibility of water, bar^-1" 		>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Periodic boundary conditions"						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "pbc"					"=" "xyz"					"; 3-D PBC" 											>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Dispersion correction" 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "DispCorr"			"=" "EnerPres"				"; account for cut-off vdW scheme" 						>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Velocity generation" 								>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "gen_vel" 			"=" "no"					"; Velocity generation is off "							>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "" "" ""											"; Extras"							 					>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "comm_mode" 			"=" "linear"				"; Remove center of mass translation "					>> ${file_conf_md}
printf '%-20s %-2s %-20s %-40s\n' "comm-grps"			"=" "${comm_grps}"			"; for center of mass motion removal, default is the whole system "	>> ${file_conf_md}
#gen_seed = -1 <-- probar
#
#	Add  iplicit solvent
#
##printf '%-20s %-2s %-20s %-40s\n' " implicit_solvent    =  GBSA "					>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " gb_algorithm        =  Still ; HCT ; OBC"		>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " nstgbradii          =  1"						>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " rgbradii            =  0   ; [nm] Cut-off for the calculation of the"	>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " Born"											>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " radii. Currently must be equal to rlist"			>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " gb_epsilon_solvent  =  80    ; Dielectric constant for the implicit"		>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " solvent"								>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " ; gb_saltconc       =  0     ; Salt concentration for implicit solvent"	>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " models, currently not used"			>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " sa_algorithm        =  Ace-approximation"	>> ${file_conf_md}
##printf '%-20s %-2s %-20s %-40s\n' " sa_surface_tension  = -1"	>> ${file_conf_md}
