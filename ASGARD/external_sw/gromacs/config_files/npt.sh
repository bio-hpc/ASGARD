
#!/usr/bin/env bash
printf '%-20s %-2s %-20s %-40s\n'  "title"				"=" "NPT `basename $file_conf_npt`"											                > ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "define"				"=" "-DPOSRES"				"; position restrain the protein"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""									  		"; Run parameters"										>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "integrator"			"=" "md"					"; leap-frog integrator"								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "tinit" 				"=" "${step_step}"			"; tiempo de inicio de la simulacion"			        >> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "nsteps"				"=" "${step_npt}"			"; 2 * 50000 = 100 ps"							        >> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "dt"		    		"=" "0.002"					"; 2 fs"												>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Output control"										>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "nstxout-compressed" "=" "${write_data}"			"; save coordinates every 1.0 ps"				        >> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "nstvout"			"=" "${write_data}"			"; save velocities every 1.0 ps"					    >> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "nstenergy"			"=" "${write_data}"			"; save energies every 1.0 ps"					        >> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "nstlog"				"=" "${write_data}"			"; update log file every 1.0 ps"				    	>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""										    "; Bond parameters"										>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "continuation"	      "=" "yes"					"; Restarting after NVT "								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "constraint_algorithm" "=" "lincs"	    		"; holonomic constraints "								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "constraints"	      "=" "all-bonds"			"; all bonds (even heavy atom-H bonds) constrained"		>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "lincs_iter"	          "=" "1"		    		"; accuracy of LINCS"									>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "lincs_order"	      "=" "4"		    		"; also related to accuracy"							>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Neighborsearching"									>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "cutoff-scheme"   	  "=" "Verlet"	    															    >> ${file_conf_npt} ##esquema
printf '%-20s %-2s %-20s %-40s\n'  "ns_type"		      "=" "grid"				"; search neighboring grid cells"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "nstlist"		      "=" "${nstlist}"	    	"; 20 fs, largely irrelevant with Verlet scheme"		>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Electrostatics"										>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "coulombtype"		  "=" "PME"					"; Particle Mesh Ewald for long-range electrostatics"	>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "pme_order" 		      "=" "4"		    		"; cubic interpolation"									>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "fourierspacing"		  "=" "0.16"				"; grid spacing for FFT"								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Temperature coupling is on"							>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "tcoupl"				  "=" "Nose-Hoover"	        "; modified Berendsen thermostat"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "tc-grps"			  "=" "${tc_grps}"			"; two coupling groups - more accurate"					>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "tau_t"				  "=" "${tau_t}"			"; time constant, in ps"								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "ref_t"				  "=" "${ref_t}"			"; reference temperature, one for each group, in K"		>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Pressure coupling is on"								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "pcoupl"		          "=" "Parrinello-Rahman"	"; Pressure coupling on in _NPT"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "pcoupltype"	          "=" "isotropic"           "; uniform scaling of box vectors"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "tau_p"		          "=" "2.0"		            "; time constant, in ps"								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "ref_p"		          "=" "${pressure_npt}"	    "; reference pressure, in bar"							>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "compressibility"      "=" "4.5e-5"	            "; isothermal compressibility of water, bar^-1"			>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "refcoord_scaling"     "=" "com"																			>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Periodic boundary conditions"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "pbc"				  "=" "xyz"					"; 3-D PBC"												>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Dispersion correction"								>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "DispCorr"			  "=" "EnerPres"			"; account for cut-off vdW scheme"						>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; Velocity generation"									>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n'  "gen_vel"			  "=" "no"				    "; Velocity generation is off "							>> ${file_conf_npt} # las generacion de velocidades viene del  el nvt
printf '%-20s %-2s %-20s %-40s\n'  "" "" ""											"; ADD"									        		>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n' "rcoulomb"	 		  "=" "${padding_grid}"	    "; short-range electrostatic cutoff (in nm)" 			>> ${file_conf_npt} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-20s %-40s\n' "rvdw"		    	  "=" "${padding_grid}"	    "; short-range van der Waals cutoff (in nm)" 			>> ${file_conf_npt} ##igual que el tamaño de &paddingGrid
#printf '%-20s %-2s %-20s %-40s\n' "energygrps"            "="  "${energy_grps}"												 				>> ${file_conf_npt}
printf '%-20s %-2s %-20s %-40s\n' "comm_mode" 			  "=" "linear"				"; Remove center of mass translation "					>> ${file_conf_npt}








#printf '%-20s %-2s %-20s %-40s\n' "gen-seed"             "=" "2015"					"; Semilla de la rpueba 					"			>> ${file_conf_npt}
##printf '%-20s %-2s %-20s %-40s\n'  "tc-grps"			  "=" "Protein Non-Protein"	"; two coupling groups - more accurate"					>> ${file_conf_npt}
##printf '%-20s %-2s %-20s %-40s\n'  "ref_t"			  "= ${tc_grps}	  ${tc_grps}	        ; reference temperature, one for each group, in K"		>> ${file_conf_npt}