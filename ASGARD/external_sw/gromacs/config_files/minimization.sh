#!/usr/bin/env bash
#
#	Configuracion para minimizar la enrgia
#	Se le pasa
#	$1 ruta deonde se guardara el fichero
#	$2 stepMinimizacion: cuantos pasos va a tener la simulacion
#	$3 paddingGird: tamaño sobresaliente de la grid para que el cutoff (rcoulomb,rvdw)sea igual
#	$4 ligNom: Nombre del ligando para el indice
#

printf '%-20s %-2s %-10s %-40s\n' "" "" ""                                  "; minim.mdp - used as input into grompp to generate em.tpr"			    > ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "integrator"	"=" "steep"		            "; Algorithm (steep = steepest descent minimization)"						>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "emtol"		"=" "1000.0"  	            "; Stop minimization when the maximum force < 1000.0 kJ/mol/nm"				>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "emstep"      "=" "0.01"                  "; Energy step size"														>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "nsteps"		"=" "${step_min}"	  	    "; Maximum number of (minimization) steps to perform"						>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "" "" ""                                  "; Parameters describing how to find the neighbors of each atom and how to calculate the interactions" >> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "nstlist"		    "=" "1"		            "; Frequency to update the neighbor list and long range forces"			    >> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "cutoff-scheme"   "=" "Verlet"																			        	>> ${file_conf_min}
#echo "cutoff-scheme   = group"																				>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "ns_type"		    "=" "grid"		        "; Method to determine neighbor list (simple, grid)"						>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "coulombtype"	    "=" "PME"		        "; Treatment of long range electrostatic interactions"					    >> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "pbc"		        "=" "xyz" 		        "; Periodic Boundary Conditions (yes/no)"								    >> ${file_conf_min}
#____________________ ADDD
printf '%-20s %-2s %-10s %-40s\n' "rcoulomb"	    "=" "${padding_grid}"	"; short-range electrostatic cutoff (in nm)" 								>> ${file_conf_min} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-10s %-40s\n' "rvdw"		    "=" "${padding_grid}"	"; short-range van der Waals cutoff (in nm)" 								>> ${file_conf_min} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-10s %-40s\n' "gen-seed"      "=" "2015"		        "; Seed "																	>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "comm_mode" 			"=" "linear"	    "; Remove center of mass translation "				                    	>> ${file_conf_min}
printf '%-20s %-2s %-10s %-40s\n' "reseed" 			"=" "${seedg}"	    		"; "				      									              	>> ${file_conf_min}
	



#
#	Add par iplicit solvent
#
##echo " implicit_solvent    =  GBSA "					>> ${file_conf_min}
##echo " gb_algorithm        =  Still ; HCT ; OBC"		>> ${file_conf_min}
##echo " nstgbradii          =  1"						>> ${file_conf_min}
##echo " rgbradii            =  0   ; [nm] Cut-off for the calculation of the"	>> ${file_conf_min}
##echo " Born"											>> ${file_conf_min}
##echo " radii. Currently must be equal to rlist"			>> ${file_conf_min}
##echo " gb_epsilon_solvent  =  80    ; Dielectric constant for the implicit"		>> ${file_conf_min}
##echo " solvent"								>> ${file_conf_min}
##echo " ; gb_saltconc       =  0     ; Salt concentration for implicit solvent"	>> ${file_conf_min}
##echo " models, currently not used"			>> ${file_conf_min}
##echo " sa_algorithm        =  Ace-approximation"	>> ${file_conf_min}
##echo " sa_surface_tension  = -1"	>> ${file_conf_min}
