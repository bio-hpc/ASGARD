#!/usr/bin/env bash

printf '%-20s %-2s %-10s %-40s\n' "" "" ""  "; ions.mdp - used as input into grompp to generate ions.tpr" 											> ${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "" "" ""  "; Parameters describing what to do, when to stop and what to save" 									>>${file_conf_tpr}

printf '%-20s %-2s %-10s %-40s\n' "integrator"  "=" "steep"		    "; Algorithm (steep = steepest descent minimization)"						    >>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "emtol"	    "=" "1000.0"  	    "; Stop minimization when the maximum force < 1000.0 kJ/mol/nmW"				>>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "emstep"      "=" "0.01"          "; Energy step size"															>>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "nsteps"	    "=" "10000"	  	    "; Maximum number of (minimization) steps to perform"						    >>${file_conf_tpr} #puede ser menor en los ejemplo viene con 10000
printf '%-20s %-2s %-10s %-40s\n' "" "" ""                          "; Parameters describing how to find the neighbors of each atom and how to calculate the interactions"	>>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "nstlist"	    "=" "1"		        "; Frequency to update the neighbor list and long range forces"	>>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "cutoff-scheme" "=" "Verlet"																		>>${file_conf_tpr} #Verlet
printf '%-20s %-2s %-10s %-40s\n' "ns_type"		"=" "grid"		    "; Method to determine neighbor list (simple, grid)"			>>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "coulombtype"	"=" "PME"		    "; Treatment of long range electrostatic interactions"			>>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n' "pbc"		    "=" "xyz" 		    "; Periodic Boundary Conditions (yes/no)"						>>${file_conf_tpr}
#_______________________- ADDD
printf '%-20s %-2s %-10s %-40s\n' "rcoulomb"	"=" ${padding_grid} "; short-range electrostatic cutoff (in nm)" 					>>${file_conf_tpr} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-10s %-40s\n'  "rvdw"		"=" ${padding_grid}	"; short-range van der Waals cutoff (in nm)" 					>>${file_conf_tpr} ##igual que el tamaño de &paddingGrid
printf '%-20s %-2s %-10s %-40s\n'  "gen-seed"   "=" "2015"			"; Semilla de la rpueba 	"	                		        >>${file_conf_tpr}
printf '%-20s %-2s %-10s %-40s\n'  "comm_mode" 	"=" "linear"		"; Remove center of mass translation "					        >> ${file_conf_tpr}
