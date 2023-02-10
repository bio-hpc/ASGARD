#!/bin/bash

#
#	Directorio del experimento
#
folder_out_ucm=${folder_experiment}soft_out/         	#the folder for the output log files
folder_molec=${folder_experiment}molecules/         # the folder where molecules files are found
folder_energies=${folder_experiment}energies/         # the folder where energy files are found
folder_inputs=${folder_experiment}inputs/           #inputs for the analysis

arrayFolders=( $folder_experiment $folder_out_jobs $folder_templates_jobs $folder_jobs_done $folder_grid $folder_out_ucm $folder_molec $folder_energies,$folder_inputs)
