#!/bin/bash

#
#	Directorio del experimento
#
#folder_error_jobs=${folder_experiment}out/            #directorio donde se guardan los errores de los jobs
#folder_out_jobs=${folder_experiment}jobs_out/              #directorio donde se guardan la salida de los jobs
#folder_templates_jobs=${folder_experiment}jobs/  	#directorio donde se guardan los templates que se lanzan de los jobs
#folder_jobs_done=${folder_experiment}jobs_done/      #directorio donde se guardan los templates que se lanzan de los jobs
#folder_grid=${folder_experiment}grids/         		#directorio donde se guardan las grid en caso neceario
folder_out_ucm=${folder_experiment}soft_out/         	#directorio para salidas auxiliares
folder_molec=${folder_experiment}molecules/         #directorio donde se guardan las moleculas de la prueba
folder_energies=${folder_experiment}energies/         #directorio donde se guardan las energias de los dockings
folder_inputs=${folder_experiment}inputs/           #inputs de la dinámica (archivo .itp)

arrayFolders=( $folder_experiment $folder_out_jobs $folder_templates_jobs $folder_jobs_done $folder_grid $folder_out_ucm $folder_molec $folder_energies,$folder_inputs)
