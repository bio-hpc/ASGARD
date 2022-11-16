#!/usr/bin/env bash

pathSL="${2}ShuttleMol/login_node/"
source ${2}ShuttleMol/login_node/debug.sh 							#Para que funcione el debug
source ${2}ShuttleMol/login_node/colors.sh			  				#colores del script
source ${2}ShuttleMol/login_node/parameters.sh 		  				#parametros
source ${CWD}ShuttleMol/login_node/folders_experiment.sh 			#variables folders de la prueba
source ${CWD}ShuttleMol/login_node/read_file_conf.sh
read_file ${CWD}/ShuttleMol/