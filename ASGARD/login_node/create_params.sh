#!/usr/bin/env bash
#
#	Comprueba si existe el parametro en al plantilla y devueleve error o no
#

source ${path_login_node}/special_params.sh

read_template_params $software

for i in `seq 1 $contadorFile`;do
	IFS='::' read -ra ADDR <<< "${file[$i]}"


	#	Si es diferente a guien significa que el parametro hay que introducirlo por defecto
	#
	if [ ${ADDR[0]} != "-" ];then
		if [[ ${ADDR[4]} !=  *$allComand* ]];then
			if [[ "$optAdicionals" != *${ADDR[4]}* ]];then #si el parametro se encuentra en las optiones no se incluye el por defecto
				if [ "${ADDR[2]}" == "Y" ];then
					optAdicionals="$optAdicionals ${ADDR[4]} ${ADDR[0]} "

				else
					optAdicionals="$optAdicionals ${ADDR[4]} "
				fi
			fi
		fi
	fi
done
echo "$optAdicionals" > ${folder_templates_jobs}parameter_aux.txt  #se guardan las optiones para esa prueba

debugB "_________________________________________Extra Params________________________________________"
debugB ""
debugB "createParams.sh  $optAdicionals"
debugB ""
