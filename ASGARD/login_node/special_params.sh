#!/usr/bin/env bash
#
#	Lee las plantlas de scriptLanzador/templateXX.txt
#   Lee los posibles parametreos de la plantilla de sofware y los devuelve en un array (file)
#
read_template_params()
{
	file=""
	contadorFile=0

	if [ -f ${path_login_node}/templateParams/template${1}.txt ];then
		auxIFS=$IFS
		while read -r line
		do
   			if [[ "$line" != \#* ]] && [[ -n $line ]];then
   				contadorFile=`expr $contadorFile + 1`
				file[$contadorFile]="$line"
			fi
		done < "${path_login_node}/templateParams/template${1}.txt"
		IFS=$auxIFS
	else
		echo "No existe plantilla para ese SW"
		exit
	fi
}
