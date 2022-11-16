#!/usr/bin/env bash
execute_script()
{
	unset DISPLAY
	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
    execute "${path_external_sw}ligandScout/idbgen -i ${CWD}${target} -M ${mem} -C ${cores} -F ${ini} -L ${fin} ${opt_aux}  -o ${out_molec}.ldb -e ${out_aux}.err &> ${out_aux}.log"
    #execute "${path_external_sw}ligandScout/idbgen -i ${CWD}${target} -d ${CWD}${query} -M ${mem} -C ${cores} -F ${ini} -L ${fin} ${opt_aux}  -o ${out_molec}.sdf"

}


