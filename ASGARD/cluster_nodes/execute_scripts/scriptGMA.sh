#!/usr/bin/env bash
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta

execute_script()
{
    execute "${path_external_sw}GMA/GMA-UNIX.exe $CWD$query $CWD$target ${out_molec}.sdf ${out_energies}.eng &>${out_aux}.ucm"
    hits=`tail -n 1 ${out_aux}.ucm |awk '{print $7}'`
    if [ "$hits" == "hits=0" ];then #si no ha encontrado hits se eliminan los datos
         rm ${out_molec}.sdf
         rm ${out_energies}.eng
    fi
}

#result=`cat ${out_aux}.ucm`
#hits=`echo ${result:$((${#resut}-1)):1}`
