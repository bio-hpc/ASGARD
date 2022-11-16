#!/usr/bin/env bash
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
execute_script()
{
 	execute "${pathSW}RF/rf-score  ${path_external_sw}RF/pdbbind-2014-refined.rf $CWD$target $CWD$query > ${out_energies}.eng 2>> ${out_aux}.ucm"
}
