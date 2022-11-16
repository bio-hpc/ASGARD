#!/usr/bin/env bash
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
execute_script()
{
    execute "score=\"`${path_external_sw}lisica/lisica -R ${CWD}${target} -T ${CWD}${query} -n ${cores} ${opt_aux} |awk '{print $1}'`\""
    if (( $(echo $score '>' ${scoreCorte} | bc -l) )); then
	   execute "echo $score > ${out_energies}.eng"
    fi
}
#_____________________________- Antonio
# #${CWD}scriptsDocking/externalSw/lisica/lisica -R $proteina -T $ligando -d 2 -w 100 -n 1 -f $salida > $salida".txt"
# head -n-2 ${salida}.txt | awk -v SALIDA=${salida} '{if ($1 ~ /[0-9]/) {system("cp "SALIDA"/*/*"$2".mol2 "SALIDA"/../VS-LI-"$2".mol2")} else {system("rm "SALIDA"/*/*"$2".mol2"); system("sed -i -e \"/"$2"$/d\" "SALIDA".txt")}}'
















#./lisica[.exe] -R <path to reference molecule> -T <path to target molecules> [parameters];
#-n <number of CPU threads> Default value: the default is to try to detect the number of CPUs or, failing that, use 1
#-d <product graph dimension> Possible input: 2 or 3; Default value: 2
#-m <maximum allowed atom spatial distance for the 3D product graph measured in angstroms> Default value: 1.0
#-s <maximum allowed shortest path size for the 2D product graph measured in the number of covalent bonds> Default value: 1
