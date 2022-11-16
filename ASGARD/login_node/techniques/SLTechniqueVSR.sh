#!/usr/bin/env bash
#lanzTechnique:VSR - Virtual screening Receptors
#______________________________________________________________________________________________________________________________________________#
pritnDebug()
{
        debugB "SLTechniqueVSR.sh: find ${CWD}${query}  -name \"*${ext_query}\" |wc -l"
        debugB "SLTechniqueVSR.sh num targets: ${numFicheros}"
        debugB "SLTechniqueVSR.sh: salida=$directory$option"_"$software"_"${name_target}-$name_query-$x-$y-$z"

        debugB "SLTechniqueVSR.sh: send_jobs"
}
#________________________________________________________________________________________________________________________________
#
#       VS
#_________________________________________________________________________________________________________________________________
function check_coords(){
    #   checkquea que los nombres de los ligandos existan en la carpeta y que tenga 4 campos el fichero coords.txt
    file_coords=`cat ${CWD}${query}/coords.txt`
    for i in $(echo $file_coords | tr " " "\n")
    do
        num_fields=`echo "$i" | awk -F':' '{ print NF }'`
        file_lig=`echo $i |awk -F: '{print $1}'`
        if [ ! -f  $query${file_lig} ] || [ $num_fields -lt 4 ] ;then
            echo "ERROR: Receptor in coords.txt ${file_lig}"
            echo "ERROR: File no exits:  $query${file_lig}"
            exit
        fi
    done
}
if [ ! -f ${CWD}${query}/coords.txt ];then
    echo "No existe el fichero de coordenadas"
    exit
fi
check_coords
source ${path_login_node}lanza_job.sh
numFicheros=`find ${CWD}${query}/ -maxdepth 1 -name "*${ext_query}" |wc -l`
salida=${option}_${software}_${name_target}_${name_query}_${x}_${y}_${z}
outGrid=${folder_grid}${salida}
outDock=${folder_grid}${salida} #solo se usa con FB por que tiene sistema de fichero en el docking
pritnDebug
send_jobs








