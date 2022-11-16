#!/bin/bash
#
#   Virtual screening en bloques
#
TAG="tecniqueVSB"
if [[ $software == "LS"* ]];then
	ini=`expr ${ini} + 1` #se suma una para LS
	executeTecnique "bash ${path_cluster_nodes}run_software.sh  -c $CWD -q $query  -nq `basename ${query%.*}_${ini}-${fin}`"

elif [ $software == "MT" ]; then
	executeTecnique "bash ${path_cluster_nodes}run_software.sh  -c $CWD -q $query -nq ${ini}-${fin} "
elif [ $software == "COSMIA" ];then
	 executeTecnique "bash ${path_cluster_nodes}run_software.sh  -c $CWD -q $query -nq ${ini}-${fin} "

else 
	source ${path_cluster_nodes}techniques/tecniqueVS.sh
fi

