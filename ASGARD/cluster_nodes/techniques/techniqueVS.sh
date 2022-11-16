#!/bin/bash

TAG="tecniqueVS"
ctr_queries=0; 									#contador para saber el rango de queries que hayq ue hacer dockin
query_aux=`basename ${query}`
for fil in ${CWD}${query}*"$ext_query"; do 				#leo todos los ficheros 1 a 1 cuando coincida con el inicio y sea menor que el final separo el lgando y llamo al comando
	if [ $ctr_queries -ge $ini ] && [ $ctr_queries -lt $fin ];then
	    name_query=${fil##$CWD}
		executeTecnique "bash ${path_cluster_nodes}run_software.sh -c ${CWD} -q ${name_query} -nq ${query_aux}" #cambiamos usammos como query el nuevo fichero y como name_query el directorio
	fi
	let ctr_queries=ctr_queries+1
done 




























