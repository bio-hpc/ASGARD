#!/bin/bash
TAG="tecniqueBD"
findCords 						#almancena las coordenadas de la proten un array
final=`expr $fin - 1` 			#le resto uno a final para que el siguiente job se corresponda
for num_execution in `seq $ini $final`;		#en la secuencia entre inicio y final
do
    if  [ $((${num_execution}%${bd_exhaustiveness})) == 0 ];then

	    extractDataLine "${coords[num_execution]}"
	    executeTecnique "bash ${path_cluster_nodes}run_software.sh  -c $CWD -q $query -nq $name_query -na ${numAminoacdo} -ch ${chain} -ej ${num_execution} "	    
	 fi
done