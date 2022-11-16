#!/bin/bash
# hay que pasarle como parametro un directorio con pmz y una libreria
if [ $# -ne 2 ]; then
	echo ""
	echo "Debe introducir un directorio con  donde se encuentren las pmz y una libreriea"
	echo "Ejemplo: "
	echo "sh ShuttleMol/extra_shuttlemol/launcher_ls.sh targets/pmz/ queries/Drugbank.ldb"
	echo ""
	exit

fi
JOBS_PER_EXECUTION=15
#library="libraries/DBALLv503.ldb/"
#library="libraries/FDB.ldb/"
library=$2
lib="$(basename $library)"
lib="${lib%.*}"
for i in $(ls ${1});do
		filename=${i%.*}
		extension="${i##*.}"
		if [ $extension == "pmz" ];then	
			for j in $(seq 0 5);do		
			#echo ${1}${i}  
				./sm.sh -t ${1}/${i} -q ${library} -s LS -o VSB -j ${JOBS_PER_EXECUTION}  -prp LS -prl LS -tj 24:00:00  -d LS_${filename}_${lib}_a_${j}   -a ${j} -TT 50m  -sf relative -el n
			done
		fi
done
