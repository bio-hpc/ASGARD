#!/usr/bin/env bash
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta

execute_script()
{
	execute "${path_external_sw}shafts/cynthia -q ${CWD}${target} -t ${CWD}${query}  -o ${out_molec} $opt_aux &>${out_aux}.ucm"
    if [ $? == 0 ];then
		execute "tail -n 1 ${out_molec}Result.list |awk '{print \$3\" '$x' '$y' '$z' '$query' \"\$4}' >${outTxt}.txt "
		execute "rm ${out_molec}Result.list"
		execute "rm ${out_aux}.ucm"
	fi
}

#_____________ Antonio _________________-
#
#       ${CWD}scriptsDocking/externalSw/shafts/cynthia -q ${CWD}${proteina} -t ${CWD}${ligando} -o ${salida} -postOpt -normalizeType tanimoto $opt_aux &>${salida}.ucm
#       if [ $? == 0 ];then
#         sed -i '1d' ${salida}Result.list
#         cd ${directorio}
#         python ${CWD}scriptsWeb/splitSF.py ${salida}Hits.mol2
#         rm ${salida}.ucm
#       fi
# ;;


