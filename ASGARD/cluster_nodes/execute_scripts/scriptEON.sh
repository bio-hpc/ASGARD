#!/usr/bin/env bash
#
#   Se hace primero ROC para "alineacion creo"
#   Seguidamente eon y saca un .eng en eneergies
#   NO ESTA PROBADO EN SHUTTLEMOL
#
execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)`

	execute "export OE_LICENSE=${path_external_sw}licenses/license${software}.txt"
	execute "${path_external_sw}rocs/bin/rocs -query ${CWD}$target -dbase ${CWD}$query -prefix  ${out_aux}ROC  2>${out_aux}.ucm"
	if [ $error == 0 ];then #por si falla que para
		execute "a=`cat ${out_aux}ROC_1.rpt | tail -1 |awk '{print $5}'`"
		execute "${path_external_sw}eon/openeye/bin/eon -dbase ${out_aux}ROC_hits_1.oeb  -prefix ${out_aux}EON $opt_aux 2>${out_aux}.ucm"
		if [ $error == 0 ];then #por si falla que para
			execute "score=`cat ${out_aux}EON.rpt | tail -1 |awk '{print \$7}'` "
			if [ -n $b ];then #por si el score es 0 que no se guarde la peruan

				execute "rm ${out_aux}ROC*"
				execute "rm ${out_aux}EON.log"
				execute "rm ${out_aux}EON.status"
				execute "rm ${out_aux}EON.parm"
				if [ $check != "N/A" ];then
					execute "echo $score $x $y $z $query 1 $a >	${out_energies}.eng"
				else
					 if  [[ $score  =~ $re ]] ; then
						execute "echo $score $x $y $z $query 1 $a >	${out_energies}.eng"

					else #si el score no es un numero quiere decir que no ha superado el score corte
						execute "rm ${out_aux}*"
						execute "rm ${out_aux}*"
					fi

				fi
			fi
		fi
	fi

}
