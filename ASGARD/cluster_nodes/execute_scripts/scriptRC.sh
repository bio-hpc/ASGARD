#!/usr/bin/env bash
#
#   No funciona
#
execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
	execute "ligVS=`echo $query | sed 's/.$//g'` "
	execute "ligVS=${ligVS##*/}"
	execute "salida2=${directorio}${numeroEjecucion}${option}-${software}-${name_target}-${ligVS}-${x}-${y}-${z}"


	execute "export OE_LICENSE=${pathSW}eon/openeye/bin/oe_license.txt"
	execute "${pathSW}rocs/bin/rocs -query ${CWD}$target -dbase ${CWD}$query -prefix  ${salida2}ROC 2>/dev/null"
	execute "a=`cat ${salida2}ROC_1.rpt | tail -1 |awk '{print $5}'`"
	execute "${pathSW}eon/openeye/bin/eon -dbase ${salida2}ROC_hits_1.oeb  -prefix ${salida2}EON 2>/dev/null"
	execute "b=`cat ${salida2}EON.rpt | tail -1 |awk '{print \$7}'`"
	execute "rm ${salida2}ROC*"
	execute "rm ${salida2}EON.log"
	execute "rm ${salida2}EON.status"
	execute "rm ${salida2}EON.parm"
	execute "echo $b $x $y $z $query 1 $a >	${salida2}.txt"

}
