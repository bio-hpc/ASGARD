#_________________________________________________________________________________________________
#
#    Autor:        Jorge de la Peña
#    Email:        jorge.dlpg@gmail.com
#    version:      0.1
#    Descripcion:  Ejecuta diferentes softwares de docking y suma sus scores
#__________________________________________________________________________________________________
#!/bin/bash
ejecutar()  #hace Docking en la posicion indicada con el query indicado
{

	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
	option=VS 							#Se define la option VS
	changeNames 1 AD pdbqt pdbqt
	execute "source ${pathSD}executeScripts/scriptAD.sh"
	ejecutar
	opt_aux=""							#hay que buscar una manera de poderle pasar los parametros
	changeNames 2 LF mol2 mol2
	execute "source ${pathSD}executeScripts/scriptLF.sh"	
	ejecutar
	readData
	
}
changeNames()
{
	#
	#	Genera el nombre para esa prueba LF, AD
	#
	software=$2
	execute "salida=${numeroEjecucion}${option}_${software}_${name_target}_${name_query}_${x}_${y}_${z}"
	target="${target%.*}".$3
	query="${query%.*}".$4
	outTxt=${folderOutUcm}${salida}
	out_grid=${folderGrid}${salida}
	out_aux=${folderOutUcm}${salida}
	out_molec=${folderOutUcm}${salida}
	out_energies=${folderOutUcm}${salida}
	txts[$1]=$outTxt
}
readData ()
{
	scores=0
	re="^-[0-9]+([.][0-9]+)?$"
	for index in ${!txts[*]}
	do
		if [ -f ${txts[$index]}".txt" ];then
			aux=`cat ${txts[$index]}".txt" |awk '{print $1}'`
			if  [[ $aux  =~ $re ]] ; then
				scores=`echo $scores + $aux | bc`
			fi
		fi
	done
	outTxt=${folderTxt}${salida}
	echo "$scores $x $y $z $query " > ${outTxt}.txt
}
