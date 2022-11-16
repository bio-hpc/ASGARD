#!/usr/bin/env bash
#lanzTechnique:BD - Blind Docking
#___________________________________________________________________________________________________________________________________________________#
#																																					#
#											Funciones BD																							#
#___________________________________________________________________________________________________________________________________________________#
#______________________________________________________________________________________________________
#
#	Convierte el query en caso de que fuera diferente a pdbqt, solo pdb o mol2 , no funciona ahoa
#_____________________________________________________________________________________________________
funcionConvertLignad()
{
	extensionLigando=$(basename "$query" | cut -d. -f2) #extenion del query
	if [ $extensionLigando == ".mol2" ] || [ $extensionLigando == ".pdb" ] ;then
		fileTmp=`expr ${query%$extensionLigando}`
		#externalSw/mgltools_x86_64Linux2_latest/bin/pythonsh externalSw/mgltools_x86_64Linux2_latest/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py
		${pathSW}mgltools_x86_64Linux2_latest/bin/pythonsh ${pathSW}mgltools_x86_64Linux2_latest/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ${CWD}$fileTmp$extensionLigando -o ${CWD}$fileTmp".pdbqt" -A 'hydrogens' -U \'\'
		query=$fileTmp".pdbqt"
	fi
}
#_______________________________________________________________________________________________________
#
#		Lee un el fichero target.tmp y lanza 1 job con esas coordenadas por cada fila
#________________________________________________________________________________________
funcionBlindDocking()
{
	#funcionConvertLignad
	numFicheros=`python ${path_extra_shuttlemol}used_by_shuttlemol/standar_file_coords.py ${CWD}${target} |grep -v "##" |grep ${bd_atom_default} |wc -l`
	x=0 			#es necesario por una modificacion y poderlo hacer compatible con VS
	y=0 			#es necesario por una modificacion
	z=0 			#es necesario por una modificacion
	numAminoacdo=0 	#es necesario por una modificacion
	send_jobs

}
source ${path_login_node}lanza_job.sh

funcionBlindDocking

















#___________________________________________________________________________________________________
#
#		Genera un fichero nombreProteina.tmp con todos los CAS de la misma
#___________________________________________________________________________________________________
#funcionFindCoordsCas()
# {
#	coordsAux=`python ${pathSL}standarFileCoords.py $target |grep -v "##" |grep CA `
#
#
# }
#_____________________________________________________________________________________________________________________
#
#			Genera la flexibilidad y la guarda en un fichero targetF.txt para que luego hacer docking con vina (No se usa)
#______________________________________________________________________________________________________________________
#funcionBDFlex()
# {
#
#	ficheroFlex=${target%$extension}F.txt
#	if [ -f $ficheroFlex ];then
#		echo "existe el fichero de flexibilidad, Si desa generar flexibilidad nueva borrelo"
#		#rm $ficheroFlex
#	elif [ -n "$dinamyc" ];then
#		if [ -d $ficheroFlex ];then
#				cat  $dirDinamyc*.txt |grep -v "#" |sort -k6 -n >$target.aux
#				while read linea
#				do
#					chain=`echo $linea | awk '{print $9}'`
#					numAmino=`echo $linea | awk '{print $8}'`
#					x=`echo $linea | awk '{print $2}'`
#					y=`echo $linea | awk '{print $3}'`
#					z=`echo $linea | awk '{print $4}'`
#				  	cadenaFlex=`ssh scriptsLanzador python scriptPymolFlexCoords.py targets/$name_target.pdbqt $x $y $z $chain |sort -n |grep ATOM |uniq -f 1 |awk ' {print "'$name_target':" $3":" $2$1}'|tr '\n' ',' |sed 's/.$//g'`
#  					echo $numAmino $chain $a $x $y $z >> $ficheroFlex
#			done <  $target.aux
#		else
#				echo "Lanzador: No existe el directorio con los txt del BD anterior"
#				ayuda
#		fi
#	else
#		scp $target scriptsLanzador:targets/
#		for line in $(cat $target".Tmp"  ); # | head -1 quitar el head-3
#		do
#			x=`echo $line | cut -d\: -f1`
#			y=`echo $line | cut -d\: -f2`
#			z=`echo $line | cut -d\: -f3`
#			numAminoacdo=`echo $line | cut -d\: -f4`
#			chain=`echo $line | cut -d\: -f5`
#			cadenaFlex=`ssh scriptsLanzador python scriptPymolFlex.py targets/$name_target.pdbqt $numAminoacdo $chain|sort -n |grep ATOM |uniq -f 1 |awk ' {print "'$name_target':" $3":" $2$1}'|tr '\n' ',' |sed 's/.$//g'`
#			echo "$numAminoacdo $chain $cadenaFlex  $x $y $z >> $ficheroFlex"
#			echo "$numAminoacdo $chain $cadenaFlex $x $y $z" >> $ficheroFlex
#		done
#		ssh scriptsLanzador rm targets/$name_target.pdbqt
#	fi
#}
#____________________________________________________________________________________________
#
#	Funciones similaridad (No se usa)
#__________________________________________________________________________________________
#funcionSimilaridad()
# {
#
#	datosVS
#	if [ "$software" == "LS" ];then
#	numFicheros=`${pathSD}}ligandscout3/idbinfo -d $query |grep Database |awk '{print $3}'`	######OJO AÑADIDO -d por LS
#	div=`expr $numFicheros \/ $num_per_job`
#	resto=`expr $numFicheros \% $num_per_job`
#	contIni=1
#	contFin=`expr $contIni \+ $num_per_job`
#	extension="ldb"
#	SALIR=0
#	final=""
#	while [  $SALIR == 0 ]
#	do
#
#		funcionJob
#		contIni=$contFin
#		contFin=`expr $contFin + $num_per_job`
#		if [ $contFin -gt `expr $numFicheros + 1` ]; then #+1 Ojo
#	  		 SALIR=1
#		fi
#	done
#	if [ $resto -ne 0 ];then
#		if [ $contFin -gt $numFicheros ];then
#			final="-ul"
#		fi
#
#		contFin=`expr $contIni + $resto`
#		funcionJob
#	fi
#	fi
# }


