#!/bin/bash
lig=""
TAG="scriptBDVS.sh"
#
#	encuenta el query para hacer docking, como son varios hay que buscarlo
#
encontrarLig()
{
	contador=1;
	for lig in $(ls ${query}*${extension})
	do
		if [ $1 -eq $contador ];then
			break
		fi
		contador=`expr $contador + 1`
	done
}
#
#	Incrementa ini
#
sumIni()
{
	ini=`expr $ini + 1`
}
#
#	Se le pasa el inicio el query y el final
#
funcionSend()
{
	extractDataLine "${coords[$1]}"												#cogemos la fila de coordenadas de su CA
	name_query=`basename ${lig}`												#nombre del query sin la ruta del directorio de ejecución
	executeTecnique "${pathSD}runCommand.sh -c ${CWD} -l ${lig} -ej ${ini} -nl ${name_query}"
}
findCords
numLigandos=`find $query -name "*${ext_query}" | wc -l` 						#numero de queries para ejecutar el BDVS
numAminoAcidoInicial=`expr $ini / $numLigandos`									#numero del aminoacido incial
numAminoAcidoFinal=`expr $fin / $numLigandos`									#numero del aminoacdio final
numLigInicial=`expr $ini % $numLigandos`										#numero del liagndo inicial
numLigInicial=`expr $numLigInicial + 1`											#se le suma 1 para que no empieze donde donde acabo el otro job
numLigFinal=`expr $fin % $numLigandos`											#numero del query final%num queries para que no se salga del rango de la siguiente tabla
#
#
#	Sependiendo del numero de jobs se van cogiendo grupos en este ejepl del 1 al 21 por ejemplo si son 3 jobs se cogen del 1 al 7 lugeo del 8 al 14 luego del 15 al 21,
#	los queries solo pueden ser 1 2 o 3 en este ejemplo
#
#	 	   L1   L2 L3			   L1 L2 L3
#		_____________		   	_____________
#Amin	0| 1  2  3	|			0| 1  2  3  |
#Amin	1| 4  5  6	|			1| 1  2  3  |
#Amin	2| 7  8  9	|			2| 1  2  3  |
#Amin	3|10 11 12	|			3| 1  2  3  |
#Amin	4|13 14 15 	|			4| 1  2  3  |
#Amin	5|16 17 18	|			5| 1  2  3  |
#Amin	6|19 20 21	|			6| 1  2  3  |
#		 |__________|			 |__________|
#

#outGrid=${directorio}grids/${option}"_"${software}"_"${name_target}_${name_query}_${x}_${y}_${z}
#echo $outGrid;




for i in `seq $numLigInicial $numLigandos`;do 										#se comienza enviando desde el incio hasta el final OJO aqui puede haber problemas si son cientos de queries
	encontrarLig $i
	funcionSend $numAminoAcidoInicial
	#echo $ini $lig "aminoacid " $numAminoAcidoInicial
	sumIni
done
diferencia=`expr $numAminoAcidoFinal - $numAminoAcidoInicial` 						#diferencias para sabercuantos reiduos hay que hacer dockin
for i in `seq 1 $diferencia`;do #se emepieza en uno por que seg teermina en el ultimo nomeuro <=
	for j in `seq 1 $numLigandos`;do
		#echo $numAminoAcidoFinal  `expr $numLigInicial + 1`
		if [ `expr $i + $numAminoAcidoInicial` -lt $numAminoAcidoFinal ];then
			encontrarLig $j
			funcionSend `expr $i + $numAminoAcidoInicial`
			#echo $ini $lig "aminoacid " `expr $i + $numAminoAcidoInicial`
			sumIni
		else
			if [ $j -le `expr $numLigFinal` ];then
				encontrarLig $j
				funcionSend `expr $i + $numAminoAcidoInicial`
				#echo $ini $lig "aminoacid " `expr $i + $numAminoAcidoInicial`
				sumIni
			fi
		fi
	done
done
