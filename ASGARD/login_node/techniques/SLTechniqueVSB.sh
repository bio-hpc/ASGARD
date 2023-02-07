#!/usr/bin/env bash
#lanzTechnique:VSB - Virtual sceening Bloque
#______________________________________________________________________________________________________________________________________________#
#																																			   #
#													Funciones VS	BBDD																	   #
#		Para softwares como LS que es mas rapido tener un solo fichero con todos los queries, no como VS que se hace con varios numFicheros   #
#______________________________________________________________________________________________________________________________________________#
pritnDebug()
{
	debugB "_________________________________________Technique VSB ________________________________________"
	debugB "vistualScriningVSD.sh: find ${CWD}${query}  -name \"*${extension}\" |wc -l"
	debugB "vistualScriningVSD.sh: numFicheros: $numFicheros"
	debugB "vistualScriningVSD.sh: salida=$directorio$option"-"$software"-"${name_target}-$name_query-$x-$y-$z" 
	debugB "vistualScriningVSD.sh: ${pathSD}generate_grid.sh -s $software -x $x -y $y -z $z -d $salida -p $target -c $CWD -o $option"
	debugB "vistualScriningVSD.sh: funcionEnviarJobs"
}
#________________________________________________________________________________________________________________________________
#
#	Cuenta los ficheros del directorio y escoge una salida
#_________________________________________________________________________________________________________________________________
source ${path_login_node}lanza_job.sh
source ${path_login_node}lanza_job.sh
if [ ${software} == "LS" ];then #Ls tiene sus bloques de librerias
	numFicheros=`${path_external_sw}ligandScout/idbinfo -d ${query} 2> /dev/null |grep Database |awk '{print $3}'`
	#echo ${path_external_sw}ligandScout/idbinfo -d ${query}
elif [ ${software} == "LSC" ];then #Ls tiene sus bloques de librerias
    numFicheros=`cat  ${target} |grep '$$$$$\|@<TRIPOS>ATOM'  |wc -l`
else
	numFicheros=`find ${CWD}${query}/*  -name "*${ext_query}" |wc -l`
fi
pritnDebug
send_jobs



#salida=${directorio}${option}-${software}-${name_target}-${name_query}-${x}-${y}-${z}	


























###____________________________________________________________________________________________________________________________________
###
###	Bucle que envia los jobs por paquetes
###_________________________________________________________________________________________________________________________________
#funcionEnviarJobs()
# {
##	###bucle para enviar jobs
##	SALIR=0
##	final=""
##	while [  $SALIR == 0 ]
##	do	
		#if [ $contFin -gt $numFicheros ];then
		#	final="-ul"
		#fi
##		funcionJob
##		contIni=$contFin
		
##		contFin=`expr $contFin + $num_per_job`
##		if [ $contFin -gt `expr $numFicheros + 1` ]; then #+1 Ojo
##	  		 SALIR=1
##		fi	
##	done
##	if [ $resto -ne 0 ];then 
		#if [ $contFin -gt $numFicheros ];then
		#	final="-ul"
		#fi
##		contFin=`expr $contIni + $resto`
##		funcionJob
##	fi
	#funcionSubirExcell 	#sube al escel
##}
