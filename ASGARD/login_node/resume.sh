#!/bin/bash

##OJO HAY PARETES DEL INFORME QUE SE GENERAN EN OTROS FICHEROS
generate_informe()
{
	
	echo "___________________________________________Resume________________________________________________________">>${salidaResume}
	echo  "--">>${salidaResume}
	echo  "--">>${salidaResume}
	informe=${informe}"-- LigandIn:\t\t${query}\n"
	informe=${informe}"-- ReceptorIn:\t\t${target}\n"
	informe=${informe}"-- protocolProtIn:\t\t${protocolP}\n"
	informe=${informe}"-- protocolLigIn:\t${protocolL}\n"
	informe=${informe}"-- X:\t\t\t$x\n"
	informe=${informe}"-- Y:\t\t\t$y\n"
	informe=${informe}"-- Z:\t\t\t$z\n"
	informe=${informe}"-- Grid:\t\t$grid\n"
	informe=${informe}"-- SizeGridX:\t\t$gridSizeX\n"
	informe=${informe}"-- SizeGridY:\t\t$gridSizeY\n"
	informe=${informe}"-- SizeGridZ:\t\t$gridSizeZ\n"
	title=`cat ${path_login_node}templateParams/template${software}.txt |grep "#LanzNombre" | cut -d : -f 2`
	informe=${informe}"-- Software:\t\t$title\n"
	#sinforme=${informe}"-- Software:\t\t$software\n"
	informe=${informe}"-- Technique:\t\t$option\n"
	informe=${informe}"-- Cores:\t\t$cores\n"
	informe=${informe}"-- Mem:\t\t\t$mem\n"
	informe=${informe}"-- NomJob:\t\tname_job\n"
	informe=${informe}"-- QueueManager:\t$gestorColas\n"
	informe=${informe}"-- Directory:\t\t$folder_experiment\n"
	informe=${informe}"-- TotalRuns:\t\t$numFicheros\n"
	informe=${informe}"-- run x Jobs:\t\t$nRuns\n"
	informe=${informe}"-- Jobs:\t\t$num_per_job\n"
	informe=${informe}"-- Time Job:\t\t$time_job\n"
	informe=${informe}"-- Time Dock:\t\t$time_experiment\n"
	cl=`hostname`
	informe=${informe}"-- Cluster:\t\t$cl\n"
	informe=${informe}"--\n--\n"
	revision=`git rev-list --reverse HEAD 2>/dev/null |tail -1 `
	codVersion=`git rev-list --reverse HEAD 2>/dev/null |tail -1`
	numVersion=`git rev-list --reverse HEAD 2>/dev/null |wc -l`
	rama=`git branch -avv 2>/dev/null|grep "\*" |awk '{print $2}'`
	informe="${informe}-- numVesrion:\t $numVersion codVersion:\t\t$codVersion  \t branch: $rama\n" #fecha para el informe
	informe="${informe}-- date:\t\t$fecha (Y-m-d)\n" #fecha para el informe
	informe=${informe}"-- Command:\t\t$allComand $optAux"
	echo -e "${informe}" >>${salidaResume}
	echo  "--">>${salidaResume}
	echo  "--">>${salidaResume}
	if [ "$grid" == "Y" ] ;then
		echo  "___________________________________________GenerateGrid Script________________________________________________________" >>${salidaResume}
		cat ${path_cluster_nodes}grids/generateGrid${software}.sh >>${salidaResume}
		echo  "--">>${salidaResume}
		echo  "--">>${salidaResume}
	fi
}
getSteemSimulationGR()
{
	
	OLDIFS=$IFS

	IFS='-' read -r -a array <<< "$1"
	inicio=0
	for element in "${array[@]}"
	do
		aux=`echo $element | cut -f1 -d' '`
		param=`echo $element | cut -f2 -d' '`
		case "$aux" in
			stepSimulacionDM) 
				#optAux=`echo "$optAux" |sed -e "s,-stepSimulacionDM $param,,g"`
				#echo $param
				sSimu=$param;; 
		esac
	done
	IFS=$OLDIFS;
	echo $sSimu
}

#
#	Para cualquier software que no sea gromacs Next se genera un informe normal
#
#	Opciones:
#	$lanzCreateResume=Y --> se crea un resumen nuevo
#	$lanzCreateResume=C --> se Continua un resumen y se actualizan datos
#	$lanzCreateResume=N --> no se toca el reumen
if [ "$lanzCreateResume" == "Y" ];then #se crea un resumen nuevo
	salidaResume=${folder_experiment}Resume_${name_job}_${name_target}_${name_query}.txt
	generate_informe
elif [ "$lanzCreateResume" == "C" ];then #se continua un resumen
	file_resume=`ls -hrt ${folder_experiment}/Resume_*txt |tail -1`				#resumen anterior (Ultimo resumen)
    num_last_command=`cat ${file_resume} |grep -n Command |tail -1 |awk -F: '{print $1}'`
    num_last_command=`expr $num_last_command + 1`   #linea siguiente el ultimo command
    commnad_n=`cat ${file_resume} | grep Command_ |tail -1 |awk -F: '{print $1}'`

    if [[ "" == "$commnad_n" ]];then
        command_number="Command_1"
    else
        n_command=`echo ${commnad_n} |awk -F-- '{print $2}'`
        n_command=`echo ${n_command} | awk -F_ '{print $2}'`
        n_command=`expr $n_command + 1`
        command_number=Command_${n_command}
    fi
    sed  -i "${num_last_command} i --\n-- ${command_number}:\t\t$allComand $optAux" ${file_resume}
fi
