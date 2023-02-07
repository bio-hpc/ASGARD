#!/usr/bin/env bash
#
#	Description: 	Obtien la ruta donde se encuentra y crea un directorio (si no existiera)
#					con la nomenclatura [option]_[software]_[nombreProteina]_[query]_[|X_Y_X]_fecha
#  					Si se le indica con -d a ShuttleMol el directorio, este no se creara
#_________________________________________________________________________________________________________________________________________
asigVar()
{
	source ${path_login_node}/folders_experiment.sh
}

create_dirs()
{
	asigVar
	for dir_name in "${arrayFolders[@]}";
	do
		# Crear el directorio si no existe 
		if [ ! -d ${dir_name} ]; then
			mkdir ${dir_name}
		fi
	done
}
#__________________________________________________________________________________________________________________________________________
#
#	pone nombre al directorio ${option}_${software}_${name_target}_${nomquery}...
#__________________________________________________________________________________________________________________________________________
search_directory()
{

	folder_experiment=${PWD}/${option}_${software}_${name_target}_${name_query}
	if [ "$x" == "0" ] && [ "$y" == "0" ] && [ "$z" == "0" ];then #Si estan a 0 significa que son BD o similaridad
		if [ "$GPU" != "N/A" ];then	
			folder_experiment=${folder_experiment}_GPU 		##ruta directori
		else
			folder_experiment=${PWD}/${option}_${software}_${name_target}_${name_query}
		fi
	else #si tienen coordenadas significa 
		folder_experiment=${folder_experiment}_${x}_${y}_${z} 	##ruta para el directorio  				#
	fi 

	folder_experiment=${folder_experiment}_${fecha}/
	
}
existeDirectory()
{
	if [ -d $folder_experiment ];then #si existe el directorio se pregunta por borralro
		while [ "$input" != "Y" ] && [ "$input" != "y" ] && [ "$input" != "N" ] && [ "$input" != "n" ] && [ "$input" != "zz" ] ; do
			echo "CreateFile: El directorio ${folder_experiment} ya existe, Para continuar debera eliminar este directorio, ¿Desea eliminarlo?"
			echo "(Y/y) borrar carpeta"
			echo "(N/n) Salir"
			read  input
		done
		if [ "$input" == "Y" ] || [ "$input" == "y" ];then
			rm -r $folder_experiment
			create_dirs
		elif [ "$input" == "n" ] || [ "$input" == "N" ];then
			exit
		elif [ "$input" == "zz" ] || [ "$input" == "ZZ" ];then
			asigVar
		fi
	else #si no se crea
		create_dirs
	fi
}
#
## 	Probar solo para BDVS que sino luego al crear el clusterizado es una historia
##
#if [ `echo  $target |grep _ |wc -l` != "0" ] && [ "$option" == "BDVS" ];then
#	echo "ERROR El fichero de target no puede tener guiones bajos _"
#fi
#if [ `echo  $query |grep _ |wc -l` != "0" ] && [ "$option" == "BDVS" ];then
#	echo "ERROR El fichero de query no puede tener guiones bajos _"
#fi

#
#	Apaño para GRN se le debe indicar el directorio de la prueba antigua si o si
#
#if [ "$software" == "GRN" ];then
#	if [ "$folder_experiment" == "N/A" ];then
#		echo "Debe especificar el directorio de la prueba con el parametro -d para poder continuar"
#		exit
#	fi
#fi

#_____________________________________________________________________________________________________________________________
#
#	Si no existe el directorio: se crea
#	Si existe se le pregunta, si le da a n ShuttleMol se detiene si le da a si lo borra y lo crea si escribe zz se salta la restriciion
#	Si a intoriducio la option -d esque se quiere un directorio personalizado y se crea
#
#______________________________________________________________________________________________________________________________
if [ "$folder_experiment" != "N/A" ];then #si la cadena no esta vacia ha usado la option -d

	folder_experiment=${PWD}/$folder_experiment
	create_dirs
else #si se busca un nombre 

	search_directory
	existeDirectory
fi		
debugB "_________________________________________Directories________________________________________"
debugB ""
for dir_name in "${arrayFolders[@]}";do
		debugB "create_dirs.sh:  $dir_name"
done



debugB ""
