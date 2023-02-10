#!/usr/bin/env bash
#
#	Description: 	Obtain where the path is and create a folder (if it does not exist)
#				
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
		if [ ! -d ${dir_name} ]; then
			mkdir ${dir_name}
		fi
	done
}
#__________________________________________________________________________________________________________________________________________
#
#	rename directory
#__________________________________________________________________________________________________________________________________________
search_directory()
{

	folder_experiment=${PWD}/${option}_${software}_${name_target}_${name_query}
	if [ "$x" == "0" ] && [ "$y" == "0" ] && [ "$z" == "0" ];then 
		if [ "$GPU" != "N/A" ];then	
			folder_experiment=${folder_experiment}_GPU 		
		else
			folder_experiment=${PWD}/${option}_${software}_${name_target}_${name_query}
		fi
	else 
		folder_experiment=${folder_experiment}_${x}_${y}_${z} 	
	fi 

	folder_experiment=${folder_experiment}_${fecha}/
	
}
existeDirectory()
{
	if [ -d $folder_experiment ];then 
		while [ "$input" != "Y" ] && [ "$input" != "y" ] && [ "$input" != "N" ] && [ "$input" != "n" ] && [ "$input" != "zz" ] ; do
			echo "CreateFile: ${folder_experiment} folder already exist. Do you want to delete it?"
			echo "(Y/y) Delete folder"
			echo "(N/n) Exit"
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

if [ "$folder_experiment" != "N/A" ];then 

	folder_experiment=${PWD}/$folder_experiment
	create_dirs
else 

	search_directory
	existeDirectory
fi		
debugB "_________________________________________Directories________________________________________"
debugB ""
for dir_name in "${arrayFolders[@]}";do
		debugB "create_dirs.sh:  $dir_name"
done



debugB ""
