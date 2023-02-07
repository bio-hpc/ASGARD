#lanzTechnique:BDVS - Blind Docking and Virtual screening (Obsoleta)
#___________________________________________________________________________________________________________________________________________________#
#																																					#
#											Funciones BDVS:																							#
#		Hace un BD con los diferentes queries de una carpeta																						#
#___________________________________________________________________________________________________________________________________________________#
pritnDebug()
{
        
        debugB "SLtecnhniqueBDVS: numLigandos: $numLigandos"
        debugB "SLtecnhniqueBDVS: Atoms Prot: ${atomsProtCA}"
        debugB "numTest: $numFicheros"
}
funcionBDVS()
{
	#funcionFindCoordsCas
	if [ "$software" == "AD" ];then
		extensionLig=".pdbqt" 		#OJO DEPENDE DEL PROGRAMA
	elif [ "$software" == "LF" ];then
		extensionLig="*.mol2"	
	else
		echo "Solo se puede hacer BDVS con LF o AD"
		ayuda
	fi	
	numLigandos=`ls ${query}*${extensionLig} |wc -l`		#numero de queries
	atomsProtCA=`python ${path_extra_shuttlemol}used_by_shuttlemol/standar_file_coords.py $target |grep -v "##" |grep CA |wc -l`	#numero de CAs (puntos donde se haran docking)
	numFicheros=`expr $numLigandos \* $atomsProtCA`        #numero total de jobs
	if [ $num_per_job -ge $atomsProtCA ];then
		echo "Script Lanzador -j no puede ser mayor o igual que el numero de carbono alfas de la target"
		ayuda
	fi
	#x="" 			#es necesario por una modificacion y poderlo hacer compatible con VS
	#y="" 			#es necesario por una modificacion
	#z="" 			#es necesario por una modificacion
	#numAminoacdo="" 	#es necesario por una modificacion
	pritnDebug
	funcionEnviarJobs

}
#_________________________________________________________________________________________________________
#
#	Ejecuta EL BDVS
#________________________________________________________________________________________________________
source ${pathSL}lanzaJob.sh
send_jobs
