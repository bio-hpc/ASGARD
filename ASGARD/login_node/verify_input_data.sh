#!/usr/bin/env bash
#
#   Description: Validate ASGARD parameters and options
# ______________________________________________________________________________________________________________________

#
#
readParam()
{

	if [ "$1" == "N/A" ] || [  -z "$1" ] ;then
		echo `cat ${path_login_node}templateParams/template${software}.txt |grep "#${2}" | cut -d : -f 2`
	else
		echo "$1"
	fi
}
readParams()
{
	
	if [ -f "${path_login_node}/templateParams/template${software}.txt" ];then
		extensionesProt=`readParam "$extensionesProt" "LanzExtensionProtQ"`
		extensionesLig=`readParam "$extensionesLig" "LanzExtensionligB"` 
		cores=`readParam "$cores" "LanzCores"`
		time_experiment=`readParam "$time_experiment" "LanzTimeExecution"`
		x=`readParam "$x" "LanzCoordX"`
		y=`readParam "$y" "LanzCoordY"`
		z=`readParam "$z" "LanzCoordZ"`
		grid=`readParam "$grid" "LanzGrid"`
		mem=`readParam "$mem" "LanzMem"`
		optionsLanz=`readParam "$optionsLanz" "lanzOptions"`
		gridSizeX=`readParam "$gridSizeX" "lanzSizeGridX"`
		gridSizeY=`readParam "$gridSizeY" "lanzSizeGridY"`
		gridSizeZ=`readParam "$gridSizeZ" "lanzSizeGridZ"`
		lanzTimeOut=`readParam "$lanzTimeOut" "lanzTimeOut"`
		lanzCreateResume=`readParam "$lanzCreateResumen" "lanzCreateResumen"`
	else
		txtErrror="-s .Not Available Software, see the help"
		f_help
	fi
}
#
#
showVersion()
{
	codVersion=`git rev-list --reverse HEAD |tail -1`
	numVersion=`git rev-list --reverse HEAD |wc -l`
	rama=`git branch -avv |grep "\*" |awk '{print $2}'`
	version="${GREEN}Version: ${BLUE}${numVersion} ${GREEN}Rev( ${BLUE}${codVersion}${GREEN} ) ${NONE}"
	echo ""
	echo -e "${GREEN}${version}${GREEN} Branch: ${BLUE}$rama${NONE}"
	echo ""
	exit
}
#
#	Check option is empty
#
isEmpty()
{
	if [ "$1" == "N/A" ] || [  -z "$1" ] ;then
		txtErrror="$2"
	 	f_help
	fi
}


existInLst()
{

	
	auxIFS=$IFS
	
	IFS=',' read -ra extx <<< "$1"    #Convert string to array
	aux=""

	for i in "${extx[@]}"; do
    	if [ "${2}" == "${i}" ];then
	   		aux=$i
	   	fi
	done
	IFS=$auxIFS
	echo "$aux"
}


verifyXYZ()
{
	if [[ ${option} == "BD"* ]] || [[ ${option} == "VSR" ]];then
		x=0;y=0;z=0;
	fi
	if [ -z $x ] || [ -z $y ] || [ -z $z ];then
		error=7
	fi
}



validateExtProt()
{   
	#
	#	Dos casos:
	#	1º Target does not exist
	#	2º ASGARD does not accept that target
	#
	if [ -f ${target} ];then
		
		extensionProtAux="."${target##*.}
		ext_target=`existInLst "$extensionesProt" "$extensionProtAux"`
		if [ "$ext_target" == "" ];then
			txtErrror="-s El Software no soporta esa extension"
			f_help
		fi 
	else
		txtErrror="-t .No existe el fichero de target"
		f_help
	fi
}


validate_ext_query()
{
	ext_query=""
	#
	#
	if [[ "$option" == *"VS"* ]];then
	    if [ "$extensionesLig" != "" ];then
            auxIFS=$IFS
            IFS=',' read -ra extx <<< "$extensionesLig"    		
            for i in "${extx[@]}"; do
                #echo $i
                aux=`find ${query}/ -name "*${i}" -maxdepth 1 2> /dev/null|wc -l` 

                if [ $aux != 0 ] ;then
                    ext_query=$i
                    break;
                fi
            done
            IFS=$auxIFS
            if [ "$ext_query" == "" ];then
                txtErrror="-q Debe introducir una carpeta con queries validos para el software"
                f_help
            fi
        fi
	else 
		if [ -f "$query" ];then
			extensionLigAux="."${query##*.}
			ext_query=`existInLst "$extensionesLig" "$extensionLigAux"`
			if [ "$ext_query" == "" ];then
				txtErrror="-s El Software no soporta la extension del liagndo"
				f_help
			fi
		else
			txtErrror="El fichero del query no existe"
			f_help
		fi
	fi
	name_query=$(basename "$query")

	name_query="${name_query%.*}"

}



OpcionesExtras()
{

	if [ $software == "S3" ] || [ $software == "s3" ];then 
		sed -i -r "s|Java Path=.+|Java Path=${path_external_sw}jre/|" ${path_external_sw}ChemAxon/JChem/bin/java.ini
	elif [ $option == "SD" ];then
				extension=".pdb"
	elif [ $software == "AD" ];then
		if [ "$ext_target" != "" ] && [ "$ext_target" != ".pdbqt" ];then
			source ${path_login_node}convert/mol2_pdbqt.sh
			convertReceptorTopdbqtInteractivo $ext_target
		fi
		if [ -n "$extensionLig" ]  && [ $extensionLig != ".pdbqt" ];then
			source ${path_login_node}convert/mol2_pdbqt.sh
			convertLigandTopdbqtInteractivo $extensionLig
		fi
	fi
}



find_name_job()
{
    if [ "${name_job}" == "" ]; then
        if [ "$secuencial" == "N/A" ] && [ "$command_show_jobs" != "N/A" ];then
            fckUser=$USER 					
            MAXJOBS=50 						
            for (( i=1; i<=MAXJOBS; i++ ))
            do
                name_job=${option}_${software}_$i 	#job name
                a=`$command_show_jobs -u $fckUser |grep -w $name_job |wc -l`
                if [ $a -eq 0 ];then
                    break
                fi
            done
        else #sequential
            name_job=${option}_${software}_sequential
        fi
    fi
}

#
#_______________________________________________________________________________________________________________________


if [ "$check" != "N/A" ];then
	source ${path_login_node}preDebug.sh
	exit
fi

if [ "${versionHelp}" != "N/A" ];then
	showVersion
fi

isEmpty "$software" 				"-s Sofstware is empty"
readParams											
isEmpty "$option" 					"-o Option is empty"

if [ "$extensionesLig"  != "" ];then
    isEmpty "$query" 					"-q Ligand is empty"
else
    query="no_query"
fi


if [ "$extensionesProt" != "" ];then 					
	isEmpty "$target" 			"-t Receptor is empty"
	validateExtProt									
	name_target=$(basename $target)
	name_target="${name_target%.*}"
fi
echo ${software}
if [ "${extensionesProt}" == ".mol2" ];then
  if [ ${software} == "LF" ] || [ ${software} == "FB" ]  ;then
    result_check=`python ${path_extra_shuttlemol}/used_by_shuttlemol/check_protein_mol2.py $target`
    if [ "${result_check}" != "" ];then
        echo -e "\n${result_check}\n"
        read -p "Press enter to continue"
    fi
  fi
fi


isEmpty "$num_per_job" 				"-j No ha indicado numero de jobs"
isEmpty "$protocolP" 				"-prp No ha indicado El protocolo usado para convertir la target"
isEmpty "$protocolL" 				"-prl No ha indicado El protocolo usado para convertir el/los queries"

#

#
if [ ! -f "${path_login_node}techniques/SLTechnique${option}.sh" ];then
	echo "Entro"
	txtErrror="-o La opcion ${option} no existe "
	f_help
fi


valid_option=false
OLD_IFS=${IFS}
IFS=',' read -ra ADDR <<< "$optionsLanz"
for i in "${ADDR[@]}"; do
    if [[ "${i}" == "${option}" ]];then
        valid_option=true
    fi
done
IFS=${OLD_IFS}
if [[ $valid_option == false ]];then
    txtErrror="Error option ($option) not valid for $software"
    f_help
fi





validate_ext_query


OpcionesExtras
verifyXYZ

find_name_job
#
#	Debug
#

debugB "_________________________________________Input Data________________________________________"
debugB ""
debugB "verifiInputData.sh: Mem: $mem"
debugB "verifiInputData.sh: cores: $cores"
debugB "verifiInputData.sh: time_experiment: $time_experiment"
debugB "verifiInputData.sh: x: $x"
debugB "verifiInputData.sh: y: $y"
debugB "verifiInputData.sh: z: $z"
debugB "verifiInputData.sh: grid: $grid"
debugB "verifiInputData.sh: name_target: $name_target"
debugB "verifiInputData.sh: ext_target: $ext_target"
debugB "verifiInputData.sh: name_query: $name_query"
debugB "verifiInputData.sh: ext_query: $ext_query"
debugB "verifiInputData.sh: Nomjob: $name_job"
debugB "verifiInputData.sh: Error: $error"
debugB ""


