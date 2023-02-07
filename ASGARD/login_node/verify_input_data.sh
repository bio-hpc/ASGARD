#!/usr/bin/env bash
#
#   Description: Codigo encargado de validar las optiones necesiarias para el funcionamiento de shuttlemol, receptor, query, jobs, protocolos
# ______________________________________________________________________________________________________________________

#
#	Lee los parametros de la plantilla si el usuario no indica otra cosa los coge por defecto
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
#	Muestra la version y la rev
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
#
#	Busca en una array ($1)si se encuentra una palabra($2)
#
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
#
#	Verifica que las x y z o numJobs no este vacia
#
verifyXYZ()
{
	if [[ ${option} == "BD"* ]] || [[ ${option} == "VSR" ]];then
		x=0;y=0;z=0;
	fi
	if [ -z $x ] || [ -z $y ] || [ -z $z ];then
		error=7
	fi
}

#
#	Valida si la extension de la target es valida poara SW
#
validateExtProt()
{   
	#
	#	Dos casos:
	#	1º No existe target
	#	2º EL software no acepta esa target
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
#
#	VAlida la extension del query
#
validate_ext_query()
{
	ext_query=""
	#
	#	2 optiones o VS o BD
	#		Si es VS se busca en la carpeta de queries si hay alguno con las extensiones indicadas en la plantilla para el software
	#		Si es BD se comprueba que el software admita esa extension
	#
	if [[ "$option" == *"VS"* ]];then
	    if [ "$extensionesLig" != "" ];then
            auxIFS=$IFS
            IFS=',' read -ra extx <<< "$extensionesLig"    			#recorremos los posibles extensiones buscando una que no de 0
            for i in "${extx[@]}"; do
                #echo $i
                aux=`find ${query}/ -name "*${i}" -maxdepth 1 2> /dev/null|wc -l` #faltaria hacer caso especifico para cuando es VS

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
	else #ocpiones BD
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
#
#	Opciones extras para alguno s softwares
#
OpcionesExtras()
{

	if [ $software == "S3" ] || [ $software == "s3" ];then #problemas con java
		sed -i -r "s|Java Path=.+|Java Path=${path_external_sw}jre/|" ${path_external_sw}ChemAxon/JChem/bin/java.ini
	elif [ $option == "SD" ];then
				extension=".pdb"
	elif [ $software == "AD" ];then
		#echo $ext_target
		if [ "$ext_target" != "" ] && [ "$ext_target" != ".pdbqt" ];then
			source ${path_login_node}convert/mol2_pdbqt.sh
			convertReceptorTopdbqtInteractivo $ext_target
		fi
		#la extension puede no existir si es un VS
		if [ -n "$extensionLig" ]  && [ $extensionLig != ".pdbqt" ];then
			source ${path_login_node}convert/mol2_pdbqt.sh
			convertLigandTopdbqtInteractivo $extensionLig
		fi
	fi
}
#	
#	Busca un nombre de job
#
find_name_job()
{
    if [ "${name_job}" == "" ]; then
        if [ "$secuencial" == "N/A" ] && [ "$command_show_jobs" != "N/A" ];then
            fckUser=$USER 					# nombre de usuario para ver los jobs
            MAXJOBS=50 						#numero de nombres maximo de jobs
            for (( i=1; i<=MAXJOBS; i++ ))
            do
                name_job=${option}_${software}_$i 	#nombre identificativo del job
                a=`$command_show_jobs -u $fckUser |grep -w $name_job |wc -l`
                if [ $a -eq 0 ];then
                    break
                fi
            done
        else #si es secuencial
            name_job=${option}_${software}_sequential
        fi
    fi
}

#
#	MAAIN								comprobacion de datos basicos de entrada, no chequeamos si son correctos solo que existen
#_______________________________________________________________________________________________________________________

#	1º Modo check de software, manda un job donde se ejecutara una vez cada software y comprobara que funciona Todavia no funciona del todo
if [ "$check" != "N/A" ];then
	source ${path_login_node}preDebug.sh
	exit
fi
#	2º primero se comprueba si a solicitado la version
if [ "${versionHelp}" != "N/A" ];then
	showVersion
fi
#
# datos necesarios para una ejecuion minima de shuttlemol
#
isEmpty "$software" 				"-s Sofstware is empty"
readParams												#leemos los parametros de la plantilla, si no existirra plantilla daria erro
isEmpty "$option" 					"-o Option is empty"

if [ "$extensionesLig"  != "" ];then
    isEmpty "$query" 					"-q Ligand is empty"
else
    query="no_query"
fi


if [ "$extensionesProt" != "" ];then 					#algunos softwares no requieren target, si en la plantilla no hay indicado extension no se usara
	isEmpty "$target" 			"-t Receptor is empty"
	validateExtProt										#validamos extension de la target en caso de que no exista
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
#	Validamos Si exsite la option
#
if [ ! -f "${path_login_node}techniques/SLTechnique${option}.sh" ];then
	echo "Entro"
	txtErrror="-o La opcion ${option} no existe "
	f_help
fi
#
#   Validamos si el software permite esa opcion
#
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




#
#	Validamos la extension del query
#
validate_ext_query
#
#	Opciones extras para softwares y coordenadas xyz
#
OpcionesExtras
verifyXYZ

find_name_job
#
#	Opciones de debug
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


