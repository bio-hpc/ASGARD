#!/usr/bin/env bash
#
#   Description: ASGARD parameters
# ______________________________________________________________________________________________________________________


#______________________________________________________________________________________________________________________
#
#			Parameters
#_________________________________________________________________________________________________________________________
read_file_template()
{
	source ${pathSL}special_params.sh
	read_template_params $software

	if [ "$error" == "0" ];then

		for i in `seq 1 $contadorFile`;do
			error=12
			IFS='::' read -ra ADDR <<< "${file[$i]}"
			if [ "$1" == "${ADDR[4]}" ];then
				if [ "${ADDR[2]}" == "Y" ];then
					error=0		
					allComand=`echo "$allComand" |sed -e "s,${ADDR[4]} $2,,g"` #remove all_comand
					optAdicionals="$optAdicionals ${ADDR[4]} $2"
					shift
					break
				else
					error=0
					if [[ "$1" != \-* ]];then #for - parameters
						error=12 # the parameters must start with -
					fi
					allComand=`echo "$allComand" |sed -e "s,${ADDR[4]},,g"` #remove all_command
					optAdicionals="$optAdicionals ${ADDR[4]}"
					break
				fi
			fi
		done
	fi

}
function empty_variable()
{
	if [ -z $2 ];then
		val="N/A"
	else
		val=$2
	fi
    if [ -z "$1" ];then
        echo $val
    else
        echo $1
    fi
}

while (( $# )) #browse all the parameters and they are assigned
 do
    if [[ "$1" == \-[a-z]* ]] || [[ "$1" == \-[A-Z]* ]] || [[ "$1" == \-\-[a-a]* ]] || [[ "$1" == \--[A-Z]* ]];then 
	   case `printf "%s" "$1" | tr '[:lower:]' '[:upper:]'`  in
      -F)   folder_analysis=${2%/};; # Folder where MD results are found
      -P)   profile=$2;; # Profile (TARGET, TARGET_PROTEIN, TEST...)
      -L)   ligand=$2;; # Ligand (Indicate what is the ligand which you want to analyze
      -D )  folder_experiment=${2%/}/;;	#directorio where the data is stored (optional)
	    esac
	fi
  shift #Modify $2 a $1, $3 a $2 ....
done

if [[ "${profile}" == "" ]];then
	profile="STANDAR_${option}"		
fi

if [[ ${software} == "GR" ]] || [[ ${software} == "GRN" ]] || [[ ${software} == "GRC" ]] ;then
	
    if [[ "$query" == "no_query" ]] || [[ "$target" == "no_target" ]] || [[ "$query" == "" ]];then
    	
        if [[ $target != "no_target" ]] || [[ $target != "" ]];then
        	
        	aux="${target%.*}"
            mode_gr=`cat ${aux}*.conf |head -1 |grep  Profile | awk -F\: '{print $2}'`

        else
            mode_gr=`cat ${query}*.conf |head -1 |grep  Profile | awk -F\: '{print $2}'`
        fi
    fi
    
 
    if [[ "${mode_gr}" == "QUERIES" ]] || [[ "${mode_gr}" == "BIPHSIC_SYSTEMS" ]];then
        target="targets/test/queries.pdb"
        if [[ ! -f $target ]];then
            echo "ERROR: Para el modo QUERIES | BIPHSIC_SYSTEMS"
            echo "ERROR: debe existir un fichero pdb placebo en "$target
            exit
        fi
    elif [[ "${mode_gr}" == "TARGET" ]];then
	
        query="ASGARD/external_sw/gromacs/config_files/target/"
        aux=`ls ShuttleMol/external_sw/gromacs/config_files/target/*.mol2 |wc -l`
        if [[ ${aux}  != 1 ]];then
            echo "ERROR: Para el modo target"
            echo "ERROR: debe existir solo un fichero .mol2  "+query
            exit
        fi
    fi
fi




