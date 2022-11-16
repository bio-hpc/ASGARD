#!/bin/bash
# 	Author: Jorge de la Peña García
#	Email:  jorge.dlpg@gmail.com
#	Description: Se encarga de llamar a los scripts de Software correspondiente
#_________________________________________________________________________________________________________________________________________________________________________

#
#	ejecuta el comando y si debug esta activo ejecuta el debug
#
#TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
TAG=""
DELIMITER="\\n"
#
#	Si ocurre algun error esta funcion hara un pequeño resumen del error lo guardara en ${directorio}ERRORS-${name_target}.errResu
#
execute()
{
	debugY "$TAG: $1"
	eval "$1"
	error=$? 															#avisa si hay error en los comandos, se usa en los scriptAD.sh .....
}
#
#	Estandariza la salida a json
#
function standar_out_file()
{
	out_put=${name_query}${DELIMITER}       	 				#	Nombre query                            (automatic)
	out_put=${out_put}${query}${DELIMITER}					    # 	Fichero query                           (automatic)
	out_put=${out_put}${target}${DELIMITER}					    # 	Fichero query                           (automatic)
	out_put=${out_put}${file_result}${DELIMITER}				#	molecula de salida                      (need) (file_result=${out_molec}.pdbqt)
	out_put=${out_put}${coords}${DELIMITER}						#	Coordenadass                            (need in BD) (X:Y:Z) get_center_ligand.py
	out_put=${out_put}${number_execution}${DELIMITER}			#	Numero de ejecucion                     (automatic)
	out_put=${out_put}${num_amino_acid}${DELIMITER}				# 	Numero de aminoacido                    (automatic)
	out_put=${out_put}${chain}${DELIMITER}						#	Cadena de la target                     (automatic)
	out_put=${out_put}${global_score}${DELIMITER}				#	Score global                            (need)
	out_put=${out_put}${gridSizeX}","${gridSizeY}","${gridSizeZ}${DELIMITER} #				                (automatic)
	out_put=${out_put}${option}${DELIMITER}                     #	opcion ("VS|BD..")                      (automatic)
	out_put=${out_put}${software}${DELIMITER}                   #	software ("AD|MT...")                   (automatic)
	out_put=${out_put}${global_score_md}${DELIMITER}            #	score md para programas de dock
	out_put=${out_put}${global_score_qu}${DELIMITER}            #	score md para programas de dock
	out_put=${out_put}${DELIMITER}                  # 2 free slots
	out_put=${out_put}${graph_global_color}${DELIMITER}			#	Colores para score global               (1º optional)
	out_put=${out_put}${graph_global_field}${DELIMITER}			#	Campos para Score global                (2º optional)
	out_put=${out_put}${graph_global_score}${DELIMITER}			#	Score global                            (3º optional)
	out_put=${out_put}${graph_atoms_color}${DELIMITER}			#	Colores para score Atom                 (1º optional)
	out_put=${out_put}${graph_atoms_field}${DELIMITER}			#	Campos para score atoms                 (2º optional)
	out_put=${out_put}${graph_atoms_type}${DELIMITER}			#	Campos para tipo atoms                  (3º optional)
	out_put=${out_put}${graph_atoms_score}${DELIMITER}			#	Score por atomos                        (4º (optional))	
	python ${path_cluster_nodes}standar_out_put.py "${out_energies}.json" "$out_put"
}



#
#	Main
#
if [ -z ${check} ];then 												# si no esta en modo checkeo se cargan los parametros
	source ${CWD}ShuttleMol/login_node/read_all_conf.sh
fi
execute "re='^[0-1]+([.][0-9]+)?$'" 									#empiezan por 0-1 para comprobar scores
#
#   Coordenadas diferentes a 0 normalmente software de  docking
#       Si el numero de ejecucion es -1 = VS
#       Si el numero de ejecucion es != -1 = BD
#  	Coordenadas = 0 software de similaridad 
#		
#
if [ "$x" != "0" ] || [  "$y" != "0" ] || [ "$z" != "0" ];then
    if [ "$number_execution" == "-1" ];then
        aux_query=`basename $query`
        aux_query="${aux_query%.*}"
        execute "out_prefix=${option}_${software}_${name_target}_${aux_query}_${x}_${y}_${z}"
    else
        execute "out_prefix=${number_execution}_${option}_${software}_${name_target}_${name_query}_${x}_${y}_${z}"
    fi

else
	: 'Comentado por fallo en LS 13/12/2018
    #aux_query=`basename $query`
    #aux_query="${aux_query%.*}"
    if [ "$number_execution" == "-1" ];then
        execute "out_prefix=${option}_${software}_${name_target}_${aux_query}"
    else
        execute "out_prefix=${number_execution}${option}_${software}_${name_target}_${name_query}"
   #fi
   '
   if [ "$software" == "OP" ];then
         echo "entro a ops"
        name_query=$(basename $query)
        name_query=${name_query%.*}
  else
      aux_query=`basename $query`
      aux_query="${aux_query%.*}"
   fi
   if [ "${option}" == "VSB" ];then
   		aux_query=${aux_query}_${ini}_${fin}   	   	
   fi
   execute "out_prefix=${option}_${software}_${name_target}_${aux_query}"
fi

#out_prefix=`echo ${out_prefix%.*}`


#
#	Carpetas de salida para los scripts de ejecucion
#
global_score_md=""
global_score_qu=""
opt_aux=`cat ${folder_templates_jobs}parameter_aux.txt` #lee si se ha introducid algun aprametro extra
out_grid=${folder_grid}${out_prefix}
out_aux=${folder_out_ucm}${out_prefix}
out_molec=${folder_molec}${out_prefix}
out_energies=${folder_energies}${out_prefix}
execute "source ${path_cluster_nodes}execute_scripts/script${software}.sh"
execute "execute_script"


