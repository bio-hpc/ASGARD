#!/usr/bin/env bash

graph_global_field="Electrostatic":"Van_der_Waals":"Hydrogen_bonds":"Scores"
graph_global_color="b:g:r:c"
graph_atoms_field="Electrostatic":"Van_der_Waals":"Hydrogen_bonds"
graph_atoms_color="b:g:r"

execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
	ulimit -s unlimited &>/dev/null     #util para para incremenar la pila pero no se sule poder usar por que requiere permisos (puede causar error)
	outDock=${folder_grid}${out_prefix} 	#al ser necesario fichero de docking se crea una variable nueva que luego al generar grids se usar #en caso de BD es la misma que out_grid
	case $option in
		VS) #virtualscreening
		    coords=${x}:${y}:${z}

		    out_aux=${option}_${software}_${name_target}_${name_query}_${x}_${y}_${z}
            execute "out_grid=${folder_grid}${out_aux}"
            coords=$x":"$y":"$z

			execute "source ${path_cluster_nodes}generate_grid.sh"	 																		#Genera grid y ficehro de dock si es VS la grid no se genera
			execute "${path_external_sw}"FlexScreen_Biotin/dock" ${outDock}_dock.inp &> ${out_aux}.ucm"
			if [ $? == 0 ];then
				execute "cat ${out_molec}.mol|tr -d \"\r\" > $out_molec.mol2" 		#se convierte a mol2
                create_out
			fi
			;;
		BD | BDC | VSBD) #blindDocking
			execute "source ${path_cluster_nodes}generate_grid.sh"	 																		#Genera grid y ficehro de dock si es VS la grid no se genera
			execute "${path_external_sw}FlexScreen_Biotin/dock ${outDock}_dock.inp &> ${out_aux}.ucm"
			if [ $? == 0 ];then
				execute "cat ${out_molec}.mol|tr -d \"\r\" > $out_molec.mol2" 		#se convierte a mol2
				execute "coords=\"`${path_extra_shuttlemol}used_by_shuttlemol/get_center_ligand.py ${out_molec}.mol2`\""
				create_out
				#   Borrar data grids
				execute "rm ${outDock}_dock.inp"
				execute "rm ${out_grid}_g.mol2"
				execute "rm ${out_grid}_grid.inp"
				execute "rm ${out_grid}_l.es"
				execute "rm ${out_grid}_r.mol2"
				execute "rm ${out_grid}_s.vdw "
			fi
			;;
	esac
}
function create_out()
{
    execute "rm ${out_molec}.mol" #se elimina el fichero mol
    file_result=${out_molec}.mol2
    execute "global_score=\"`cat ${out_energies}.dat |head -1 |awk '{print $7}'`\""
    graph_global_score
    graph_atom_score
    execute "standar_out_file" #normaliza la salida en un jsocn
}
function graph_global_score()
{
    execute "graph_global_score=\"`cat ${out_energies}.eng | grep  ^\- |head -1 |awk '{print $1\":\"$2\":\"$3\":\"$4}'`\""
}
function graph_atom_score()
{
    start_ener=`grep -n \  ${out_energies}.eng  |head -1  | awk -F: '{print $1}'` #OJO no hay que sumarle 1
    #start_ener=`expr $start_ener + 1`
    end_ener=`grep -n \^-  ${out_energies}.eng | head -1 |awk -F: '{print $1}'`
    end_ener=`expr $end_ener - 1`
    graph_atoms_score=`sed -n "$start_ener,$end_ener p" ${out_energies}.eng |  awk '{ print $2":"$3":"$4   "\\\n"}'`
    graph_atoms_score=`echo $graph_atoms_score |sed 's/\ //g'`
    graph_atoms_type=`sed -n "$start_ener,$end_ener p"  ${out_energies}.eng | awk '{print $1":"}'`
	graph_atoms_type=`echo $graph_atoms_type |sed 's/\ //g'`
}