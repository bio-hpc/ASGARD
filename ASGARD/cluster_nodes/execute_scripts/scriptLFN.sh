#!/usr/bin/env bash
#
# Definicion de las graficas
#
graph_global_field="Van Der Wals":"E(sol)":"E(H-bonds)":"E(elect)":"E(internal)":"E(tors)":"E(Metal)":"E(dihedral)":"E(penalty)":"dG_of_binding,_kcal/mol"
graph_global_color="b":"g":"r":"c":"m":"#eeefaa":"y":"#aaefff":"#bbefaa":"k"
graph_atoms_field="E(VDW)":"E(sol)":"E(H-bonds)":"E(elect)":"E(internal)":"E(metal)":"E(total)"
graph_atoms_color="b":"g":"r":"c":"m":"#eeefaa":"#bbefaa":"k"


execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)` 												#En todos los comandos debe existir esta etiqueta
    lead_finder_run=${path_external_sw}lead_finder_2018/leadfinder_x86_64
	case $option in
		VS)
            out_aux=${option}_${software}_${name_target}_${name_query}_${x}_${y}_${z}
            out_grid=${folder_grid}${out_aux}
            coords=$x":"$y":"$z
			execute "${lead_finder_run} -q --load-grid=${out_grid}_grid.bin\
			--protein=${CWD}${target} --ligand=${CWD}${query} --output-tabular=${out_energies}.log \
			--output-poses=${out_molec}.mol2 --text-report=${out_energies}.eng --max-poses=${numPoses}   $opt_aux &>  $out_aux.ucm "
			create_out
		;;
		BD  )
			execute "source ${path_cluster_nodes}generate_grid.sh"                  	#siempre se genera grid
			execute "${lead_finder_run} -q --load-grid=${out_grid}_grid.bin \
			--protein=${CWD}${target} --ligand=${CWD}${query} --output-tabular=${out_energies}.log \
			--output-poses=${out_molec}.mol2 --text-report=${out_energies}.eng --max-poses=${numPoses}   $opt_aux &>  $out_aux.ucm "
            execute "coords=\"`${path_extra_shuttlemol}used_by_shuttlemol/get_center_ligand.py ${out_molec}.mol2`\""
			create_out
			#execute "rm ${out_grid}_grid.bin"
		;;
	esac

}
function create_out()
{
    file_result=${out_molec}.mol2
    execute "global_score=\"`cat ${out_energies}.log |grep '^/' |awk '{print \$3}'`\""
    graph_global_score
    graph_atom_score
    execute "standar_out_file" #normaliza la salida en un jsocn
}
function  graph_global_score
{
    g_ener=(`cat ${out_energies}.eng |grep -A 10 "Total Energy components :" |tail -10`)
    #"Van Der Wals":"E(sol)":"E(H-bonds)":"E(elect)":"E(internal)":"E(tors)":"E(Metal)":"E(dihedral)":"E(penalty)":"dG of binding, kcal/mol"
    graph_global_score=${g_ener[1]}:${g_ener[5]}:${g_ener[7]}:${g_ener[9]}:${g_ener[11]}:${g_ener[13]}:${g_ener[3]}:${g_ener[15]}:${g_ener[17]}:${g_ener[22]}
}
function graph_atom_score
{
    start_ener=`grep -n Name ${out_energies}.eng |head -1 | awk -F: '{print $1}'`
    start_ener=`expr $start_ener + 1`
    end_ener=`grep -n Total  ${out_energies}.eng| head -1 |awk -F: '{print $1}'`
    end_ener=`expr $end_ener - 1`
    #"E(VDW)":"E(sol)":"E(H-bonds)":"E(elect)":"E(internal)":"E(metal)":"E(total)"
    graph_atoms_score=`sed -n "$start_ener,$end_ener p" ${out_energies}.eng |awk '{print $5":"$7":"$8":"$9":"$10":"$11 "\\\n"}'`
    graph_atoms_score=`echo $graph_atoms_score |sed 's/\ //g'`
    graph_atoms_type=`sed -n "$start_ener,$end_ener p"  ${out_energies}.eng | awk '{print $1"_"$2":"}'`
	graph_atoms_type=`echo $graph_atoms_type |sed 's/\ //g'`
}

#
#	BDVS )
#
#			check=$ligVS
#            execute "ligVS=`echo $ligVS | sed 's/.$//g'`" 																			#LigVs indica el  el directorio donde estan los queries, se usa para VS usar la misma grid
#			execute "ligVS=${ligVS##*/}"
#			execute "out_grid=${directorio}grids/${option}"_"${software}"_"${name_target}_${ligVS}_${x}_${y}_${z}" 					#indica donde esta la grid 1 por coordenada
#			execute "source ${pathSD}generate_grid.sh"																			#si no existiera la grid por algun casual la crearia
#			execute "${pathSW}leadFinder/leadfinder -q --load-grid=${out_grid}grid.bin \
#                     -np ${cores} --protein=${CWD}$target --ligand=$query --output-tabular=${out_energies}.log \
#                     --output-poses=${out_molec}.mol2 --text-report=${out_energies}.eng --max-poses=${numPoses}  $opt_aux >  $out_aux.ucm 2>> $out_aux.ucm"
#            execute "a=\"`cat ${out_energies}".log" |grep '^/'|head -1| awk -F"," '{print "'$out_molec' " $3" '$query' '$numeroEjecucion' " $2}'`\""
#            execute "python ${pathSD}getCOM.py $a 0 0"
#            numOriMole=`ls ${CWD}${check}/*.mol2 |wc -l`;
#            babel -m ${out_molec}.mol2 ${out_molec}_.mol2
#			numDestMole=`find  ${directorio}energies/ -name "*${x}_${y}_${z}.log" |wc -l`
#			#
#			#	esto funciona bien en secuencial
#			#
#			if [ $numOriMole -eq $numDestMole ];then #si los ficheros de enerigias son iguales a los ficheros de queries se borra la grid (Se da por supuesto que esa grid no se va a usar mas)
#
#				execute "rm ${out_grid}grid.bin"
#			fi
#		;;