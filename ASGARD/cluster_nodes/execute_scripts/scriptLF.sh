#!/usr/bin/env bash
#
# Definicion de las graficas
#
graph_global_field="Van_Der_Wals":"E(sol)":"E(H-bonds)":"E(elect)":"E(internal)":"E(tors)":"E(Metal)":"E(dihedral)":"E(penalty)":"dG_of_binding,_kcal/mol"
graph_global_color="b":"g":"r":"c":"m":"#eeefaa":"y":"#aaefff":"#bbefaa":"k"
graph_atoms_field="E(VDW)":"E(sol)":"E(H-bonds)":"E(elect)":"E(internal)":"E(metal)":"E(total)"
graph_atoms_color="b":"g":"r":"c":"m":"#eeefaa":"#bbefaa":"k"


execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)` 												#En todos los comandos debe existir esta etiqueta
if [ ${nom_serv} ==  "leftraru4" ];then
		export HOME=/mnt/flock/hperez
	fi
	#
	#	En LF existe un prametro especial para el refiunado de score mediante dinamica, o normal
	#
	OLDIFS=$IFS
	IFS='-' read -r -a array <<< "$opt_aux"
	refined_energy=""
	for element in "${array[@]}"
	do
		aux=`echo $element | cut -f1 -d' ' | tr a-z A-Z`
		param=`echo $element | cut -f2 -d' ' | tr a-z A-Z`
		case "$aux" in
			REFINED_ENERGY)	refined_energy=$param;;
		esac
	done
	IFS=$OLDIFS;

	if [ $refined_energy != "LF" ] && [ $refined_energy != "GR" ];then
		  echo "ERROR: bad -refined_energy $refined_energy"
		  exit
	fi

	opt_aux=`echo "${opt_aux/-refined_energy $refined_energy/}"` #removemos la variable
	case $option in
		VS)
            out_aux=${option}_${software}_${name_target}_${name_query}_${x}_${y}_${z}
            out_grid=${folder_grid}${out_aux}
            coords=$x":"$y":"$z
			execute "${path_external_sw}leadFinder/leadfinder -q --load-grid=${out_grid}_grid.bin\
			--protein=${CWD}${target} --ligand=${CWD}${query} --output-tabular=${out_energies}.log \
			--output-poses=${out_molec}.mol2 --text-report=${out_energies}.eng --max-poses=${numPoses}   $opt_aux &>  $out_aux.ucm "
			create_out
		;;
		BD  )
			execute "source ${path_cluster_nodes}generate_grid.sh"                  	#siempre se genera grid
			execute "${path_external_sw}leadFinder/leadfinder -q --load-grid=${out_grid}_grid.bin \
			--protein=${CWD}${target} --ligand=${CWD}${query} --output-tabular=${out_energies}.log \
			--output-poses=${out_molec}.mol2 --text-report=${out_energies}.eng --max-poses=${numPoses}   $opt_aux &>  $out_aux.ucm "
            execute "coords=\"`${path_extra_shuttlemol}used_by_shuttlemol/get_center_ligand.py ${out_molec}.mol2`\""
			create_out
			execute "rm ${out_grid}_grid.bin"
		;;
	esac

}
function create_out()
{
  	if [ $refined_energy  == "GR" ];then
	    refined_dm
    fi
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

refined_dm()
{
  echo "python2.7 ${path_extra_shuttlemol}/used_by_shuttlemol/get_energy_gromcas.py ${target} ${out_molec}.mol2 "
	python2.7 ${path_extra_shuttlemol}/used_by_shuttlemol/get_energy_gromcas.py ${target} ${out_molec}.mol2 > ${out_energies}".enmd"
  aux=`cat ${out_energies}".enmd" |tail -1`
	global_score_md=`echo $aux | cut -d ":" -f 2`
}