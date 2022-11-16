#!/usr/bin/env bash
#___________________________________________________________________________________
#
#	Hace docking con Autodock4
#	1º Genera la grid tanto para VS como para BD por que depende del query
#	2º hace docking
#___________________________________________________________________________________

graph_global_field="vdw:elec:ds:Free_energy"
graph_global_color="g:r:b:y"
graph_atoms_field="vdw:elec:ds:Free_energy"
graph_atoms_color="g:r:b:y"


execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
	execute "source ${path_cluster_nodes}generate_grid.sh"
	if [ "$error" == "0" ];then 																								#si la ejecucion es correcta se continua para obtener el score
		execute "${path_external_sw}autodock4/autogrid4 -p ${out_aux}_grid.gpf -l ${out_grid}_grid.gpf.out 2> ${out_aux}.ucm &"  				#se genera la grid

		
		source ${path_cluster_nodes}grids/spinner.sh
		if [ "$error" == "0" ];then 																							#si la ejecucion es correcta se continua para obtener el score
			execute "${path_external_sw}autodock4/autodock4 -p ${out_aux}_dock.dpf -l ${out_molec}.dock $opt_aux 2> ${out_aux}.ucm"
			if [ "$error" == "0" ];then #si la ejecucion es correcta se continua para obtener el score
				execute "conform=\"`cat ${out_molec}.dock|grep -A 3 Rank | head -3 |tail -1 |awk -F\| '{print $3}' |tr -d ' '`\"" #busco la mejor conformacion seg´un el rank de AK
				execute "${path_external_sw}mgltools/bin/pythonsh ${path_external_sw}mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/write_conformations_from_dlg.py  -d ${out_molec}.dock -o ${out_molec}"
				if [ "$error" == "0" ];then
					execute "cp ${out_molec}_${conform}.pdbqt ${out_molec}.pdbqt"
					execute "rm ${out_molec}_*.pdbqt"
					execute "${path_external_sw}mgltools/bin/pythonsh ${path_external_sw}mgltools/MGLToolsPckgs/AutoDockFR/bin/ade  -lig ${out_molec}.pdbqt -rec ${CWD}${target} -gridMaps ${out_aux}*.map > ${out_energies}.eng"
					if [ "$option" == "BD" ] ;then
					    execute "coords=\"`${path_extra_shuttlemol}used_by_shuttlemol/get_center_ligand.py ${out_molec}.pdbqt`\""
					else
					    coords=${x}:${y}:${z}
					fi
                    create_out
					#Se borran todos los ficheros creados para la prueba
					##execute "rm ${out_aux}*.map "
					##execute "rm ${out_grid}_grid.gpf* "
					##execute "rm ${out_aux}.maps.xyz "
					##execute "rm ${out_aux}.maps.fld "
					##execute "rm ${out_molec}.dock "
					##execute "rm ${out_aux}_dock.dpf"
				fi
			fi
		fi
	fi
}
function create_out()
{
    file_result=${out_molec}.pdbqt
    execute "global_score=\"`cat ${out_energies}.eng | grep \"Free energy of binding in kcal/mol\"| awk '{print  $NF }'`\""
    graph_global_score
    graph_atom_score
    execute "standar_out_file" #normaliza la salida en un jsocn
}
function graph_global_score()
{

    execute "graph_global_score=\"`cat ${out_energies}.eng| grep "Sum:" |awk '{print $10\":\"$11\":\"$12\":\"$9}'`\""
}
function graph_atom_score()
{
    start_ener=`grep -n Num ${out_energies}.eng |head -1 | awk -F: '{print $1}'`
    start_ener=`expr $start_ener + 1`
    end_ener=`grep -n "=========="  ${out_energies}.eng| tail -1 |awk -F: '{print $1}'`
    end_ener=`expr $end_ener - 1`
    graph_atoms_score=`sed -n "$start_ener,$end_ener p" ${out_energies}.eng |  awk '{ print $10":"$11":"$12":"$10+$11+$12   "\\\n"}'`
    graph_atoms_score=`echo $graph_atoms_score |sed 's/\ //g'`
    graph_atoms_type=`sed -n "$start_ener,$end_ener p"  ${out_energies}.eng | awk '{print $8"_"$9":"}'`
	graph_atoms_type=`echo $graph_atoms_type |sed 's/\ //g'`
}