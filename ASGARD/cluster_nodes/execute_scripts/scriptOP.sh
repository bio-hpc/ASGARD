#!/usr/bin/env bash
execute_script()
{
	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
	lib=`echo $(basename $query)`
    opt_aux=${opt_aux/-SS/-S} #como el parametreo -s esta asignado a sw se hace un mapeo a -SS y aqui se vuelve a mapear a -S
    #run_optifarm ${CWD}${target}  ${CWD}${query} "${opt_aux}" "${out_energies}.en" 





	mol_name=`cat $query |grep -A1 MOLNAME |tail -1`
	mol_name=${mol_name// /_}



	execute "out_prefix=${mol_name}"
	out_molec=${folder_molec}${out_prefix}

	#execute "ShuttleMol/external_sw/ChemAxon/JChem/bin/molconvert -3:S{fine}[mmff94]L{3} mol2 $query -o ${out_molec}.mol2 -F"
	#sed  -i "2s/.*/${mol_name}/" ${out_molec}.mol2



	target=projects/b_maigret/queries/Q1.mol2
	name_target=Q1
	execute "out_prefix=${option}_${software}_${name_target}_${mol_name}"
	out_energies=${folder_energies}${out_prefix}
    run_optifarm ${CWD}${target}  ${out_molec}.mol2 "${opt_aux}" "${out_energies}.en"

    target=projects/b_maigret/queries/Q2.mol2
	name_target=Q2
	execute "out_prefix=${option}_${software}_${name_target}_${mol_name}"
	out_energies=${folder_energies}${out_prefix}
    run_optifarm ${CWD}${target}  ${out_molec}.mol2  "${opt_aux}" "${out_energies}.en"

    target=projects/b_maigret/queries/Q3.mol2
	name_target=Q3
	execute "out_prefix=${option}_${software}_${name_target}_${mol_name}"
	out_energies=${folder_energies}${out_prefix}
    run_optifarm ${CWD}${target}  ${out_molec}.mol2  "${opt_aux}" "${out_energies}.en"


}
run_optifarm()
{
    execute "cd ${path_external_sw}/OptiPharm/"
    execute "./OPShapeSimilarity_centos  -q ${1} -d ${2}  ${3} |grep -v  optipharm: > ${4}"

    score=`cat  ${4} |awk '{print $NF}'`
    score=`echo $score |awk '{printf "%3.8f\n", $1}'`
    echo $score
    echo $(echo "$score < 0.5" |bc -l)
    if [ $(echo "$score < 0.5" |bc -l) == 1 ]; then
    	echo "entro a rm"
    	rm ${4}
    fi
}







