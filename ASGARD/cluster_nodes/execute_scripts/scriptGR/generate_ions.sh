#!/usr/bin/env bash

function ejecutar()
{
	execute "echo \"GenerateIons\""
	add_ions
    execute "echo Fin GenerateIons"
}
#
#	Añade los iones
#
function add_ions()
{
    name_solv="CA"
	execute "${prefix_gromacs} grompp${mpi} -f ${out_grid}_tpr.mdp  -c ${out_molec}_complex_solv.gro -p ${out_molec}_complex.top -o ${out_aux}_ions.tpr -po ${out_aux}_grompp_IONS.mdp -maxwarn ${MAX_WARNINGS}" #remplaza iones por atomos de agua tantos sean necesarios para estabilizar la carga en la target
	num_ions=${system_charge%.*}

	if [ "$num_ions" == "0" ];then 
	    execute "cp ${out_molec}_complex_solv.gro ${out_molec}_complex_solv_ions.gro"
		execute "mv ${out_molec}_complex_solv.gro ${folder_out_ucm}"
	else

		if [ $num_ions -gt 0 ];then
		    echo "System_charge:  $system_charge add negative ions  $num_ions CL"
			execute "echo ${solvatation} |${prefix_gromacs} genion${mpi} -s ${out_aux}_ions.tpr -o ${out_molec}_complex_solv_ions.gro -p ${out_molec}_complex.top -nname CL  -nn $num_ions -seed ${seedg}" # -np numero ions positivos
		elif [ $num_ions -lt 0 ];then
 		    num_ions=`expr $num_ions \* -1`
 		    echo "System_charge:  $system_charge add positive ions  $num_ions NA"
			execute "echo ${solvatation} | ${prefix_gromacs} genion${mpi} -s ${out_aux}_ions.tpr -o ${out_molec}_complex_solv_ions.gro -p ${out_molec}_complex.top -pname NA  -np $num_ions -seed ${seedg}" # -nn num ioons negativos
        else
            echo "Erro Añadiendo ioenes"
            exit
		fi
		execute "mv ${out_molec}_complex_solv.gro ${folder_out_ucm}"
	fi

}
