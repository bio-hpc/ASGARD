#!/usr/bin/env bash
#
#

ejecutar()
{
    equilibration_nvt
    equilibration_npt
}
#NVT
equilibration_nvt() #calentar la conformacion
{
	# We restart the processes just in case, they should be restarted
	change_porse_files 5000 5000 5000 5000

    # option -r is added for gromacs 2018
	execute "${prefix_gromacs} grompp${mpi} -f ${out_grid}_nvt.mdp -c ${out_molec}_complex_min.gro -r ${out_molec}_complex_min.gro -p ${out_molec}_complex.top -o ${out_molec}_complex_nvt.tpr -n  ${out_grid}_complex_index.ndx -maxwarn $MAX_WARNINGS"

	execute "${prefix_gromacs} mdrun${mpi} -deffnm ${out_molec}_complex_nvt ${gpu} ${thrads} -g ${out_aux}_complex_mdrun_nvt.log"
	execute "mv ${out_molec}_complex_nvt.edr ${folder_energies}"
	execute "mv ${out_molec}_complex_nvt.tpr ${folder_out_ucm}"

	execute "mv ${out_molec}_complex_min.gro ${folder_out_ucm}"



}
#NPT
equilibration_npt() #se ejecuta 4 veces incrementando la posicion de restriccion
{
	ini=${step_nvt}
	posre_fin=4000
	for ((j=1; j<5; j++))
	do
	    change_porse_files ${posre_fin} ${posre_fin} ${posre_fin}
		posre_fin=`expr $posre_fin - 1000`

		execute "$prefix_gromacs grompp${mpi} -f ${out_grid}_npt_${ini}.mdp -c ${out_molec}_complex_nvt.gro -r ${out_molec}_complex_nvt.gro -t ${out_molec}_complex_nvt.cpt -p ${out_molec}_complex.top -o ${out_molec}_complex_npt_${ini}.tpr  -n ${out_grid}_complex_index.ndx -maxwarn $MAX_WARNINGS"
		
		execute "$prefix_gromacs mdrun${mpi} -deffnm ${out_molec}_complex_npt_${ini} ${gpu} ${thrads} -g ${out_aux}_comlex_npt_${ini}.log"

		execute "mv ${out_molec}_complex_npt_${ini}.edr ${folder_energies}"
		execute "mv ${out_molec}_complex_npt_${ini}.tpr ${folder_out_ucm}"
		execute "mv ${out_molec}_complex_npt_${ini}.trr ${folder_out_ucm}"

		if [ "$j" != "4" ];then
		    execute "mv ${out_molec}_complex_npt_${ini}.cpt ${folder_out_ucm}"
			execute "mv ${out_molec}_complex_npt_${ini}.gro ${folder_out_ucm}"
		else
			last_equlibration=${out_molec}_complex_npt_${ini}
		fi
		ini=`expr $ini + $step_npt`
	done
    execute "mv ${out_molec}_complex_nvt.gro ${folder_out_ucm}"
    execute "mv ${out_molec}_complex_nvt.cpt ${folder_out_ucm}"
    execute "mv ${out_molec}_complex_nvt.trr ${folder_out_ucm}"
	echo "Fin equilibrado"

}