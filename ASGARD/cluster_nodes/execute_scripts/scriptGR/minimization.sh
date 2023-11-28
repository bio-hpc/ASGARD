#!/usr/bin/env bash
ejecutar()
{
	execute "echo \"Minimization ${step_min} steps\""
    minimizeEnergy
    execute "echo \"Fin Minimization\""
}
#
#
minimizeEnergy()
{	

	execute "${prefix_gromacs} grompp${mpi} -f ${out_grid}_min.mdp -c ${out_molec}_complex_solv_ions.gro -p ${out_molec}_complex.top -o ${out_molec}_complex_min.tpr -n  ${out_grid}_complex_index.ndx -maxwarn $MAX_WARNINGS"	#minimización de la energía a través de los GROMACS motor MD, mdrun
	execute "${prefix_gromacs} mdrun${mpi} -v -deffnm ${out_molec}_complex_min ${gpu}  ${thrads} -g ${out_aux}_mdrun_min.log"
														#-deffnm em es como poner em*
	execute "mv ${out_molec}_complex_min.edr ${folder_energies}"
	execute "mv ${out_molec}_complex_min.trr ${folder_out_ucm}"
	execute "mv ${out_molec}_complex_solv_ions.gro ${folder_out_ucm}"
}
