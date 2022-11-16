#!/usr/bin/env bash
#
#   Lanza la simulacion de DM el timepo que se le indique
#

ejecutar()
{
	execute "echo \"generateSimulacion\""
	simulacion
    execute "echo \"End Simulation\""
}
#
#
simulacion()
{
	execute "${prefix_gromacs} grompp${mpi} -f ${out_grid}_md.mdp -c ${last_equlibration}.gro -t ${last_equlibration}.cpt -p ${out_molec}_complex.top -o ${out_molec}_complex_md.tpr  -n ${out_grid}_complex_index.ndx  -maxwarn ${MAX_WARNINGS}"
	execute "${prefix_gromacs} mdrun${mpi} -deffnm ${out_molec}_complex_md ${gpu} ${thrads} -g ${out_aux}_complex_md.log"
	#execute "mv ${last_equlibration}.gro ${folder_out_ucm}"
    #execute "mv ${last_equlibration}.cpt ${folder_out_ucm}"
}

