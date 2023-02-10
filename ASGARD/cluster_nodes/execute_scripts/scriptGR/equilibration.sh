#!/usr/bin/env bash
#
#	El fichero _posre.itp esta inciiado a 1000 por defecto, según los conocimientos de Hugo lo ideal es que empieze con 5000 en equilibracionA
#	Seguidamente 4000,3000,2000,1000 en equilibracionB
#

ejecutar()
{
    equilibration_nvt
    equilibration_npt
}
#NVT
equilibration_nvt() #calentar la conformacion
{
	#	Reiniciamos los porses por si acaso, deberian esta reiniciados
	change_porse_files 5000 5000 5000 5000

    #se añade opcion -r por gromacs 2018
	execute "${prefix_gromacs} grompp${mpi} -f ${out_grid}_nvt.mdp -c ${out_molec}_complex_min.gro -r ${out_molec}_complex_min.gro -p ${out_molec}_complex.top -o ${out_molec}_complex_nvt.tpr -n  ${out_grid}_complex_index.ndx -maxwarn $MAX_WARNINGS"

	execute "${prefix_gromacs} mdrun${mpi} -deffnm ${out_molec}_complex_nvt ${gpu} ${thrads} -g ${out_aux}_complex_mdrun_nvt.log"
	execute "mv ${out_molec}_complex_nvt.edr ${folder_energies}"
	execute "mv ${out_molec}_complex_nvt.tpr ${folder_out_ucm}"

	execute "mv ${out_molec}_complex_min.gro ${folder_out_ucm}"



}
#NPT
equilibration_npt() #se ejecuta 4 veces incrementando la posicion de restriccion
{
	#
	#	Se hace 4 veces
	#
	ini=${step_nvt}
	posre_fin=4000
	for ((j=1; j<5; j++))
	#for j in `seq 1 4`
	do
	    change_porse_files ${posre_fin} ${posre_fin} ${posre_fin}
		#echo "EL NUMERO  DE ORSE INI  ES: "${posreIni}
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

	#
	#Moviendo ficheros no necesarios
	#
    execute "mv ${out_molec}_complex_nvt.gro ${folder_out_ucm}"
    execute "mv ${out_molec}_complex_nvt.cpt ${folder_out_ucm}"
    execute "mv ${out_molec}_complex_nvt.trr ${folder_out_ucm}"
	echo "Fin equilibrado"

}