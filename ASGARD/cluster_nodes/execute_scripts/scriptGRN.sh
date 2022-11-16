#!/bin/bash
#
#	Script para Continuar Ejecucion de Gromacs Automaticamente
#	informacion sacada de:
#	https://cluster.earlham.edu/wiki/index.php/Checkpoint_and_Restarting
#
BLUE='\033[0;34m'
NONE='\033[00m' #no color
topologias="N/A"

#Fatal error:
#Too many LINCS warnings (5896)
#If you know what you are doing you can adjust the lincs warning threshold in your mdp file
#or set the environment variable GMX_MAXCONSTRWARN to -1,
#but normally it is better to fix the problem
#For more information and tips for troubleshooting, please check the GROMACS
#website at http://www.gromacs.org/Documentation/Errors
export GMX_MAXCONSTRWARN=-1
FACTOR_DT=0.002
execute()
{

	echo -e "${BLUE}$1${NONE}"
	eval "$1"
	error=$?
	echo -e "${BLUE}ERROR: $error ${NONE}"
	if [ "$error" != "0" ];then
	 		exit;
	fi
}
execute_script()
{
    #
    #   Leemos parametros
    #
    path_scripts="${path_cluster_nodes}execute_scripts/scriptGR/"
    prefix="_complex"
    source ${path_scripts}read_params.sh
    source ${path_scripts}read_config.sh
    #
    #   cambiamos los pasos a ps
    #
    step_md=`echo "${step_md}*${FACTOR_DT}" | bc`
    #
    #   Cambiamos el prefix GRN por GR
    #
    out_prefix=${out_prefix%.*}
    out_prefix=`echo "${out_prefix}" | sed -e "s,GRN,GR,g"`
    #
    #   Asignamos nuevas varaiebls con el prefix modificado
    #
    out_grid=${folder_grid}${out_prefix}
    out_ucm=${folder_out_ucm}${out_prefix}${prefix}
    out_molec=${folder_molec}${out_prefix}${prefix}
    echo ${out_ucm}
    #
    #   movemos el tpr antiguo para crear uno nuevo
    #
    new_tpr=${out_molec}_md.tpr
    old_tpr=${out_molec}_md_old.tpr

    cp $new_tpr $old_tpr #guardo el tpr antiguo # cp es más seguro
    #
    #   Extendemos simulación
    #
    execute "${prefix_gromacs} convert-tpr -s ${old_tpr} -extend ${step_md} -o ${new_tpr}  ## -e ${out_molec}_md.edr -f ${out_molec}_md.trr"
    #
    #   Lanzamos simulacion
    #
    execute "${prefix_gromacs} mdrun -s ${new_tpr} -cpi ${out_molec}_md.cpt ${gpu} ${thrads} -append -deffnm ${out_molec}_md -g ${out_ucm}_md.log"
    #
    #   Centramos simulacion
    #
    bash ${path_scripts}center_simulation.sh ${out_molec}_md.tpr ${out_molec}_md.xtc ${mode_gr} ${out_molec}_md_center.xtc ${prefix_gromacs}
    dump=-1
    bash ${path_scripts}create_pdb.sh ${out_molec}_md_center.xtc ${out_molec}_md.gro ${out_molec}_md.pdb $dump $prefix_gromacs
    #
    #   Analizar resultados
    #
    echo "${python_run}  ${path_analize_results}Simulation_gromacs/analyze_trajectory/get_results.py ${out_molec} ${mode_gr} $prefix_gromacs"
    ${python_run}  ${path_analize_results}Simulation_gromacs/analyze_trajectory/get_results.py ${out_molec} ${mode_gr} $prefix_gromacs


}





    #echo $out_molec
    ##echo $prefix_gromacs
    ##echo $step_md
    ##exit
    ##querySnExt=` echo ${query%.*}`
    #
    #	recojo parametros especiales, demmoemtno solo estos pocos
    #
    ##patIn=$(dirname "$query")"/"



    ##salida=${salida%.*}
    ##salida=`echo "${salida}" | sed -e "s,GRN,GR,g"`

    ##out_grid=${folderGrid}${salida}
    ##outTxt=${folderTxt}${salida}
    ##out_aux=${folderOutUcm}${salida}
    #out_molec=${folderMolec}${salida}
    #out_energies=${folderEnergy}${salida}
    ##newTpr=${out_molec}_mdrum_simulation.tpr
    ##oldTpr=${out_molec}_mdrum_simulation_old.tpr
    ##mv $newTpr $oldTpr #guardo el tpr antiguo
    #
    #	Extiendo el nuevo tpr a lo que se indique
    #

    ## execute "$gromacs5 convert-tpr${mpi} -s ${oldTpr} -extend $stepSimulacionDM -o ${newTpr} -e ${out_molec}_mdrum_simulation.edr -f ${out_molec}_mdrum_simulation.trr"


    ## execute "$gromacs5 mdrun${mpi} -s ${out_molec}_mdrum_simulation.tpr -cpi ${out_molec}_mdrum_simulation.cpt ${gpu} ${thrads} -append -deffnm ${out_molec}_mdrum_simulation -g ${out_aux}_mdrum_simulation.log"
    #
    #	Centrar simulacion y genera el pdb
    #
    ##printf "Protein\nSystem\n" | gmx_mpi trjconv -s ${out_molec}_mdrum_simulation.tpr -f ${out_molec}_mdrum_simulation.xtc -center -ur compact -pbc mol -o ${out_molec}_mdrum_simulationOriginal_1.xtc
    ##printf "Backbone\nSystem\n" | gmx_mpi trjconv -s ${out_molec}_mdrum_simulation.tpr  -f ${out_molec}_mdrum_simulationOriginal_1.xtc -fit rot+trans -o ${out_molec}_mdrum_simulation_center.xtc
    ##execute "rm  ${out_molec}_mdrum_simulationOriginal_1.xtc"
    ##execute "$gromacs5 editconf${mpi} -f ${out_molec}_mdrum_simulation.gro -o ${out_molec}_mdrum_simulation.pdb" #no hace falta pero convierto la simulacion a pdb

    ##echo "El modo es: "$mode
    ##if [ "${mode}" != "" ];then
    ##    echo  "python ${pathSW}gromacs/analizarResults/getResults.py ${out_molec} ${mode}  $versionG"
    ##    python ${pathSW}gromacs/analizarResults/getResults.py ${out_molec} ${mode} $versionG
    ##fi