#!/bin/bash
#	Script para ejecutar Gromacs automaticamente,
BLUE='\033[0;34m'
NONE='\033[00m' #no color
topologias="N/A"
path_scripts="${path_cluster_nodes}execute_scripts/scriptGR/"
MAX_WARNINGS=-1
export GMX_MAXCONSTRWARN=$MAX_WARNINGS
export GMX_MAXBACKUP=-1        #se eliminan backups de gromacs
#Fatal error:
#Too many LINCS warnings (5896)
#If you know what you are doing you can adjust the lincs warning threshold in your mdp file
#or set the environment variable GMX_MAXCONSTRWARN to -1,
#but normally it is better to fix the problem
#For more information and tips for troubleshooting, please check the GROMACS
#website at http://www.gromacs.org/Documentation/Errors


execute()
{
echo "--"
echo $1
    #2&> t.txt
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
    #	recojo parametros especiales, demmoemtno solo estos pocos
    #
    source ${path_scripts}read_params.sh
    source ${path_scripts}read_config.sh

    cp ${input_topology_data}.gro ${out_molec}_complex.gro
    cp ${input_topology_data}.top ${out_molec}_complex.top

    system_charge=`cat ${input_topology_data}.eng |grep Suma |awk '{print $2}'` #carga del sistema

    source ${path_scripts}search_itps.sh
    source ${path_scripts}read_name_queries.sh
    source ${path_scripts}reset_porses.sh
    #source ${path_scripts}version_gromacs.sh
    source ${path_scripts}resources.sh
    source ${path_scripts}generate_config_files.sh

    ejecutar

    source ${path_scripts}resume_input.sh


    #echo "bash ${path_scripts}center_simulation.sh ${out_molec}_complex_md.tpr ${out_molec}_complex_md.xtc ${mode_gr} ${out_molec}_complex_md_center.xtc ${prefix_gromacs}"
    #bash ${path_scripts}center_simulation.sh ${out_molec}_complex_md.tpr ${out_molec}_complex_md.xtc ${mode_gr} ${out_molec}_complex_md_center.xtc ${prefix_gromacs}

    #
    # 1
    #	Generamos
    #		Grid
    #		y rellenamos grid con disolvente
    source ${path_scripts}generate_grid.sh
    ejecutar


    #
    # 2
    #	Add ions para compensar la carga de la topologia
    #
    source ${path_scripts}generate_ions.sh
    ejecutar

    #   2b
    #   Se buscan los grupos de ligandos cuando se han agregado los iones y solven
    #
    source ${path_scripts}get_groups_number.sh
    #
    # 2
    #	Generate index
    #
    source ${path_scripts}generate_index.sh
    ejecutar


    #
    #4
    #	Minimizar energia para empezar con una primera structura y ahorrar tiempo en la simulacion
    #	debemos asegurarnos de que el sistema no tiene choques estéricos o la geometría inadecuada
    source ${path_scripts}minimization.sh
    ejecutar

    #
    # 5 Equilibracion 1 nvt 4npts cambiando las restricciones
    #
    last_equlibration="" #IMPORTANTE ultima equilibracion realizada
    source ${path_scripts}equilibration.sh
    ejecutar


    #6
    #	Simulacion
    #
    source ${path_scripts}simulation.sh
    ejecutar

    #
    #   Centrar simulacion y convertir a pdb
    #
    #execute "echo ${path_scripts}center_simulation.sh  ${out_molec}_complex_md.tpr ${out_molec}_complex_md.xtc ${mode_gr} ${out_molec}_complex_md_center.xtc ${prefix_gromacs}"
    bash ${path_scripts}center_simulation.sh ${out_molec}_complex_md.tpr ${out_molec}_complex_md.xtc ${mode_gr} ${out_molec}_complex_md_center.xtc ${prefix_gromacs}
    dump=-1
    bash ${path_scripts}create_pdb.sh ${out_molec}_complex_md_center.xtc ${out_molec}_complex_md.gro ${out_molec}_complex_md.pdb $dump $prefix_gromacs
    #dump=$(echo $step_md*0.002 | bc) #0.002 es el paso -1 igual last frame
    #if [ $dump -lt 1 ];then
    #	dump=$((dump - 1))
    #fi
    ## dump=-1
    ## execute "echo "bash ${path_scripts}create_pdb.sh ${out_molec}_complex_md_center.xtc ${out_molec}_complex_md.gro ${out_molec}_complex_md.pdb $dump $prefix_gromacs""
    #
    #   Analizar resultados
    #
    echo "${python_run}  ${path_analize_results}Simulation_gromacs/analyze_trajectory/get_results.py ${out_molec}_complex ${mode_gr} $prefix_gromacs"
    ${python_run}  ${path_analize_results}Simulation_gromacs/analyze_trajectory/get_results.py ${out_molec}_complex ${mode_gr} $prefix_gromacs

}