#!/usr/bin/env bash
#
# Create a custom configuration in ASGARD/config.conf
# This file contains:
#   ASGARD installation path
#   External SW path
#   ASGARD scripts/software directory
#   Python used
#   GROMACS modules and configuration

function write_comment()
{
	echo "#" >> ${file_conf}
	echo "# ${1}" >> ${file_conf}
	echo "#" >> ${file_conf}
}

file_conf=${path_ASGARD}/config.cfg

function get_python_run()
{	

	python_run=`which python`
	python_version=$(${python_run} --version 2>&1)
	python_version=`echo $python_version | cut -d\  -f 2`
	python_run="python"
	echo ${python_run} ${python_version}

}

if [ ! -f $file_conf ] ; then

	nom_serv=`echo $HOSTNAME`
	aux=`get_python_run`
	python_run=`echo $aux | cut -d\  -f 1`
	python_version=`echo $aux | cut -d\  -f 2`
	write_comment "Path"
	echo "path_ASGARD: ${path_ASGARD}" >> ${file_conf}
	echo "path_login_node: ${path_ASGARD}login_node/" >> ${file_conf}
	echo "path_analize_results: ${path_ASGARD}analyze_results/" >> ${file_conf}
	echo "path_external_sw: ${path_ASGARD}external_sw/" >> ${file_conf}
    echo "path_extra_ASGARD: $path_ASGARD}extra_ASGARD/" >> ${file_conf}
	write_comment "Python"
	echo "python_run: ${python_run}" >> ${file_conf}
	echo "python_version: ${python_version}" >> ${file_conf}

	write_comment "Molecular Dynamics"
	echo "g_mmpbsa: ${path_ASGARD}analyze_trajectory/extra/g_mmpbsa" >> ${file_conf}
	echo "g_gmx: gmx_mpi" >> ${file_conf}
	echo "g_amber_home: ${path_ASGARD}external_sw/amber14/" >> ${file_conf}



	write_comment "ASGARD Extra modules Escribir a continuacion si se necesita cargar algún modulo necesario apra shuttlemol"
	echo "#module load imagemagick/6.9.0-4" >> ${file_conf}
	write_comment "ASGARD Extra modules Escribir a continuacion si se necesita cargar algún modulo para get histogram"
	echo "#module laod openbabel/2.3.2" >> ${file_conf}
	write_comment "ASGARD end_modules"

	source ${path_ASGARD}login_node/read_file_conf.sh
	read_file ${path_ASGARD}	
else
	source ${path_ASGARD}login_node/read_file_conf.sh
	read_file ${path_ASGARD}	
fi




debugB "_________________________________________Global Config________________________________________"
debugB ""
debugB "Create_config.sh path_ASGARD: ${path_ASGARD}"
debugB "Create_config.sh path_login_node: ${path_login_node}"
debugB "Create_config.sh path_cluster_nodes: ${path_cluster_nodes}"
debugB "Create_config.sh path_analize_results: ${path_analize_results}"
debugB "Create_config.sh path_external_sw: ${path_external_sw}"
debugB "Create_config.sh path_extra_shuttlemol: ${extra_shuttlemol}"

debugB "Create_config.sh python_run: ${python_run}"
debugB "Create_config.sh python_version: ${python_version}"
debugB ""