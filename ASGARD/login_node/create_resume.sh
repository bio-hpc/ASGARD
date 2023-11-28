#!/bin/bash

generate_informe()
{   
	echo "___________________________________________Resume________________________________________________________">>${salidaResume}
	echo  "--">>${salidaResume}
	echo  "--">>${salidaResume}
	informe=${informe}"-- LigandIn:\t\t${query}\n"
	informe=${informe}"-- ReceptorIn:\t\t${target}\n"
	informe=${informe}"-- protocolProtIn:\t\t'GROMACS MD'\n"
	informe=${informe}"-- protocolLigIn:\t'GROMACS MD'\n"
	informe=${informe}"-- X:\t\t\t'0'\n"
	informe=${informe}"-- Y:\t\t\t'0'\n"
	informe=${informe}"-- Z:\t\t\t'0'\n"
	informe=${informe}"-- Grid:\t\t'N'\n"
	informe=${informe}"-- SizeGridX:\t\t'0'\n"
	informe=${informe}"-- SizeGridY:\t\t'0'\n"
	informe=${informe}"-- SizeGridZ:\t\t'0'\n"
	informe=${informe}"-- Software:\t\t'GR - Gromacs'\n"
	informe=${informe}"-- Technique:\t\t'GR'\n"
	informe=${informe}"-- Cores:\t\t'-'\n"
	informe=${informe}"-- Mem:\t\t\t'-'\n"
	informe=${informe}"-- NomJob:\t\t'-'\n"
	informe=${informe}"-- QueueManager:\t'-'\n"
	informe=${informe}"-- Directory:\t\t$folder_experiment\n"
	informe=${informe}"-- TotalRuns:\t\t'1'\n"
	informe=${informe}"-- run x Jobs:\t\t'1'\n"
	informe=${informe}"-- Jobs:\t\t'1'\n"
	informe=${informe}"-- Time Job:\t\t'-'\n"
	informe=${informe}"-- Time Dock:\t\t'-'\n"
	cl=`hostname`
	informe=${informe}"-- Cluster:\t\t$cl\n"
	informe=${informe}"--\n--\n"
	revision=`git rev-list --reverse HEAD 2>/dev/null |tail -1 `
	codVersion=`git rev-list --reverse HEAD 2>/dev/null |tail -1`
	numVersion=`git rev-list --reverse HEAD 2>/dev/null |wc -l`
	rama=`git branch -avv 2>/dev/null|grep "\*" |awk '{print $2}'`
	informe="${informe}-- numVesrion:\t $numVersion codVersion:\t\t$codVersion  \t branch: $rama\n" #fecha para el informe
	informe="${informe}-- date:\t\t$fecha (Y-m-d)\n" #fecha para el informe
    allComand=$(echo -e ./ASGARD.sh -t ./$target -q ./$query -o VS -s GR -j 1 -prp na -prl na -gp 1 -co 1 -td 259200 -el n'                '-step_min 10000 -write_data $(cat $mdp | grep nstlog | awk '{print $3}') -step_npt 200000 -step_nvt 200000 -step_md $(cat $mdp| grep nsteps | awk '{print $3}') -solvent tip3p -force_field amber99sb -solvatation SOL -temp $(cat $mdp | grep ref_t | awk '{print $3}') -bt dodecahedron -padding_grid 0.9 -seedg 2015 -prefix_gromacs gmx -pressure_npt 1.0)
	informe=${informe}"-- Command:\t\t$allComand $optAux"
	echo -e "${informe}" >>${salidaResume}
	echo  "--">>${salidaResume}
	echo  "--">>${salidaResume}
}

folder_experiment=$1
name=$2
mdp=$3
target=$4
query=$5


salidaResume=${folder_experiment}/Resume_VS_GR_1_"$name".txt
generate_informe
