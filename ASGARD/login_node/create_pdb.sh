#!/usr/bin/env bash
if [ "$#" -ne 5 ]; then
    echo "ERROR: The following parameters are needed:"
    echo "1º xtc"
    echo "2º gro"
    echo "3º out_pdb"
    echo "4 num step"
    echo "5º prefix_gromacs"
    exit
fi

xtc=$1
gro=$2
out_pdb=$3
step=$4
prefix_gromacs=$5
bind=$(pwd | cut -d/ -f1-2)/
singularity="${PWD}/singularity/"
gmx=$(echo singularity exec --bind $bind "$singularity"/ASGARD.simg)
#echo $gmx
if [ "$step" == "-1" ];then 
	##echo "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -e $step"
	##echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -e $step #Just for some gromacs versions
	echo "${gmx} ${prefix_gromacs} check  -f /home/jpg/workspace/shuttlemol/shuttlemol/VS_GR_1le0_test_dm_2018-12-18/molecules/VS_GR_1le0_test_dm_complex_md_center.xtc  2>&1 |grep Step |awk '{print $2}'"
	dump=`${gmx} ${prefix_gromacs} check  -f ${xtc}  2>&1 |grep Step |awk '{print $3}'`
	echo "echo 0 |${gmx} ${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $dump"
	echo 0 |${gmx} ${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $dump
else
	echo "echo 0 |${gmx} ${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $step"
	echo 0 |${gmx} ${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $step
fi


