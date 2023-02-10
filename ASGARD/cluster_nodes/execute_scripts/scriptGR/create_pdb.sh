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
#variable steps no se usa
if [ "$step" == "-1" ];then 
	##echo "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -e $step"
	##echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -e $step #Esto funciona solo con alguna version de gromacs
	echo "${prefix_gromacs} check  -f /home/jpg/workspace/shuttlemol/shuttlemol/VS_GR_1le0_test_dm_2018-12-18/molecules/VS_GR_1le0_test_dm_complex_md_center.xtc  2>&1 |grep Step |awk '{print $2}'"
	dump=`${prefix_gromacs} check  -f ${xtc}  2>&1 |grep Step |awk '{print $3}'`
	echo "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $dump"
	echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $dump
else
	echo "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $step"
	echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $step
fi





#dump=$(echo $step_md*0.002 | bc) #0.002 es el paso
#execute "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${out_molec}_complex_md.xtc -s ${out_molec}_complex_md.gro -o ${out_molec}_complex_md_no_center.pdb -dump $dump"
#execute "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${out_molec}_complex_md_center.xtc -s ${out_molec}_complex_md.gro -o ${out_molec}_complex_md.pdb -dump $dump"