#!/usr/bin/env bash
if [ "$#" -ne 5 ]; then
    echo "ERROR: The following parameters are needed:"
    echo "1º tpr simulation"
    echo "2º xtc simulation"
    echo "3º mode_gr  [	TARGET_QUERY | TARGET | TARGET_QUERIES | QUERIES | DNA_QUERY ] "
    echo "4º output xtc center"
    echo "5º prefix gromacs [ gmx | gmx_mpi ]"
    exit
fi
tpr=$1
xtc=$2
mode_gr=$3
out_xtc_center=$4
prefix_gromacs=$5

tmp_xtc="${xtc%.*}_complex_md_ori_tmp.xtc"
tmp_index="${xtc%.*}_complex_index.ndx.tmp"
tmp_gro="${xtc%.*}.gro"


#
#	Centrar la simulacion simulacion
#
echo \"Center simulation\"
if [ "${mode_gr}" == 'DNA_QUERY' ];then
    #echo "printf \"DNA\nSystem\n\" | ${prefix_gromacs} trjconv -s ${tpr} -f ${xtc} -center -ur compact -pbc mol -o ${tmp_xtc}"
    #exit
    printf "DNA\nSystem\n" | ${prefix_gromacs} trjconv -s ${tpr} -f ${xtc} -center -ur compact -pbc mol -o ${tmp_xtc}
    printf "DNA\nSystem\n" | ${prefix_gromacs} trjconv -s ${tpr}  -f ${tmp_xtc} -fit rot+trans -o ${out_xtc_center}
elif [ "${mode_gr}" == 'QUERIES' ] ;then #si es solo query e sobresscriben los parametros anteriores

    echo 'q' |${prefix_gromacs} make_ndx -f ${tmp_gro} > ${tmp_index}2>/dev/null
    aux="L01" # L0q siempre va a existir cuando se trate de queries ${name_queries[0]}
    aux=`cat ${tmp_index}  |grep ${aux} |awk '{print $1}' |head -1`
    printf "${aux}\nSystem\n" | ${prefix_gromacs} trjconv -s ${tpr} -f ${xtc} -center -ur compact -pbc mol -o ${tmp_xtc}
    printf "${aux}\nSystem\n" | ${prefix_gromacs} trjconv -s ${tpr} -f ${tmp_xtc} -fit rot+trans -o ${out_xtc_center}
    rm $tmp_gro $tmp_index
elif [ "${mode_gr}" == 'BIPHSIC_SYSTEMS' ];then
 printf "0\n0\n" | ${prefix_gromacs} trjconv -s ${tpr} -f ${xtc} -center -ur compact -pbc mol -o ${out_xtc_center}

else

    printf "Protein\nSystem\n" | ${prefix_gromacs} trjconv -s ${tpr} -f ${xtc} -center -ur compact -pbc mol -o ${tmp_xtc}
    printf "Backbone\nSystem\n" | ${prefix_gromacs} trjconv -s ${tpr} -f ${tmp_xtc} -fit rot+trans -o ${out_xtc_center}
fi
if [ -f ${tmp_xtc} ] ;then
    rm  ${tmp_xtc}
fi
































#
#	Una vez centrado el xtc se genera un pdb
#


#ejecutar()
#{
#    #
#    #	Centrar la simulacion simulacion
#    #
#    execute "echo \"Center simulation\""
#    if [ "${mode_gr}" == 'DNA_QUERY' ];then
#        execute "printf \"DNA\nSystem\n\" | ${prefix_gromacs} trjconv -s ${out_molec}_complex_md.tpr  -f ${out_molec}_complex_md.xtc -center -ur compact -pbc mol -o ${out_molec}_complex_md_ori.xtc"
#        execute "printf \"DNA\nSystem\n\" | ${prefix_gromacs} trjconv -s ${out_molec}_complex_md.tpr  -f ${out_molec}_complex_md_ori.xtc -fit rot+trans -o ${out_molec}_complex_md_center.xtc"
#    elif [ "${mode_gr}" == 'QUERIES' ] ;then #si es solo query e sobresscriben los parametros anteriores
#        echo 'q' |${prefix_gromacs} make_ndx${mpi} -f ${out_molec}_complex_md.gro > ${out_aux}_complex_index.ndx.tmp 2>/dev/null
#        aux=${name_queries[0]}
#        aux=`cat ${out_aux}_complex_index.ndx.tmp |grep ${aux} |awk '{print $1}' |head -1`
#        printf "${aux}\nSystem\n" | ${prefix_gromacs} trjconv -s ${out_molec}_complex_md.tpr -f ${out_molec}_complex_md.xtc -center -ur compact -pbc mol -o ${out_molec}_complex_md_ori.xtc
#        printf "${aux}\nSystem\n" | ${prefix_gromacs} trjconv -s ${out_molec}_complex_md.tpr  -f ${out_molec}_complex_md_ori.xtc -fit rot+trans -o ${out_molec}_complex_md_center.xtc
#
#    else
#
#        printf "Protein\nSystem\n" | ${prefix_gromacs} trjconv -s ${out_molec}_complex_md.tpr -f ${out_molec}_complex_md.xtc -center -ur compact -pbc mol -o ${out_molec}_complex_md_ori.xtc
#        printf "Backbone\nSystem\n" | ${prefix_gromacs} trjconv -s ${out_molec}_complex_md.tpr  -f ${out_molec}_complex_md_ori.xtc -fit rot+trans -o ${out_molec}_complex_md_center.xtc
#    fi

#    execute "rm  ${out_molec}_complex_md_ori.xtc"
    #
    #	Una vez centrado el xtc se genera un pdb
    #

#    dump=$(echo $step_md*0.002 | bc) #0.002 es el paso
    #execute "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${out_molec}_complex_md.xtc -s ${out_molec}_complex_md.gro -o ${out_molec}_complex_md_no_center.pdb -dump $dump"
#    execute "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${out_molec}_complex_md_center.xtc -s ${out_molec}_complex_md.gro -o ${out_molec}_complex_md.pdb -dump $dump"
#}
