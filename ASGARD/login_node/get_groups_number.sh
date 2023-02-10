#!/usr/bin/env bash

#
#   Return group number
#
get_num_group( ){

    echo `cat ${out_aux}_complex_index.ndx.tmp |grep -v "There" |grep -v "Analysing"|grep ${1} |awk '{print $1}' |head -1`
}



echo 'q' |${prefix_gromacs} make_ndx${mpi}   -f ${out_molec}_complex_solv_ions.gro  -o ${out_aux}_complex_index.ndx > ${out_aux}_complex_index.ndx.tmp  2>/dev/null


number_queries=""
for query in "${name_queries[@]}"; do
    aux=`get_num_group $query`

    number_queries="$number_queries  ${aux} |"
done
number_queries=${number_queries::-1}
number_protein=`get_num_group "Protein"`
number_dna=`get_num_group "DNA"`
number_sol=`get_num_group "SOL"`

number_all_groups=""


if [ "$number_protein" != "" ];then  number_all_groups=${number_protein} ;fi
if [ "$number_dna" != "" ];then number_all_groups="${number_all_groups} ${number_dna} |";fi
if [ "$number_sol" != "" ];then number_all_groups="${number_all_groups} ${number_sol}  |";fi
if [ "$number_queries" != "" ];then number_all_groups="${number_all_groups} ${number_queries} |";fi
number_all_groups=${number_all_groups::-1}

i=`get_num_group "Water_and_ions"` 
number_all_groups=$i" | "$number_queries
