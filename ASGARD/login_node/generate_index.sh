#!/usr/bin/env bash

#
#   Se le pasa un nombre (Protein, RNA.. ) y devuelve el numero de grupo
#

get_num_group( ){

    echo `cat ${ndx}.tmp |grep -v "There" |grep -v "Analysing"|grep ${1} |awk '{print $1}' |head -1`
}

gro=$1
ndx=$2
query=$3
bind=$(pwd | cut -d/ -f1-2)/
singularity="${PWD}/singularity/"
gmx=$(echo singularity exec --bind $bind "$singularity"/ASGARD.simg)

echo 'q' | $gmx gmx make_ndx -f $gro -o $ndx > ${ndx}.tmp  2>/dev/null


#query=$(cat ../../queries/md1/4eje_preprocess_pduev3_pduev3_preprocess_complex_query.gro | head -3 | awk '{print $2}' | tail -1)
query=$(cat $query | head -3 | awk '{print $2}' | tail -1)


number_query=`get_num_group "$query"`
number_protein=`get_num_group "Protein"`
number_dna=`get_num_group "DNA"`
number_sol=`get_num_group "SOL"`

echo $number_sol
echo $number_query

number_all_groups=""

if [ "$number_protein" != "" ];then  number_all_groups=${number_protein} ;fi
if [ "$number_dna" != "" ];then number_all_groups="${number_all_groups} ${number_dna} |";fi
if [ "$number_sol" != "" ];then number_all_groups="${number_all_groups} ${number_sol}  |";fi
if [ "$number_query" != "" ];then number_all_groups="${number_all_groups} ${number_queries} |";fi
number_all_groups=${number_all_groups%.*}

i=`get_num_group "Water_and_ions"` #antes se hacia [ proteina | DNA | SOL Queries ] ahora solo water_and_ions_lig
number_all_groups=$i" | "$number_queries

ejecutar()
{
	echo 'Generate index\'
  generate_index
  echo 'Fin Generate index\'

}

generate_index()
{
    echo "${number_all_groups}" >  $ndx.tmp.tmp
    echo $ndx.tmp.tmp
    echo  $ndx.tmp.tmp
    echo "q" >>$ndx.tmp.tmp
    $gmx gmx make_ndx -f $gro -o $ndx < $ndx.tmp.tmp
#    rm  $ndx.tmp
#    rm  $ndx.tmp.tmp
}

ejecutar