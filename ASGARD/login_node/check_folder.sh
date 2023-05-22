#!/bin/bash
# Inputs: 
# MD Result folder (xtc, tpr, gro,pdb_protein,mol2_ligand)

DIR=$1
bind=$(pwd | cut -d/ -f1-2)/
singularity="${PWD}/singularity/"
gmx=$(echo singularity exec --bind $bind "$singularity"/ASGARD.simg gmx)


if [[ ${DIR: -1} == / ]]; then
DIR=$(echo $DIR | sed 's/.$//')
fi

ORIGIN=${DIR##*/}
NAME='VS_GR_'${DIR##*/}

RESULTS=$NAME'_results'-"$(date +%Y-%d-%m)"

if [[ -d $RESULTS ]]; then 
                while [ "$input" != "Y" ] && [ "$input" != "y" ] && [ "$input" != "N" ] && [ "$input" != "n" ] && [ "$input" != "zz" ] ; do
                        echo "Analysis folder already exists. Do you want to delete it?"
                        echo "(Y/y) Delete folder"
                        echo "(N/n) Exit"
                        read  input
                done
                if [ "$input" == "Y" ] || [ "$input" == "y" ];then
                        rm -r $RESULTS
                elif [ "$input" == "N" ] || [ "$input" == "n" ];then
                        exit
                fi  
fi
