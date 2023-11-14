#!/bin/bash
# Inputs: 
# MD Result folder (xtc, tpr, gro,pdb_protein,mol2_ligand)

DIR=$1
bind=$(pwd | cut -d/ -f1-2)/
singularity="${PWD}/singularity/"
gmx=$(echo singularity exec --bind $bind "$singularity"/ASGARD.simg gmx)
echo $gmx
#
#PDB=$2
#MOL2=$3
#
#TRAJ=$1
#TOP=$2 # tpr
#GRO=$3
#PDB=$4
#MOL2=$5
#NAME=$6 # optional
#
#
#
##PARA ONLY TARGET
#
## if NAME is empty
#

############################
echo 'Creating analysis folder...'
###########################

if [[ ${DIR: -1} == / ]]; then
DIR=$(echo $DIR | sed 's/.$//')
fi

ORIGIN=${DIR##*/}
NAME='VS_GR_'${DIR##*/}

#NAME="${f%.*}"
#  NAME=$(echo $TRAJ | cut -f 1 -d '/' | cut -f 1 -d '.')
#fi
RESULTS=$NAME'_results'-"$(date +%Y-%d-%m)"
echo $RESULTS

if [[ -d $RESULTS ]]; then 
                while [ "$input" != "Y" ] && [ "$input" != "y" ] && [ "$input" != "N" ] && [ "$input" != "n" ] && [ "$input" != "zz" ] ; do
                        echo "Analysis folder already exists. Do you want to delete it?"
                        echo "(Y/y) Delete folder"
                        echo "(N/n) Exit"
                        read  input
                done
                if [ "$input" == "Y" ] || [ "$input" == "y" ];then
                        rm -r $RESULTS
                elif [ "$input" == "n" ] || [ "$input" == "N" ];then
                        exit
                fi  
fi
