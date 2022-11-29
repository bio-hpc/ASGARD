#!/bin/bash
# Inputs: 
# Carpeta de resultados de MD (xtc, tpr, gro,pdb_protein,mol2_ligand)

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
mkdir -p "$RESULTS"/{molecules,grids,energies,jobs,results} # o mejor mkdir -p "$RESULTS"-"$(date +%Y-%d-%m-%H:%M:%S)"
#cp $NAME/*.xtc $NAME/*.tpr $NAME/*.top  $RESULTS/molecules

for i in $(ls $ORIGIN/* | cut -d'.' -f2); do
      if [ $i  = 'top' ]; then
            cp $ORIGIN/*.$i $RESULTS/molecules/$NAME'.'$i                
      elif [ $i  = 'pdb' ]; then
            cp $ORIGIN/*.$i targets/
            PDB=$(ls $ORIGIN/*.pdb | sed s/^.*\\/\// | cut -d'.' -f1)
      elif [ $i  = 'mdp' ]; then
            cp $ORIGIN/*.$i $RESULTS/grids/$NAME'_md.'mdp
      elif [ $i  = 'mol2' ]; then
            cp $ORIGIN/*.mol2 queries/
            MOL2=$(ls $ORIGIN/*.mol2 | sed s/^.*\\/\// | cut -d'.' -f1)
      elif [ $i  = 'gro' ]; then
            for j in $(ls $ORIGIN/*$i); do
              if [[ $j = *"npt"* ]];then
                cp $j $RESULTS/molecules/$NAME'_npt_1.gro_md.'$i
              else
                cp $j $RESULTS/molecules/$NAME'_md.'$i
              fi
            done
      else
           cp $ORIGIN/*.$i $RESULTS/molecules/$NAME'_md.'$i
      fi
done

#cp $NAME/* $RESULTS/molecules
#cp $NAME/*.top $RESULTS/molecules/$NAME.top
mkdir targets/$NAME
cp $ORIGIN/*.pdb targets/$NAME
mkdir queries/$NAME
cp $ORIGIN/*.mol2 queries/$NAME

###########################
echo 'Centering trajectory...'
###########################

CENTER="$RESULTS"/molecules/"$NAME"_center.xtc

echo 1 0 | $gmx trjconv -s $RESULTS/molecules/*.tpr -f $RESULTS/molecules/*.xtc -center -ur compact -pbc mol -o "$RESULTS"/molecules/"$NAME"_ori_tmp.xtc
echo 4 0 | $gmx trjconv -s $RESULTS/molecules/*.tpr -f $RESULTS/molecules/"$NAME"_ori_tmp.xtc -fit rot+trans -o $CENTER

###########################
echo 'Generating pdb file...'
###########################

#echo 0 |$gmx trjconv -f $CENTER -s $RESULTS/molecules/*.tpr -o $RESULTS/molecules/"$NAME".pdb -tu ns -e 100 # ultimo frame

sh ASGARD/login_node/create_pdb.sh $CENTER $RESULTS/molecules/*.tpr $RESULTS/molecules/"$NAME".pdb -1 gmx


#############################
echo 'Generating topology'   
#############################

if [ -z "$(ls -A queries/$NAME)" ]; then
	singularity exec --bind $bind singularity/ASGARD.simg ASGARD/external_sw/gromacs/topology/generate_topology.py -t targets/$NAME/*.pdb -p TARGET
	cp targets/$NAME/*.top $RESULTS/molecules/$NAME.top
fi
singularity exec --bind $bind singularity/ASGARD.simg ASGARD/external_sw/gromacs/topology/generate_topology.py -t targets/$NAME/*.pdb -q queries/$NAME/
#
sh ASGARD/login_node/edit_topology.sh $RESULTS/molecules/$NAME queries/$NAME #prefijo
sh ASGARD/login_node/edit_include.sh queries/$NAME

##############################
echo 'Creating index files'
##############################

#echo 'q' | $gmx make_ndx -f $RESULTS/molecules/$NAME'_md.gro' -o $RESULTS/molecules/grids/VS_GR_4ejeB_preprocess_ebolaB_preprocess_complex_index.ndx

sh ASGARD/login_node/generate_index.sh $RESULTS/molecules/$NAME'_md.gro' $RESULTS/grids/"$NAME"_index.ndx queries/$NAME/*query.gro

############################## 
echo 'Creating resume'
##############################

TITLE=$PDB'_'$MOL2

sh ASGARD/login_node/create_resume.sh $(pwd)'/'$RESULTS $TITLE $RESULTS/grids/$NAME'_md.mdp' ${ORIGIN}/${PDB}.pdb ${ORIGIN}/${MOL2}.mol2




