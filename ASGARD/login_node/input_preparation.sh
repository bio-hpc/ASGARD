#!/bin/bash
# Inputs: 
# MD Result folder (xtc, tpr, gro,pdb_protein,mol2_ligand)

DIR=$1
bind=$(pwd | cut -d/ -f1-2)/
singularity="${PWD}/singularity/"
gmx=$(echo singularity exec --bind $bind "$singularity"/ASGARD.simg gmx)

############################
echo 'Creating analysis folder...'
###########################

if [[ ${DIR: -1} == / ]]; then
DIR=$(echo $DIR | sed 's/.$//')
fi

ORIGIN=${DIR##*/}
NAME='VS_GR_'${DIR##*/}
RESULTS=$NAME'_results'-"$(date +%Y-%d-%m)"

mkdir -p "$RESULTS"/{molecules,grids,energies,jobs,results} # better mkdir -p "$RESULTS"-"$(date +%Y-%d-%m-%H:%M:%S)"

for i in $(ls $ORIGIN/* | cut -d'.' -f2); do
      if [ $i  = 'top' ]; then
            cp $ORIGIN/*.$i $RESULTS/molecules/$NAME'.'$i
      elif [ $i  = 'tpr' ]; then
            cp $ORIGIN/*.$i $RESULTS/molecules/$NAME'_md.'$i
      elif [ $i  = 'xtc' ]; then
            cp $ORIGIN/*.$i $RESULTS/molecules/$NAME'_md.'$i                     
      elif [ $i =  'edr' ]; then
            cp $ORIGIN/*.$i $RESULTS/molecules/$NAME'_md.'$i
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
        cp $ORIGIN/*.$i $RESULTS/molecules/
      fi
done

if [ -d $ORIGIN/*'.ff' ]; then
    cp -r $ORIGIN/*.'ff' $RESULTS/molecules/
fi
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

sh ASGARD/login_node/create_pdb.sh $CENTER $RESULTS/molecules/*.tpr $RESULTS/molecules/"$NAME".pdb -1 gmx

#############################
echo 'Generating topology'   
#############################

if [ -z "$(ls -A queries/$NAME)" ]; then
	singularity exec --bind $bind singularity/ASGARD.simg ASGARD/external_sw/gromacs/topology/generate_topology.py -t targets/$NAME/*.pdb -p TARGET
fi

##############################
echo 'Creating index files'
##############################

sh ASGARD/login_node/generate_index.sh $RESULTS/molecules/$NAME'_md.gro' $RESULTS/grids/"$NAME"_index.ndx queries/$NAME/*query.gro

##############################
echo 'Creating tpr file for g_mmpbsa'
##############################

/opt/gromacs/bin/gmx grompp -f $RESULTS/grids/$NAME'_md.mdp' -c $RESULTS/molecules/$NAME'_md.'gro -p $RESULTS/molecules/$NAME'.top' -o $RESULTS/molecules/$NAME'_pre_md.tpr' -n $RESULTS/grids/$NAME'_index.ndx' -maxwarn 5

############################## 
echo 'Creating resume'
##############################

TITLE=$PDB'_'$MOL2

sh ASGARD/login_node/create_resume.sh $(pwd)'/'$RESULTS $TITLE $RESULTS/grids/$NAME'_md.mdp' ${ORIGIN}/${PDB}.pdb ${ORIGIN}/${MOL2}.mol2




