#!/bin/bash
# Inputs: 
# MD Result folder (xtc, tpr, gro,pdb_protein,mol2_ligand)

DIR=$1
LIGAND=$2
bind=$(pwd | cut -d/ -f1-2)/
singularity="${PWD}/singularity/"
gmx=$(echo singularity exec --bind $bind "$singularity"/ASGARD.simg gmx)
echo $gmx

############################
echo 'Creating analysis folder...'
###########################

if [[ ${DIR: -1} == / ]]; then
DIR=$(echo $DIR | sed 's/.$//')
fi

ORIGIN=${DIR##*/}
NAME='VS_GR_'${DIR##*/}

RESULTS=$NAME'_results'-"$(date +%Y-%d-%m)"
echo $RESULTS

mkdir -p "$RESULTS"/{molecules,grids,energies,jobs,results} 

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
        if ls $ORIGIN/*.$i 1> /dev/null 2>&1; then
          cp $ORIGIN/*.$i $RESULTS/molecules/
        fi
      fi
done

if [ -d $ORIGIN/*'.ff' ]; then
    cp -r $ORIGIN/*.'ff' $RESULTS/molecules/
fi

###########################
echo 'Centering trajectory...'
###########################

CENTER="$RESULTS"/molecules/"$NAME"_center.xtc

echo 1 0 | $gmx trjconv -s $RESULTS/molecules/*.tpr -f $RESULTS/molecules/*.xtc -center -ur compact -pbc mol -o "$RESULTS"/molecules/"$NAME"_ori_tmp.xtc
echo 4 0 | $gmx trjconv -s $RESULTS/molecules/*.tpr -f $RESULTS/molecules/"$NAME"_ori_tmp.xtc -fit rot+trans -o $CENTER

###########################
echo 'Generating pdb file...'
###########################
 
sh ASGARD/login_node/create_pdb.sh $CENTER $RESULTS/molecules/*.tpr $RESULTS/molecules/"$NAME".pdb '-1' gmx

###########################
echo 'Checking topology and forcefield...'
###########################

FILES=""
TOP=$RESULTS/molecules/$NAME'.'top
for i in $(cat $TOP | grep '#include' | cut -d'"' -f 2); do
  if [ -f "$i" ]; then 
      echo "File found"
      for x in $(cat $i | grep '#include' | cut -d'"' -f 2); do
          if [ -f "$x" ]; then
            echo "File found"
          elif  [ -f $ORIGIN/$(basename $x) ]; then
            REPLACE=$(realpath $ORIGIN/$(basename $x))
            sed -i "s|$x|$REPLACE|g" $ORIGIN/$(basename $i)
          else
            FILES="${FILES} $x"
        fi
      done
  elif [ -f $ORIGIN/$(basename $i) ]; then
      REPLACE=$(realpath $ORIGIN/$(basename $i))
      sed -i "s|$i|$REPLACE|g" $TOP
      for j in $(cat $REPLACE | grep '#include' | cut -d'"' -f 2); do
        if [ -f $ORIGIN/$(basename $j) ]; then
          REPLACE=$(realpath $ORIGIN/$(basename $j))
          sed -i "s|$j|$REPLACE|g" $ORIGIN/$(basename $i)
        else 
          FILES="${FILES} $j"
        fi
      done
  else 
      FILES="${FILES} $i" 
  fi
done

if [ ! -z "$FILES" ]; then
  echo "Please include the following forcefield and/or topology files in the input folder:"
  for j in $FILES; do
    echo $(basename $j)
  done
fi

if [ ! -z "$FILES" ]; then
    exit
fi



##############################
echo 'Creating index files'
##############################

sh ASGARD/login_node/generate_index.sh $RESULTS/molecules/$NAME'_md.gro' $RESULTS/grids/"$NAME"_index.ndx "$LIGAND"

############################## 
echo 'Creating resume'
##############################

TITLE=$PDB'_'$MOL2

sh ASGARD/login_node/create_resume.sh $(pwd)'/'$RESULTS $TITLE $RESULTS/grids/$NAME'_md.mdp' ${ORIGIN}/${PDB}.pdb ${ORIGIN}/${MOL2}.mol2




