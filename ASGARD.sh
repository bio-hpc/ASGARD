#!/bin/bash
error=0
paraMError=""
txtErrror=""
informe=""
fecha=$(date +"%Y-%m-%d")
pathSL="ASGARD/login_node/"
ASGARD_analysis="ASGARD/analyze_trajectory/get_results.py"
path_ASGARD="${PWD}/ASGARD/"
singularity="${PWD}/singularity/"
bind=$(pwd | cut -d/ -f1-2)/

source ${pathSL}colors.sh                                               #colores del script
source ${pathSL}parameters.sh                                           #parametros
if [[ ! -f ${path_ASGARD}/config.cfg ]]; then 
  source ${pathSL}create_conf.sh  2>/dev/null                           #crea oi carga configuracion global
  echo "Config file created"                                  
fi
                            
if [[ -z $folder_analysis ]]; then 
  echo $folder_analysis
  source ${pathSL}help.sh
  echo "Indicate absolute path where input folder is found"
  exit
fi

if [[ -z $profile ]]; then 
  echo $profile
  source ${pathSL}help.sh
  echo "You have to choose a valid profile (TARGET_QUERY or TARGET)"
  exit
fi

echo 'Folder: '$folder_analysis
echo 'Profile: '$profile

date

echo "Preparing analysis folder"

source ${pathSL}input_preparation.sh $folder_analysis    #crea una carpeta con todos los archivos necesarios

folder=$RESULTS

echo "Running analysis..."

prefix=$(ls $RESULTS"/molecules"|grep ".top"|cut -d. -f1-1)
echo $PWD/$RESULTS"/molecules/"$prefix

#singularity exec --bind $bind singularity/ASGARD.simg python $ASGARD_analysis $folder_analysis $profile gmx > "$folder_analysis"_"profile"_"$fecha".err 2>&1 # Análisis de MD
#python $ASGARD_analysis $PWD/$RESULTS"/molecules/"$prefix $profile gmx > "$folder_analysis"_"$profile"_"$fecha".err 2>&1 # Análisis de MD
singularity exec --bind $bind "$singularity"/ASGARD.simg python $ASGARD_analysis $PWD/$RESULTS"/molecules/"$prefix $profile gmx > "$folder_analysis"_"$profile"_"$fecha".err 2>&1 # Análisis de MD
rm mdout.mdp area.xvg
#singularity exec --bind $bind singularity/ASGARD.simg python $ASGARD_analysis"/molecules/(ls $ASGARD_analysis"/molecules/ -top)" $folder_analysis $profile gmx > "$folder_analysis"_"profile"_"$fecha".err 2>&1 # Análisis de MD
cd $RESULTS/results
latex=$(ls *documnet.tex) 
echo $latex
singularity exec --bind $bind "$singularity"/ASGARD.simg xelatex -interaction nonstopmode -file-line-error $latex
cp *.pdf ..

echo "Report "$latex" completed and generated in "$RESULTS"" 

date