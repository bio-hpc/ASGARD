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

source ${pathSL}colors.sh                                               #script colors
source ${pathSL}parameters.sh                                           #parameters
if [[ ! -f ${path_ASGARD}/config.cfg ]]; then 
  source ${pathSL}create_conf.sh  2>/dev/null                           # create the global config (config.cfg) 
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

echo "Preparing analysis folder..."

source ${pathSL}check_folder.sh $folder_analysis # check if the folder is already exists
echo "Generating the working folder..."
echo "Input preparation" > "$folder_analysis"_"$profile"_"$fecha".err
singularity exec --bind $bind "$singularity"/ASGARD.simg ./${pathSL}input_preparation.sh $folder_analysis >> "$folder_analysis"_"$profile"_"$fecha".err 2>&1  #create a folder with all the required files

if [ "$input" == "n" ] || [ "$input" == "N" ];then
  exit
fi

folder=$RESULTS

echo "Running analysis..."

prefix=$(ls $RESULTS"/molecules"|grep ".top"|cut -d. -f1-1)
echo "ASGARD analysis" >> "$folder_analysis"_"$profile"_"$fecha".err

singularity exec --bind $bind "$singularity"/ASGARD.simg python $ASGARD_analysis $PWD/$RESULTS"/molecules/"$prefix $profile gmx >> "$folder_analysis"_"$profile"_"$fecha".err 2>&1 # MD Analysis

echo "Generating PDF report"
rm mdout.mdp area.xvg
cd $RESULTS/results
latex=$(ls *documnet.tex) 
singularity exec --bind $bind "$singularity"/ASGARD.simg xelatex -interaction nonstopmode -file-line-error $latex >> "$folder_analysis"_"$profile"_"$fecha".err 2>&1 # Report
cp *.pdf ..

echo "Report "$latex" completed and generated in "$RESULTS"!" 

date
