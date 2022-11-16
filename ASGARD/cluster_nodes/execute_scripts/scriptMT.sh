#!/usr/bin/env bash
#
#
#   No FUNCIONA
#
execute_script()
{
  TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
  config_file_bd=${path_external_sw}metaDock/config-bd.met
  config_file_vs=${path_external_sw}metaDock/config-vs.met
  export LD_LIBRARY_PATH=/nfs/colab/bimbernon/Jorge/shuttlmol/shuttlemol/ShuttleMol/external_sw//metaDock/:$LD_LIBRARY_PATH
  config_file_bd=${path_external_sw}metaDock/config-bd.met
  config_file_vs=${path_external_sw}metaDock/config-vs.met
  output_dir=${folder_experiment}
  receptor=${CWD}${target}
  ligand=${CWD}${query}
  input_dir_vs=${CWD}${query}
  output_dir_vs=${folder_experiment}
  receptor_vs=${CWD}${target}
  pdbqt_dir=${path_external_sw}metaDock/pdbqt_dir
  param_field=${path_external_sw}metaDock/force_field.dat

  case "${option}" in
    VSB)
	#./energy -a $input_dir_vs $output_dir_vs 2> log

      n_rigids=`ls $input_dir_vs/rigid | wc -l`
      n_flex=`ls $input_dir_vs/flex | wc -l`

      if [ $n_rigids -gt 0 ]; then 
        input_rig=$input_dir_vs/rigid
        #./energy -c $2 -o $4 -i $3 -r $5 -d $input_rig -v 1 -f 0 2> log
        ${path_external_sw}/metaDock/energy -c $config_file_vs -o $output_dir_vs -p $receptor_vs -i $input_rig -q $pdbqt_dir -x $x -y $y -z $z -D $param_field -v 1 -f 0 2> ${out_aux}_R.log
      # #rm -f $input_flex
      fi
      if [ $n_flex -gt 0 ]; then
        input_flex=$input_dir_vs/flex
      #  # ./energy -c $2 -o $4 -i $3 -r $5 -d $input_flex -v 1 -f 1 2> log
        ${path_external_sw}/metaDock/energy -c $config_file_vs -o $output_dir_vs -p $receptor_vs -i $input_flex -q $pdbqt_dir -x $x -y $y -z $z -D $param_field -v 1 -f 1 2> ${out_aux}_F.log

      #   echo "Generating energy files"

      #  ./process_results/energy_report/energy_report -i $output_dir_vs -D ./process_results/ -m 1 2>log_er
        #rm -f $input_flex/*
      fi

      #echo "${path_external_sw}/metaDock/energy -c $config_file_bd -o $output_dir -p $receptor -l $ligand -D $param_field -q $pdbqt_dir 2> log"
      #${path_external_sw}/metaDock/energy -c $config_file -o $output_dir -p $receptor -l $ligand -D $param_field -q $pdbqt_dir 2> log
      # Generar resultados 
      #echo "Generating energy files"
      #./process_results/energy_report/energy_report -i $output_dir -l $ligand -D ./process_results/ -m 0
  ;;
    VS)

      ${path_external_sw}/metaDock/energy -c $config_file_bd -o $output_dir -p $receptor -l $ligand -D $param_field -q $pdbqt_dir 2> ${out_aux}.log
      #echo "${path_external_sw}/metaDock/energy -c $config_file_bd -o $output_dir -p $receptor -l $ligand -D $param_field -q $pdbqt_dir 2> ${out_aux}.log"
      ;;
    *) 
      echo "Invalid OPTION"
      echo "Usage BD: ./metadock.sh -BD config_file force_field_file pdbqt_dir output_dir receptor ligand"
            echo "Usage VS: ./metadock.sh -VS config_file force_field_file pdbqt_dir input_ligs output_dir receptor x_coord y_coord z_coord"
      ;;
  esac




































 
}






#case ${option} in
#VS)
  #echo $directorio
  #echo outTxt $outTxt
  #echo out_aux $out_aux
  #echo out_molec $out_molec
  #echo out_energies $out_energies
  #exit
  # echo `dirname /home/bimbernon/docking/VS_MT_2bsm_rec_testMT_GPU/txt/VS_MT_2bsm_rec_222bsm_lig.mol2`/1-`basename /home/bimbernon/docking/VS_MT_2bsm_rec_testMT_GPU/txt/VS_MT_2bsm_rec_222bsm_lig.mol2 ` 
  #  module load cuda
  #  module load openbabel
#  echo "Config File: $config_file_bd * Output directory: $directorio * Receptor: $target * Ligand: $query"  
  #echo "Config File: $2 * Input directory: $3 * Output directory: $4 * Receptor: $5 * Ligand: $6"
  #mkdir -p ${directorio}/best_ligs
  #mkdir -p ${directorio}/pymol
#  echo "Executing Metadock BD mode Now...OK"
#  execute "${pathSW}metaDock/energy -c $config_file_bd -o $directorio -p $target -l $query -D ${pathSW}metaDock 2> ${out_aux}.log"
  # Generar resultados 
#  echo "Generating results and optimization results"
#  execute "${pathSW}metaDock/energy_report -i $directorio -l $query -D ${pathSW}metaDock/ -m 0 "
#  ;;

#VSB)


#  echo "Config File: $config_file_vs * input_directory: $query * output directory: $directorio * Receptor: $target "
#        mkdir -p $query/rigid
#        mkdir -p $query/flex

#  echo "Executing Metadock VS mode Now...OK"
#  
#        execute "${pathSW}metaDock/energy -a $query $directorio 2> ${out_aux}.log"#


#        n_rigids=`ls $query/rigid | wc -l`
#        n_flex=`ls $query/flex | wc -l`
#
#         if [ $n_rigids -gt 0 ]; then
#                input_rig=$query/rigid
#                # ./energy -c $2 -o $4 -i $3 -r $5 -d $input_flex -v 1 -f 1 2> log
#                execute "${pathSW}metaDock/energy -c $config_file_vs -o $directorio -p $target -i $input_flex -x $x -y $y -z $z -D ${pathSW}metaDock/ -v 1 -f 0 2> ${out_aux}.log"
#        fi
#        if [ $n_flex -gt 0 ]; then
#               input_flex=$query/flex
#                # ./energy -c $2 -o $4 -i $3 -r $5 -d $input_flex -v 1 -f 1 2> log
#                execute "${pathSW}metaDock/energy -c $config_file_vs -o $directorio -p $target -i $input_flex -x $x -y $y -z $z -D ${pathSW}metaDock/ -v 1 -f 1 2> ${out_aux}.log"
#        fi
#  echo "Generating results and optimization results"
#        execute "${pathSW}metaDock/energy_report -i $directorio -D ${pathSW}metaDock/ -m 1 "
#  ;;
#*) 
#  echo "Invalid OPTION"
#  ;;
#esac
  
