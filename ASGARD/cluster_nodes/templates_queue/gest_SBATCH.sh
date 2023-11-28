#!/usr/bin/env bash
echo "#!/bin/sh" 										> $name_template_job
if [ "$renice" != "N/A" ];then
	echo "#SBATCH --nice=$renice"						>>$name_template_job 
fi
if [ "${GPU}" != "N/A" ];then
	echo "#SBATCH --gres=gpu:${GPU}" 					>>$name_template_job 
fi
if [ "${mem}" != "0" ];then 
	echo "#SBATCH --mem=${mem}" 						>>$name_template_job
fi
if [ "$cluster_queue" != "N/A" ];then
	echo "#SBATCH -p "${cluster_queue} 							>>$name_template_job
fi

echo "#SBATCH --output=${folder_out_jobs}${outJob}.out"	>>$name_template_job
echo "#SBATCH --error=${folder_out_jobs}${outJob}.err"	>>$name_template_job
echo "#SBATCH -J ${name_job}"							>>$name_template_job
echo "#SBATCH --time=$time_job"							>>$name_template_job
echo "#SBATCH --cpus=$cores"							>>$name_template_job
echo "#SBATCH --nodes=${nodos}"							>>$name_template_job
source ${path_cluster_nodes}templates_queue/codigo.sh















































#sed "-e s,DDDDD,${directorio},g" $CWD""scriptsDocking/cabecerasCola/sbatchMalga.sh  \
#| sed -e "s,SSSSS,${software},g" \
#| sed -e "s,PPPPP,${target},g"\
#| sed -e "s,LLLLL,${query},g" \
#| sed -e "s,XXXXX,${x},g" \
#| sed -e "s,YYYYY,${y},g"\
#| sed -e "s,ZZZZZ,${z},g"  \
#| sed -e "s,IIIII,${confAdicionalI},g" \
#| sed -e "s,FFFFF,${confAdicionalF},g"\
#| sed -e "s,NNNNN,${name_target},g"\
#| sed -e "s,KKKKK,${NombreJob},g" \
#| sed -e "s,CCCCC,${comando},g" \
#| sed -e "s,WWWWW,${CWD},g" \
#| sed -e "s,UUUUU,${ultimo},g" \
#| sed -e "s,GGGGG,${num_amino_acid},g" \
#| sed -e "s,TTTTT,${flex},g" \
#| sed -e "s,CHCHCH,${chain},g" \
#| sed -e "s,CPCPCP,${cpus},g" \
#| sed -e "s,OOOOO,${option},g"	> $nomTemplateJob	
#execute="sbatch"
