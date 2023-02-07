#!/usr/bin/env bash
echo "echo \"start job\" 1>&2" 																			>>$name_template_job
echo "date +\"%s   %c\" 1>&2"																			>>$name_template_job
#echo "${path_cluster_nodes}techniques/baseTechniques.sh -c ${CWD} \
#-d ${folder_experiment} -s ${software} -p ${target} -l ${query} -x ${x} -y ${y} -z ${z} \
#-np ${name_target} -o ${option} -in ${confAdicionalI} -fn ${confAdicionalF} -cm ${comando} \
#-nc ${numPoses} -fl ${flexFile} -nj ${NombreJob} -fx ${flex} -na ${num_amino_acid} -ch ${chain} \
#-nl ${name_query} -de ${debug} -td ${time_experiment} -co ${cores} -gp ${GPU} -mm ${mem} -tr ${torsion} \
#-sc ${scoreCorte} -grid ${grid} -EXP ${ext_target} -EXL ${ext_query} \
#-GX ${gridSizeX} -GY ${gridSizeY} -GZ ${gridSizeZ} -nn  ${nodos} -lto ${lanzTimeOut} \
#-nj $jobName"																							>>$name_template_job

echo $modules_shuttlemol >>$name_template_job
echo "${path_cluster_nodes}techniques/baseTechniques.sh -c ${CWD} \
-d ${folder_experiment} -s ${software} -t ${target} -q ${query} -x ${x} -y ${y} -z ${z} \
-nt ${name_target} -o ${option} -in ${contIni} -fn ${contFin} \
-nc ${numPoses} -fl ${flexFile} -nj ${NombreJob} -fx ${flex} -na ${num_amino_acid} -ch ${chain} \
-nq ${name_query} -de ${debug} -td ${time_experiment} -co ${cores} -gp ${GPU} -mm ${mem} -tr ${torsion} \
-sc ${scoreCorte} -grid ${grid} -EXP ${ext_target} -EXL ${ext_query} \
-GX ${gridSizeX} -GY ${gridSizeY} -GZ ${gridSizeZ} -nn  ${nodos} -lto ${lanzTimeOut} \
-nj ${name_job} -BE ${bd_exhaustiveness} -BDA ${bd_atom_default}"		 														>>$name_template_job
echo "mv $name_template_job ${folder_jobs_done}"                                                        >>$name_template_job
echo "echo \"end job\" 1>&2" 																			>>$name_template_job
echo "date +\"%s   %c\" 1>&2"																			>>$name_template_job


