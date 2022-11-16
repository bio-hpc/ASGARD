#!/usr/bin/env bash
execute_script()
{	
    
	unset DISPLAY
	TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
	lib=`echo $(basename $query)`
	if [ ${option} == "VSB" ];then
		opt_aux=${opt_aux/-TT/-T} #como el parametreo -T esta asignado a target, se hace un mapeo a -TT y aqui se vuelve a mapear a -t
		opt_aux=${opt_aux/-sf/-S} #como el parametreo -S esta asignado a software, se hace un mapeo a -sf y aqui se vuelve a mapear a -s
		execute "${path_external_sw}ligandScout/iscreen -q ${CWD}${target} -d ${CWD}${query} -M ${mem} -C ${cores} -F ${ini} -L ${fin} -o ${out_molec}.sdf ${opt_aux} &>${out_aux}.ucm"
		#
		#   Genera ficheros de caracteristicas y los json apra listas cruzadas
		#
		if [ -f "${out_molec}.sdf" ] && [ -s "${out_molec}.sdf" ]; then

			execute "${path_external_sw}tools/babel -m ${out_molec}.sdf ${out_molec}_.sdf &> /dev/null" # &>salida.ucam
			#execute "${path_external_sw}tools/babel -m ${out_molec}.sdf ${out_molec}_.sdf" # &>salida.ucam
			for fich in ${out_molec}_*.sdf; do				
				name=`head -n 1 "${fich}"`						
				name="$(echo -e "${name}" | tr -d '[:space:]')"				
		    	execute "filename=`echo $(basename $fich)`"
		    	filename=${filename%.*}
		    	filename=${filename%_*}_${name}
		    	out_energies=${folder_energies}/${filename} #ojo se pisa la salida de energias por defecto no eliminar
				# extract info about matching features

				execute "H=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"H:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				execute "HBD=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"HBD:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				execute "HBA=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"HBA:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				execute "AR=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"AR:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				execute "PI=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"PI:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				execute "NI=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"NI:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				execute "ZNB=`grep \"<Matching Features>\" -A1 $fich | tail -1 | grep -o \"ZNB:[0-9]*\" | uniq -c | awk -F: '{ feature=$2"_"feature} END {print feature}'`"
				#total number of matches
				#execute "totalFeatMatch=$((H + HBD + HBA + AR + PI + NI + ZNB))"
				#extract scoring value
				execute "score=`grep \"<Score\" -A1 $fich | tail -1`"				
				# write files
				#execute "echo \"${name} ${score} ${totalFeatMatch} "H:"$H "HBD:"$HBD "HBA:"$HBA "AR:"$AR "PI:"$PI "NI:"$NI "ZNB:"$ZNB\" > ${folder_energies}/${filename}.feat"
				execute "echo \"${name} ${score} "H:"$H "HBD:"$HBD "HBA:"$HBA "AR:"$AR "PI:"$PI "NI:"$NI "ZNB:"$ZNB\" > ${folder_energies}/${filename}.feat"
				execute "global_score=${score}"
                execute "query=/${name}"
				execute "file_result=$fich"
                execute "standar_out_file" #normaliza la salida en un json
                #                
                mv ${fich} ${folder_grid}
			done
			#rm ${out_molec}.sdf

		fi

		
	else
		echo "LS solo soporta VSB"
		exit
	fi
}

















#echo  el fichero existe  ${salida}.sdf
		#execute "nomp=`cat ${CWD}${target} | grep \"pharmacophore name=\" | awk -F"\"" '{print $2}'`"
#grep -Hrn '|.' *.sdf |awk -F: {'if (NF>=7)  print "cp "$1 " ..LSCHEMBL252557-4/"'}
#if [ error == 0 ];then #si la ejecucion del software es correcta se sigue
	#	l=`cat $salida.sdf |head -1` #nombre del query en la BBDD ZINC
	#	s=`cat $salida.sdf  | grep -A1 "<Score>" |grep -v ^" " |grep -v "Score" |grep -v ^--` #Score obteniado
	#	contador=0;
	#	for  i in $s; do
	#		debugC "ScriptLS: echo $i >>$salida""$contador"".txt"
	#		echo $i >>$salida""$contador"".txt
	#		contador=`expr $contador + 1`
	#	done
	#else
	#	cat ${salida}.ucam	|grep -i 'ERROR\|licence' > ${salida}.txt
	#fi
#/mnt/home/users/ac_001_um/helena/docking/externalSw/ligandscout4/iscreen -q /mnt/home/users/ac_001_um/helena/docking/pharmacophores/cox2_cel_str.pml -d /mnt/home/users/ac_001_um/helena/docking/queries/vp_2015_LS/vp_2015_LS_vpcode_library.ldb -a 2 -S relative -M 4 -C 4 -F 6501 -L 7000 -o /mnt/home/users/ac_001_um/helena/docking/queries/vp_2015_LS/screen_vp_2015_cox2STR_a2/screen_vp_2015_cox2STR_a2_6501_7000.sdf
#for i in `< $alida".txt"`; do echo -n ${i}" ";done; echo ""> $salida""Tmp.txt
#echo "split -l 1 $salida""Tmp.txt $salida_".txt""
#echo $l","$target","$s > $salida.sd
#echo $s $x $y $z $query 1 > $salida.txt

