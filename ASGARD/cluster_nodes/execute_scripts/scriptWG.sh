#!/usr/bin/env bash
execute_script()
{
  TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
  execute "${path_external_sw}wega/wega20160112 -outsd ${out_molec}.sd -qsd ${CWD}${target} -sd ${CWD}${query} $opt_aux  "    #ejecucion de wega, genera 3 ficheros sd, txt , bestalign
  if [ ${option} == "VS" ];then
      execute "score=`cat ${out_molec}.txt |grep -v \"target molecule\" |sort -t',' -n -k3|tail -1 | awk -F"," '{print $3}'`"     #fichero de score ${out_molec}.txt
      execute "namelig=`cat ${CWD}${query} | head -1 | awk -F"_" '{print $1}'`"
      execute "query=`cat ${CWD}${target} | head -1 | awk -F"_" '{print $1}'`"
      execute "rm ${out_molec}.txt"

      if  [[ $score  =~ $re ]] ; then
          execute "echo $score $x $y $z ${query} ${query}_${namelig} > ${out_energies}.en"
          round=`echo $score | python -c "print round(float(raw_input()),4)"` #este comando no deja debuygger
          if [ -f ${out_molec}_Aligned.sd ];then
              execute "${path_external_sw}tools/babel -m ${out_molec}_Aligned.sd ${out_molec}_Aligned_.sd"   #dividimos por si tiene conformaciones el best align ojo NO tiene hidrogenos
              execute "rm ${out_molec}_Aligned.sd"                               #eliminamo el fihceor
              execute "rm ${out_molec}.sd"                                       #eliminamos el otro fichero de salida de wg
              for i in ${out_molec}_Aligned_*.sd; do
              #echo "entro en el for para buscar la molecula alineada de menor score"
                  execute "tempscore=`grep \"<WEGAScore_query-\" -A1 -Hr $i | tail -1 | grep -Po 'sd-\K.*'`"
                  if (( $(echo $tempscore '==' $round | bc -l) )); then
                    execute  "mv $i ${out_molec}_BestAligned.sd"
                  else
                    execute "rm $i"
                  fi
                done
          fi
          create_out
      else
         execute "rm ${out_molec}*sd" #si ha fallado la prueba se eliminan sus ficheros
      fi
  elif [ ${option} == "VSB" ];then
    echo "VSB"
    #execute "${pathSW}wega/wega20160112 -outsd ${out_molec}.sd -qsd ${CWD}${target} -sd ${CWD}${query} $opt_aux  &> /dev/null "    #ejecucion de wega, genera 3 ficheros sd, txt , bestalign
    #mv ${out_molec}.txt ${outTxt}.txt

 fi

}
function create_out()
{

	file_result=${out_molec}_BestAligned.sd
	execute "global_score=\"` cat ${out_energies}.en| awk '{print $1}'`\""
	q=`basename $target`
	aux_target=${folder_molec}${q}
	if [ ! -f $aux_target ];then
	    cp ${target} ${folder_molec}
	fi
	target=$aux_target


	echo $target
    execute "standar_out_file" #normaliza la salida en un jsocn
}

#_______________________________ Antonio
#${CWD}scriptsDocking/externalSw/wega/wega2016 -outsd ${salida}.sd -qsd ${CWD}${proteina} -sd ${CWD}${ligando} -M 0 -ST 0 -CC true -SO2 true &> /dev/null
#cd ${directorio}
#python ${CWD}scriptsWeb/splitWG.py ${salida}_Aligned.sd
#cd ${CWD}


