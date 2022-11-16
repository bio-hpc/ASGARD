#!/usr/bin/env bash
execute_script()
{
   TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta

   execute "export OE_LICENSE=${path_external_sw}licenses/license${software}.txt"
   execute "${pathSW}posit/bin/posit \
        -receptor ${CWD}${target} \
        -in $CWD$query \
        -prefix ${out_molec}_ &> ${out_aux}.ucm"
   if [ $? == 0 ];then
      execute "tail -n 1 ${out_molec}_score.txt |awk '{print \$2\" '$x' '$y' '$z' '$query' \"1}' > ${outTxt}.txt"
      execute "${path_external_sw}tools/oe_convert ${out_molec}_docked.oeb.gz  ${out_molec}.mol2"
      execute "mv ${out_molec}_score.txt ${out_energies}.eng"
      execute "rm -f ${out_molec}_{settings.param,status.txt,report.txt}"
      ####execute "cat ${out_aux}.ucm  |grep \"Below minimum probability\" |awk '{print \$4}' > ${salida}.txt"
      ####execute "rm -f ${salida}_{settings.param,status.txt,report.txt}"
   fi
}
#   outTxt=${folderTxt}${salida}
#   out_grid=${folderGrid}${salida}
#   out_aux=${folderOutUcm}${salida}
#   out_molec=${folderMolec}${salida}
#   out_energies=${folderEnergy}${salida}


