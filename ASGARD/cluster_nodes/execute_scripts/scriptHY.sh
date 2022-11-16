#!/usr/bin/env bash
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
graph_global_field="CG3_Steric":"CG3_Clash":"CG3_ProDesolv":"CG3_LigDesolv":"CG3_LigDesolvHB":"CG4_HB":"Chemgauss4_score"
graph_global_color="b":"g":"r":"c":"y":"m":"k"

execute_script()
{
    execute "export OE_LICENSE=${path_external_sw}licenses/license${software}.txt"
    execute "${path_external_sw}oedocking/bin/hybrid \
        -receptor ${CWD}${target} \
        -dbase ${CWD}${query} $opt_aux \
        -prefix ${out_molec}_ &>${out_aux}.ucm"
    if [ $? == 0 ];then # si no falla la ejecucion
        #execute "tail -n 1 ${out_molec}_score.txt |awk '{print \$2\" '$x' '$y' '$z' '$query' \"1}' > ${outTxt}.txt"
        execute "${path_external_sw}tools/oe_convert ${out_molec}_docked.oeb.gz  ${out_molec}.mol2"
        execute "mv ${out_molec}_score.txt ${out_energies}.eng"
        execute "rm -f ${out_molec}_{settings.param,status.txt,report.txt}"

        execute "global_score=\"`cat ${out_energies}.eng |tail -1 |awk '{print $2}'`\""
        #"CG3_Steric":"CG3_Clash":"CG3_ProDesolv":"CG3_LigDesolv":"CG3_LigDesolvHB":"CG4_HB":"Chemgauss4_score"
        execute "graph_global_score=\"`cat ${out_energies}.eng | tail -1 |awk '{print $3\":\"$4\":\"$5\":\"$6\":\"$7\":\"$8\":\"$2  }'`\""
        coords=${x}":"${y}":"${z}
        execute "standar_out_file" #normaliza la salida en un json
    fi

}


