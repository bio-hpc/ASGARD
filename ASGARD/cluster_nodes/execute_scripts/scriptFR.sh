#!/usr/bin/env bash
graph_global_field="CG3_Steric":"CG3_Clash":"CG3_ProDesolv":"CG3_LigDesolv":"CG3_LigDesolvHB":"CG4_HB":"Chemgauss4_score"
graph_global_color="b":"g":"r":"c":"y":"m":"k"
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
execute "export OE_LICENSE=${path_external_sw}licenses/license${software}.txt"
execute_script()
{
	#execute "coords=\"`${path_extra_shuttlemol}used_by_shuttlemol/get_center_ligand.py ${out_molec}.mol2`\""
    #create_out
    #exit
    #   Si la option es BD tienes que genrar la target indicando el numero de residuo tipo y cadena
    if [ $option == "BD" ];then # option BD
        convert_receptor
        execute "${path_external_sw}oedocking/bin/fred \
           -receptor ${CWD}${target} \
            -dbase ${CWD}${query} $opt_aux \
            -prefix ${out_molec}_ &> ${out_aux}.ucm"
        execute "${path_external_sw}tools/oe_convert ${out_molec}_docked.oeb.gz  ${out_molec}.mol2"
        execute "coords=\"`${path_extra_shuttlemol}used_by_shuttlemol/get_center_ligand.py ${out_molec}.mol2`\""

    else                    # option VS
        execute "${path_external_sw}oedocking/bin/fred \
           -receptor ${CWD}${target} \
            -dbase ${CWD}${query} $opt_aux \
            -prefix ${out_molec}_ &> ${out_aux}.ucm"
        execute "${path_external_sw}tools/oe_convert ${out_molec}_docked.oeb.gz  ${out_molec}.mol2"
        coords=$x":"$y":"$z
    fi   
    create_out

}
function convert_receptor()
{
    aux=${out_aux}_receptor.oeb  #   out del nuevo receptor
    if [ ${chain} == 'N/A' ];then
        chain=''
    fi
    execute "${path_external_sw}oedocking/bin/apopdb2receptor  -pdb ${CWD}${target} -site_residue:${resName}${num_amino_acid}${chain} -receptor ${aux}"
    target=${aux#${CWD}}    #   Se cambia el target al nuevo receptor
}

function create_out()
{
    if [ $? == 0 ];then # si no falla la ejecucion
        execute "mv  ${out_molec}_docked.oeb.gz ${out_aux}_docked.oeb.gz"                           #  movemos la molecula para que solo se queden los mol2 en el directorio molecula
        execute "mv ${out_molec}_score.txt ${out_energies}.eng"
        execute "rm -f ${out_molec}_{settings.param,status.txt,report.txt}"
        file_result=${out_molec}.mol2
        columns=`cat ${out_energies}.eng |tail -1  |wc -w`
        #
        #   El fichero de salida si no tiene nombre el ligando solo contiene enrgias y si tiene continene nombre mas energias
        #
        if [ $columns  -eq 8 ];then
            execute "global_score=\"`cat ${out_energies}.eng |tail -1 |awk '{print $2}'`\""
        else
            execute "global_score=\"`cat ${out_energies}.eng |tail -1 |awk '{print $2}'`\""
        fi
        graph_global_score
        execute "standar_out_file" #normaliza la salida en un json
    fi
}
function graph_global_score()
{
     #"CG3_Steric":"CG3_Clash":"CG3_ProDesolv":"CG3_LigDesolv":"CG3_LigDesolvHB":"CG4_HB":"Chemgauss4_score"
     columns=`cat ${out_energies}.eng |tail -1  |wc -w`
     if [ $columns  -eq 8 ];then
        execute "graph_global_score=\"`cat ${out_energies}.eng | tail -1 |awk '{print $3\":\"$4\":\"$5\":\"$6\":\"$7\":\"$8\":\"$2  }'`\""
     else
        execute "graph_global_score=\"`cat ${out_energies}.eng | tail -1 |awk '{print $2\":\"$3\":\"$4\":\"$5\":\"$6\":\"$7\":\"$1  }'`\""
     fi

}

#
#	Add para obtener las moleculas en sdf
#
#execute "${path_external_sw}tools/oe_convert ${out_aux}_docked.oeb.gz  ${out_molec}.sdf"
#sed  -i '$ i\> <SCORE>' ${out_molec}.sdf
#sed  -i "$ i\\${score}\n" ${out_molec}.sdf
#fi
#if [ $? == 0 ];then # si no falla la ejecucion
# execute "tail -n 1 ${out_molec}_score.txt |awk '{print \$2\" '$x' '$y' '$z' '$query' '${numAuxEjecucion}' \"}' > ${outTxt}.txt"
# execute "${path_external_sw}tools/oe_convert ${out_molec}_docked.oeb.gz  ${out_molec}.mol2"
# execute "mv ${out_molec}_score.txt ${out_energies}.en"
# execute "rm -f ${out_molec}_{settings.param,status.txt,report.txt}"
#
#   Add para guardar el score dentro del sdf
#
#execute "${path_external_sw}tools/oe_convert ${out_aux}_docked.oeb.gz  ${out_molec}.sdf"
#execute "score=\"`head -n 2  ${out_energies}.en | tail -n 1|awk '{print $2}'`\" "
#sed  -i '$ i\> <SCORE>' ${out_molec}.sdf
#sed  -i "$ i\\${score}\n" ${out_molec}.sdf
#fi
