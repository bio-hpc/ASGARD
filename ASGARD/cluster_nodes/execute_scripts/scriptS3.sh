#!/usr/bin/env bash
TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
execute_script()
{


	execute "cd ${directorio}molecules/"
	# Preparación del query, si no lo tenemos ya: incluimos flexibilidad y generación de ambos métodos
	execute "cp $CWD$target ." #se copia el query a la carpeta de la prueaa
	execute "cp $CWD$query ." #Se copia el queries a la carpeta de la prueba
	# Screenings
	execute "${path_external_sw}ChemAxon/JChem/bin/screen3d s  -q ${name_target}.ser -t ${name_query}.ser  $opt_aux   &> /dev/null " # &> $salida.ucam
	# Archivo txt para generación de resultados
	execute "a=`tail -n1 ${name_target}_${name_query}_shape_*_screenOut.txt | awk '{ print \$3 }'`"
	#mv ${local_prot/.ser/}_${local_lig/.ser/}_shape_*_screenOut.txt ${local_prot/.ser/}_${local_lig/.ser/}_shape_screenOut.out
	execute "mv ${name_target}_${name_query}_shape_*aligned.mol2 ${salida}.mol2" 	#se mueve su salida al formato estandar de ShuttleMol
	execute "rm ${name_target}_${name_query}_shape_*_screenOut.txt"  #eliminamos basura
	execute "rm ${name_query}.ser"
	execute "echo ${a/,/.} $x $y $z $query 1 0 > ${out_energies}eng"
}

#_____________________________ Antonio
#traducir() {
#        echo "$(awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR-1]}1' ${directorio}lig_ser/${nomLigando}.mol2.txt *${nomLigando}*_screenOut.txt)" > *${nomLigando}*_screenOut.txt
#}
#    cd ${directorio}
#    # Preparación del query, si no lo tenemos ya: incluimos flexibilidad y generación de ambos métodos
#    cp ${CWD}${ligando} .                   # Se copia el ligando a la carpeta de la prueba
#    # Screenings
#    ${CWD}scriptsDocking/externalSw/ChemAxon/JChem/bin/screen3d s -q ${nomProteina}.ser -t ${nomLigando}.ser -oformat mol2 -writeQuery
#    traducir
#    input=${directorio}*${nomLigando}*_aligned.mol2
#    map=${directorio}lig_ser/${nomLigando}.mol2.txt
#    python ${CWD}scriptsWeb/splitS3.py ${input} ${map} ${programa}
#    ;;

