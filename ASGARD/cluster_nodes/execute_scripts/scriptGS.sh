#!/usr/bin/env bash
execute_script()
{

	TAG=`echo $(basename $BASH_SOURCE)` #e
    if [ "${ext_query}" != ".com" ];then # if ligand extension  is different to .com,  convert ligand
        query_out_ext="${query%.*}"
        query_conf=$(basename $query)
        query_conf="${query_conf%.*}"
        lines[0]="%chk=${query_conf}.chk"
        lines[1]="%mem=$mem"
        lines[2]="%nprocshared=${cores}"
        aux=`echo $opt_aux |sed 's/-header //' |sed 's/_/ /g'`
        lines[3]="#"${aux}
        #lines[4]=`echo $opt_aux |sed 's/-header //' |sed 's/_/ /g'`
        lines[4]='\\'                       # New line characeter
        lines[5]="${query_conf}"
        lines[6]='\\'
        #echo "babel $query ${query_out_ext}.com"
        babel $query ${query_out_ext}.com
        sed -i  1,4d ${query_out_ext}.com
        for ((i=${#lines[@]}-1; i>=0; i--)); do
            sed -i -e "1i${lines[$i]}" ${query_out_ext}.com
        done
        query=${query_out_ext}.com
        echo "" >> ${query_out_ext}.com
    fi
    g09 < ${query} > ${out_energies}.log

    if [ `cat ${out_energies}.log |wc -l` != "0" ];then
        babel ${out_energies}.log ${out_molec}.mol2
    fi



}


#smi="${query%.*}".smi
#	if [ -f $smi ];then
#
#		#
#		#	Convertimos el query original y lo guardamos en las parpetas de  la prueba
#		#
#		execute "dirLig=$(dirname "${out_molec}")"
#		execute "baseLig=`basename ${out_molec}`"
#		execute "babel $query ${out_molec}.com"
#		#
#		#	Comprobamos si existe el smi con la carga
#		#
#		execute "cargaTotal=\"`cat ${CWD}${smi} |awk '{print $1}'`\""
#  		execute "pos=\"`grep -o "+" <<< "$cargaTotal" | wc -l`\""
#  		execute "neg=\"`grep -o "-" <<< "$cargaTotal" | wc -l`\""
#  		execute "charge=\"`expr $pos - $neg`\""
#		#echo "babel ${out_molec}.com -omolreport "
#		sed -i '1,5d' ${out_molec}.com
#		sed -i "1s/^/${charge} 1\n/" ${out_molec}.com
#		sed -i '1s/^/\n/' ${out_molec}.com
#		sed -i "1s/^/$baseLig\n/" ${out_molec}.com
#		sed -i '1s/^/\n/' ${out_molec}.com
#		sed -i '1s/^/   POP=MK\n/' ${out_molec}.com
#		sed -i '1s/^/   B3LYP\/6-31G(d)\n/' ${out_molec}.com
#		sed -i '1s/^/\#p OPT FREQ\n/' ${out_molec}.com
#		sed -i "1s/^/\%nprocshared=${cores}\n/" ${out_molec}.com
#		sed -i "1s/^/\%mem=4000MB\n/" ${out_molec}.com
#		asx=`basename $out_molec`
#		sed -i "1s/^/\%chk=${asx}.chk\n/" ${out_molec}.com
#		echo "">>${out_molec}.com
#		cd $dirLig
#		g09 < ${baseLig}.com > ${out_energies}.en
#
#		check=`cat ${out_energies}.en |grep -n  "Fitting point charges to electrostatic potential" |wc -l`
#		if [ $check -gt 1 ];then
#			numIni=`cat ${out_energies}.en |grep -n  "Fitting point charges to electrostatic potential" |tail -1 |awk -F: '{print $1}'`
#			numIni=`expr $numIni + 4` #el fihcero siempre sale =
#			numFin=`cat ${out_energies}.en |grep -n  "Charges from ESP fit with hydrogens summed into heavy atoms:" |tail -1 |awk -F: '{print $1}'`
#			numFin=`expr $numFin - 1`
#
#			sed -n "${numIni},${numFin} p" ${out_energies}.en >${out_molec}Tmp.com
#			echo "python ${path_external_sw}gaussian/modCharges.py ${out_molec}Tmp.com ${query}"
#			execute "python ${path_external_sw}gaussian/modCharges.py ${out_molec}Tmp.com ${CWD}/${query}"
#		else
#			echo "ERROR La prueba no ha terminado bien"
#		fi
#
#
#
#	else
#		echo "Error debe existir un segundo fichero con tanos + y - como la carga de la target, o vacio en caso de 0"
#	fi

	#echo "cd $dirLig"
	#echo "g09 < ${baseLig}.com"
	#cd $dirLig
	#g09 < ${baseLig}.com
	#mv ${query}.com ${out_molec}.com
	#echo "cd $dirLig"
	#echo "g09 ${out_molec}.com"
