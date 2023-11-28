#!/usr/bin/env bash
#
#
contador=0
while read line
do
    if [[  ${line:0:1} !=  ";" ]];then 

        if [[ $line == "#include"* ]] && [[ $line != *"force_field"* ]] && [[ $line != *"_porse"* ]];then 
            aux=`echo $line |awk '{print $2}'`
            aux=${aux%?}
            aux=${aux:1}
            itps[$contador]=$aux
            porse=`cat $aux |grep \#include |awk '{print $2}'`
            porse=${porse%?}
            porse=${porse:1}
            porse_itps[contador]=$porse
            contador=`expr $contador + 1`
        fi
        if [[ $line == *"_porse"* ]];then 
            echo $line
            aux=`echo $line |awk '{print $2}'`
            aux=${aux%?}
            aux=${aux:1}
            porse_itps[contador]=$aux
            contador=`expr $contador + 1`
        fi
     fi
done < ${out_molec}_complex.top
