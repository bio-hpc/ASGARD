#!/usr/bin/env bash
#
#	Pequeño resumen de los parametros de entrada de la pruebas
#
#echo "Output prefix: ${CWD}${patIn}${name_target}"_"${ligOut}"

printf '%50s\n' 'Resume input data'
echo ""
if [ -z ${gpu} ];then
    printf '%-20s %-50s\n' 'gpu:' No
else
    printf '%-20s %-50s\n' 'gpu:' ${gpu}
fi
printf '%-20s %-50s\n' 'ntomp:' ${cores}
printf '%-20s %-50s\n' 'cores:' ${cores}

echo ""
printf '%-20s %-50s\n' 'Gro:' ${out_molec}_complex.gro
printf '%-20s %-50s\n' 'Top:' ${out_molec}_complex.top

for itp in "${itps[@]}"; do
	printf '%-20s %-50s\n' 'itp:' ${itp}
done
for porse in "${porse_itps[@]}"; do
    printf '%-20s %-50s\n' 'porse:' ${porse}

done

for query in "${name_queries[@]}"; do
    printf '%-20s %-50s\n' 'Query:' ${query}
done
echo ""
printf '%-20s %-50s\n' 'file_conf_tpr:' ${file_conf_tpr}
printf '%-20s %-50s\n' 'file_conf_min:' ${file_conf_min}
printf '%-20s %-50s\n' 'file_conf_nvt:' ${file_conf_nvt}
ini=`expr $step_npt `
for ((j=1; j<5; j++));do
    printf '%-20s %-50s\n' "file_conf_npt_${ini}:" ${file_conf_npt}${ini}.mdp
	ini=`expr $ini + $step_npt`
done
printf '%-20s %-50s\n' 'file_conf_md:' ${file_conf_md}



echo ""
printf '%-20s %-50s\n' 'step_npt:' ${step_npt}
printf '%-20s %-50s\n' 'step_nvt:' ${step_nvt}
printf '%-20s %-50s\n' 'step_md:' ${step_md}
printf '%-20s %-50s\n' 'step_min:' ${step_min}
printf '%-20s %-50s\n' 'write_data:' ${write_data}
echo ""
printf '%-20s %-50s\n' 'Solvent:' ${solvent}
printf '%-20s %-50s\n' 'force_field:' ${force_field}
printf '%-20s %-50s\n' 'solvatation:' ${solvatation}
printf '%-20s %-50s\n' 'type_grid:' ${type_grid}
printf '%-20s %-50s\n' 'padding_grid:' ${padding_grid}
printf '%-20s %-50s\n' 'temp:' ${temp}
printf '%-20s %-50s\n' 'seedg:' ${seedg}
printf '%-20s %-50s\n' 'system_charge:' ${system_charge}

printf '%-20s %-50s\n' 'gromacs_version:' ${prefix_gromacs}
printf '%-20s %-50s\n' 'presure_npt:' ${pressure_npt}
echo ""


#file_conf_tpr
#printf '%-20s %-50s\n' 'query:' ${query}
#printf '%-20s %-50s\n' 'target:' ${target}
#printf '%-20s %-50s\n' 'mode:' ${mode}
#printf '%-20s %-50s\n' 'ph:' ${ph}







