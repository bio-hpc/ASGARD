#!/usr/bin/env bash
OLDIFS=$IFS
IFS='-' read -r -a array <<< "$opt_aux"
for element in "${array[@]}"
do

	aux=`echo $element | cut -f1 -d' '`
	param=`echo $element | cut -f2 -d' '`
	case "$aux" in
		solvent)				solvent=$param;;
		force_field)			force_field=$param;;
		solvatation)			solvatation=$param;;
		step_npt)		        step_npt=$param;;
		step_nvt)		        step_nvt=$param;;
		step_md)		        step_md=$param;;
		step_min)		        step_min=$param;;
		write_data)				write_data=$param;;
		bt)						type_grid=$param;;
		padding_grid)			padding_grid=$param;;
		temp)				    temp=$param;;
		pressure_npt)           pressure_npt=$param;;
		seedg)					seedg=$param;;
		prefix_gromacs)		prefix_gromacs=$param
		                    groamcs=$prefix_gromacs;;


		#mode_gr)            mode_gr=${param};;

		#mode)					mode=$param;;
		#ph)						ph=$param;;
	esac
done
IFS=$OLDIFS;

