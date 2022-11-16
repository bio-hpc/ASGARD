#!/usr/bin/env bash

#____________________________________________________________________
#	Genera el grid para LF
#__________________________________________________________________
if [ ! -f ${out_grid}"_grid.bin" ];then 		# si no existe el fichero de grid se crea
    lead_finder_run=${path_external_sw}lead_finder_2018/leadfinder_x86_64
	echo "Generate Grid ($x $y $z) wait a few minuts "
	#datos que necesita el fichero.Par 	#grid-center=33,13,-6 	#grid-size=15,15,15
	debugC "generateGridLF: grid-center=${x},${y},${z} >${out_grid}.par"
	debugC "generateGridLF: grid-size=${gridSizeX},${gridSizeY},${gridSizeZ} >>${out_grid}.par"
	echo "grid-center=${x},${y},${z}" >${out_grid}.par
	echo "grid-size=${gridSizeX},${gridSizeY},${gridSizeZ}" >>${out_grid}.par

	debugC "generarteGridLF: ${lead_finder_run} \
	--grid-only \
	--protein=${CWD}${target} \
	--save-grid=${out_grid}grid.bin \
	-np ${cores} \
	--parameters=${out_grid}.par >> ${out_grid}_G.ucm"

	${lead_finder_run}\
	--grid-only \
	--protein=${CWD}${target} \
	--save-grid=${out_grid}_grid.bin \
	-np ${cores} \
	--parameters=${out_grid}.par >> ${out_grid}_G.ucm
	#source ${pathSD}/grids/spinner.sh
	rm ${out_grid}".par"
fi

