#!/usr/bin/env bash

#____________________________________________________________________
#	Genera el grid para LF
#__________________________________________________________________

if [ ! -f ${out_grid}"_grid.bin" ];then 		# si no existe el fichero de grid se crea
if [ ${nom_serv} ==  "leftraru4" ];then
		export HOME=/mnt/flock/hperez
	fi

	echo "Generate Grid ($x $y $z) wait a few minuts "
	#datos que necesita el fichero.Par 	#grid-center=33,13,-6 	#grid-size=15,15,15
	debugC "generateGridLF: grid-center=${x},${y},${z} >${out_grid}.par"
	debugC "generateGridLF: grid-size=${gridSizeX},${gridSizeY},${gridSizeZ} >>${out_grid}.par"
	echo "grid-center=${x},${y},${z}" >${out_grid}.par
	echo "grid-size=${gridSizeX},${gridSizeY},${gridSizeZ}" >>${out_grid}.par

	debugC "generarteGridLF: ${path_external_sw}leadFinder/leadfinder \
	--grid-only \
	--protein=${CWD}${target} \
	--save-grid=${out_grid}grid.bin \
	-np ${cores} \
	--parameters=${out_grid}.par >> ${out_grid}_G.ucm"

	${path_external_sw}leadFinder/leadfinder\
	--grid-only \
	--protein=${CWD}${target} \
	--save-grid=${out_grid}_grid.bin \
	-np ${cores} \
	--parameters=${out_grid}.par >> ${out_grid}_G.ucm
	#source ${pathSD}/grids/spinner.sh
	rm ${out_grid}".par"
fi

