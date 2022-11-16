#!/bin/bash
#_________________________________________________________________________________________
#
#	Si el software necesita grid se genera
#_________________________________________________________________________________________
if [ ${grid} == "Y" ];then
	source  ${path_cluster_nodes}grids/generateGrid${software}.sh
fi 






















