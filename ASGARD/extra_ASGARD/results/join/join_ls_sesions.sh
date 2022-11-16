#!/bin/bash
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#   Description: Une los diferentes sdf cuando lanzas LS de 0 a 5
# ______________________________________________________________________________________________________________________


if [ "$#" -lt 1 ]; then
	echo "Error debe indicar las carpetas de la prueba"
	echo "sh join_ls_sesions.sh  LS_LB_D-LA_wexcl_v1_DBALLv503_a_"
        echo ""
        exit

	exit
fi
echo > ${1}all.sdf
for i in $(ls -d ${1}* |grep -v .sdf);do
	echo "cat ${i}/molecules/*.sdf >> ${1}all.sdf"
	cat ${i}/molecules/*.sdf >> ${1}all.sdf
done


