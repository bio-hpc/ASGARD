#!/bin/bash
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#   Description: Devuelve las nergias de la prueba
# ______________________________________________________________________________________________________________________



if [ "$#" -ne 1 ]; then
    echo "Error debe indicar la carpeta de laprueba"
    exit
fi

cat $1/energies/*.json |grep '\"global_score\"\|name' |sed 'N;s/\n/ /'  |sort -k 2
