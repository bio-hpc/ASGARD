#!/bin/bash

DIR=$1 # Nombre del directorio donde se encuenta la carpeta de VS_GR (mas tarde quiza unicamente el nombre de la topologia)

actual_dir=$(pwd)/ASGARD/external_sw/gromacs/force_field/amber99sb.ff/
topology_file=$(ls $DIR.top)


topology_include=$(cat "$topology_file" | grep include | head -n1)
pre='#include "'
post='forcefield.itp'

topology_include=${topology_include#${pre}*}
topology_include=${topology_include%%${post}*}
  
sed -i "s|"$topology_include"|"$actual_dir"|g" $topology_file


# ver que hacer con query y posre query.itp porse.itp. y cambiar de aqui tambien la ruta de query