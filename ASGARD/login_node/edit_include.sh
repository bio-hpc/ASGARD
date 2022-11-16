#!/bin/bash

DIR=$1 # Nombre del directorio donde se encuenta la carpeta de VS_GR
#INPUT=$2 # Nombre de la carpeta donde se encuentran las queries que la topologia redirige


actual_dir=$(pwd)/
topology_file=$(ls "$DIR"/*.top)
topology_include=$(cat "$topology_file" | grep include | head -n1)


pre='#include "'
post='shuttlemol/'

topology_include=${topology_include#${pre}*}
topology_include=${topology_include%%${post}*}"$post"
#
#echo $topology_include
#echo $actual_dir

for i in $(ls "$DIR")
do 
echo $i
sed -i "s|"$topology_include"|"$actual_dir"|g" "$DIR"/$i
done 
