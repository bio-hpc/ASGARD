#!/bin/bash

DIR=$1 # Folder name where VS_GR results folder is found 

actual_dir=$(pwd)/
topology_file=$(ls "$DIR"/*.top)
topology_include=$(cat "$topology_file" | grep include | head -n1)


pre='#include "'
post='shuttlemol/'

topology_include=${topology_include#${pre}*}
topology_include=${topology_include%%${post}*}"$post"


for i in $(ls "$DIR")
do 
sed -i "s|"$topology_include"|"$actual_dir"|g" "$DIR"/$i
done 
