#!/bin/bash

DIR=$1 # Folder name where VS_GR folder is found  

actual_dir=$(pwd)/ASGARD/external_sw/gromacs/force_field/amber99sb.ff/
topology_file=$(ls $DIR.top)


topology_include=$(cat "$topology_file" | grep include | head -n1)
pre='#include "'
post='forcefield.itp'

topology_include=${topology_include#${pre}*}
topology_include=${topology_include%%${post}*}
  
sed -i "s|"$topology_include"|"$actual_dir"|g" $topology_file

