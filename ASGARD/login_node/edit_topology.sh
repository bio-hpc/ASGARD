#!/bin/bash

DIR=$1 # Folder name where VS_GR result folder is found
QUERIES=$2 # Folder where the topology complex is

actual_dir=$(pwd)/ASGARD/external_sw/gromacs/force_field/amber99sb.ff/
topology_file=$(ls $DIR.top)

old_topology_include=$(cat $topology_file | grep '/queries' | head -n1)

topology_include=$(cat "$topology_file" | grep include | head -n1)
pre='#include "'
post='forcefield.itp'

topology_include=${topology_include#${pre}*}
topology_include=${topology_include%%${post}*}
  
sed -i "s|"$topology_include"|"$actual_dir"|g" $topology_file


# For posre (queries)

actual_dir=$(pwd)'/'

pre='#include "'
post='queries'

old_topology_include=${old_topology_include#${pre}*}
old_topology_include=${old_topology_include%%${post}*}

echo $old_topology_include
echo $actual_dir
sed -i "s|"$old_topology_include"|"$actual_dir"|g" $topology_file

# QUERIES INCLUDES

old_topology_include=$(cat $topology_file | grep '/queries' | head -n1)
old_topology_include=${old_topology_include#${pre}*}
post='"'
old_topology_include=${old_topology_include%%${post}*}
query_topology_include=$(echo $old_topology_include | rev | cut -d'/' -f2- | rev)
query_new_topology_include=$(echo $(pwd)/${QUERIES})

sed -i "s|"$query_topology_include"|"$query_new_topology_include"|g" $topology_file

echo $query_topology_include
echo $query_new_topology_include
echo $QUERIES