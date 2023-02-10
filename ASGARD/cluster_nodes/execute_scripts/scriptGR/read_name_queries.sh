#!/usr/bin/env bash
#
#   Encuentra el nombre de los queries
#
count=0
for itp in "${itps[@]}"; do
    if [[ $itp != *_[a-z].itp ]];then
        name_queries[$count]=`cat $itp |grep  "\[ atoms \]" -A3  |tail -n 1 |awk '{print $4}'`
        count=`expr $count + 1`
    fi
done