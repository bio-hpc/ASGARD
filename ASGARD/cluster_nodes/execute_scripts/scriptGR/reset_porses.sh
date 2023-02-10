#!/usr/bin/env bash
#
#	Antes de iniciar la simulacion nos aseguramos que estos ficheros esten correctos
#
change_porse_files()
{
    for porse in "${porse_itps[@]}"; do
		sed -i "s,[0-5]000  [0-5]000  [0-5]000,$1  $2  $3,g" $porse
    done
}