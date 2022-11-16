#!/usr/bin/env bash
SCRPIT_CONV=${path_extra_shuttlemol}used_by_shuttlemol/convert_to.py
execute_script()  #hace Docking en la posicion indicada con el query indicado
{
    TAG=`echo $(basename $BASH_SOURCE)` #En todos los comandos debe existir esta etiqueta
    file_name_query="${query%.*}"
    ext_out=${opt_aux//"-eout "/}
    ext_out=$(echo $ext_out| cut -d'-' -f 1)
    ext_out="$(echo -e "${ext_out}" | sed -e 's/^[[:space:]]*//')"
    ext_out=`trim $ext_out`
    out=${file_name_query}"."$ext_out

    if [ "$ext_out" == "ldb" ];then

        opt_aux=`echo $opt_aux | sed 's/-eout//g'`
        opt_aux=`echo $opt_aux | sed 's/ldb//g'`
        opt_aux=`echo $opt_aux | sed 's/-multiconformation//g'`
        execute "python2  $SCRPIT_CONV $query $out $opt_aux"
     else
        execute "python2  $SCRPIT_CONV $query $out"
      fi
    execute "mv $out ${folder_molec}"
}
trim()
{
    local trimmed="$1"
    # Strip leading space.
    trimmed="${trimmed## }"
    # Strip trailing space.
    trimmed="${trimmed%% }"
    echo "$trimmed"
}