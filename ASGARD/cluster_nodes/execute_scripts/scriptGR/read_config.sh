path_in=$(dirname "$query")"/" # queries/test_dm/
path_target=$(dirname "$target")"/" # queries/test_dm/
mode=""
name_query="${query%.*}"
name_query=`basename ${name_query}`

if [ -f  ${CWD}${path_in}${name_target}"_"${name_query}'_complex.conf' ]; then  # target query/ies
     input_topology_data=${CWD}${path_in}${name_target}"_"${name_query}'_complex'
elif [ -f ${CWD}${path_target}${name_target}'_target_complex.conf' ];then       # solo target
    patIn=$(dirname "$target")"/" # queries/test_dm/
    input_topology_data=${CWD}${patIn}${name_target}'_target_complex'

elif [ -f ${CWD}${patIn}$"biphsic_systems_"${name_query}'_complex.conf'  ];then  # biphsic_systems
    input_topology_data=${CWD}${patIn}$"biphsic_systems_"${name_query}'_complex'
else
    echo "ERROR: no se enecuentran los ficheros de topologia"
    exit
fi

mode_gr=`cat ${input_topology_data}.conf |grep Profile |awk -F\: '{print $2}'`

