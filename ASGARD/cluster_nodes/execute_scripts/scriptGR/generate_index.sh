#!/usr/bin/env bash
ejecutar()
{
	execute "echo \"Generate index\""
    generate_index
    execute "echo \"Fin Generate index\" "

}

generate_index()
{
    echo "${number_all_groups}" >  ${out_aux}_complex_index.ndx.tmp.tmp
    echo ${out_aux}_complex_index.ndx.tmp.tmp
    echo  ${out_aux}_complex_index.ndx.tmp.tmp
    echo "q" >>${out_aux}_complex_index.ndx.tmp.tmp
    execute "${prefix_gromacs} make_ndx${mpi} -f ${out_molec}_complex_solv_ions.gro  -o ${out_grid}_complex_index.ndx < ${out_aux}_complex_index.ndx.tmp.tmp "
    execute "rm  ${out_aux}_complex_index.ndx.tmp"
    execute "rm  ${out_aux}_complex_index.ndx.tmp.tmp"
}