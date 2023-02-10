#!/usr/bin/env bash
#
#Genera un cubo
#
ejecutar()
{	
echo $mode_gr 
    if [ $mode_gr !=  "BIPHSIC_SYSTEMS" ];then
        generate_cube
    else
        generate_biphsic_system
    fi
}
generate_cube()
{
    #
    #   The dodecahedron is re-wrapped within a triclinic representation. You can unwrap the PBC after you have created a .tpr file using trjconv -pbc mol -center but everything here looks fine. 
    #   How do I get the protein in the center of the box in gromacs? - ResearchGate. Available from: https://www.researchgate.net/post/How_do_I_get_the_protein_in_the_center_of_the_box_in_gromacs [accessed Mar 17, 2017].
    #
    execute "${prefix_gromacs} editconf${mpi} -f ${out_molec}_complex.gro -o ${out_molec}_complex_box.gro  -d ${padding_grid} -bt ${type_grid} -c "    # Generacion grid


	execute "${prefix_gromacs} solvate${mpi} -cp ${out_molec}_complex_box.gro -o ${out_molec}_complex_solv.gro -p ${out_molec}_complex.top "   				#rellena la grid  de disolvente
    #
    #   Movemos los ficheros
    #
    execute "mv ${out_molec}_complex.gro ${folder_out_ucm}"
    execute "mv ${out_molec}_complex_box.gro ${folder_out_ucm}"
}
generate_biphsic_system()
{

    boxsize=5               #esto debe podeser cambiar
    add_mols=1200
    file_gro_disolvent=${itps[0]}
    file_gro_disolvent="${file_gro_disolvent%.*}.gro"
    file_gro_mol=${itps[1]}
    file_gro_mol="${file_gro_mol%.*}.gro"
    topol=${out_molec}_complex.top

    file_radius="ShuttleMol/external_sw/gromacs/force_field/vdwradii.dat"

    solvent="spc216.gro"
    out_file_gro_mol_edit_conf=${file_gro_mol}_edit_conf.gro
    out_file_gro_mol_join=${file_gro_mol}_join.gro
    out_put_i_m=${file_gro_disolvent%.*}_insert_molecules.gro
    out_put_editconf_new=${file_gro_disolvent%.*}_editconf_new.gro
    out_solvent=${out_molec}_complex_solv.gro


    #
    #   copio fichero de radios para el carbono
    #
    #cp ${file_radius} .
    execute "${prefix_gromacs} insert-molecules -ci ${file_gro_disolvent} -nmol ${add_mols} -box ${boxsize} ${boxsize} ${boxsize} -o ${out_put_i_m} -seed ${seedg}"
    
    
    #execute "${prefix_gromacs} genconf -f ${file_gro_disolvent} -nbox ${boxsize} ${boxsize} ${boxsize} -o ${out_put_i_m} -seed ${seedg}"
    #gmx genconf -f chx.gro -nbox 8 8 8 -o chx_box.gro
    #
    #	Aqui a lo mejor hay que hacer minimizazion
    #
    execute "${prefix_gromacs}  editconf -f ${out_put_i_m} -o ${out_put_editconf_new} -box ${boxsize} ${boxsize} $(($boxsize * 2)) -center $(echo "$boxsize/2" | bc -l) $(echo "$boxsize/2" | bc -l) $(echo "$boxsize/2" | bc -l)"

    #
    #	Llenamos de solvent
    execute "${prefix_gromacs}  editconf -f ${file_gro_mol} -o ${out_file_gro_mol_edit_conf} -box ${boxsize} ${boxsize} $(($boxsize * 2)) -center $(echo "$boxsize/2" | bc -l)  $(echo "$boxsize/2" | bc -l) $(echo "$boxsize*2*0.75" | bc -l)"
    #
    #	Generamos un cubo e insertamos el sisolvente
    #
    ${prefix_gromacs}  solvate -cp ${out_file_gro_mol_edit_conf} -cs ${out_put_editconf_new} -o ${out_file_gro_mol_join}
    #
    #	Añadimos aguas
    #
    #rm "./vdwradii.dat"
    ${prefix_gromacs}  solvate -cp ${out_file_gro_mol_join} -cs ${solvent} -p ${topol} -o ${out_solvent}

    #num_disolvent=`cat ${out_solvent} |grep  L01 |grep CAA |wc -l`
    num_disolvent=`cat ${out_solvent}| grep L01 | awk '{print $(NF-4)}' |sort |uniq -c |head -1 |awk '{print $1}'`
    num_solvent=`cat ${out_solvent} |grep SOL | grep OW |wc -l`

    L01=`cat ${topol}|grep "; Compound        #mols" -A1 |tail -1`
    name_L01=`echo $L01 |awk '{print $1}'`

    L02=`cat ${topol}|grep "; Compound        #mols" -A2 |tail -1`
    echo "cat ${topol}|grep \"; Compound        #mols\" -A2 |tail -1  "
    name_L02=`echo $L02 |awk '{print $1}'`

    sed -i '$d'  ${topol}
    sed -i '$d'  ${topol}
    sed -i '$d'  ${topol}

    echo "$name_L02         1">> ${topol}
    echo "$name_L01       ${num_disolvent}"  >> ${topol}
    echo "SOL       ${num_solvent}"  >> ${topol}





}
