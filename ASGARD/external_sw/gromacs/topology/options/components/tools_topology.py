#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from os.path import join

from ..ClassQuery import ClassQuery
from .tools import *


def create_target_queries_gro(lst_class_queries,  prefix_out, cfg):
    """
        Crea un gro con la proteina y el ligando llamado target-query_complex.gro

    """
    if len(lst_class_queries) == 0:
        print (cfg.target_file_gro)
        shutil.move(cfg.target_file_gro, prefix_out+cfg.prefix_sumulation+".gro")
        return
    else:
        sum_atoms_queries = 0
        for class_query in lst_class_queries:
            sum_atoms_queries = sum_atoms_queries + len(class_query.lst_query_gro)
        
        if cfg.target_file_gro != None:
            last_line_gro, lst_target_gro = get_lst_gro(cfg.target_file_gro)
           
            last_rest = int(lst_target_gro[len(lst_target_gro) - 1][0]) + 1  # ulytimo de la lista +1
            next_number_index = int(lst_target_gro[len(lst_target_gro) - 1][3]) + 1  # ulytimo de la lista +1
            
        else:
            last_line_gro = "   20.57000    24.39800    24.52000"
            lst_target_gro = []
            last_rest = 1
            next_number_index = 1

        number_atoms = sum_atoms_queries + len(lst_target_gro)   # suma de los atomos prot lig
        complex_gro = []  # Fichero gro en conjunto target y query
        complex_gro.append("Grunge ROck MAChoS")  # linea de placebo
        complex_gro.append(number_atoms)

        lst_aux = lst_target_gro
        for class_query in lst_class_queries:
            for i in range(len(class_query.lst_query_gro)):
                class_query.lst_query_gro[i][0] = last_rest
                class_query.lst_query_gro[i][3] = next_number_index
                next_number_index += 1
            last_rest = int(class_query.lst_query_gro[len(class_query.lst_query_gro) - 1][0]) + 1  # ulytimo de la lista +1
            next_number_index = int(class_query.lst_query_gro[len(class_query.lst_query_gro) - 1][3]) + 1  # ulytimo de la lista +1
            lst_aux = lst_aux + class_query.lst_query_gro
        for i in lst_aux:
            complex_gro.append((cfg.format_gro.format(i[0], i[1], i[2], i[3], float(i[4]), float(i[5]), float(i[6]))))
        complex_gro.append(last_line_gro)
        write_file(complex_gro, prefix_out + cfg.prefix_sumulation + ".gro")


def include_itps_atoms_query(lst_class_queries):
    lst = []
    lst_incluyed_atoms = []
    if len(lst_class_queries)>0:
        lst.append("; Atoms type queries")
        for class_query in lst_class_queries:
            for i in class_query.lst_atoms_type:
                if i not in lst_incluyed_atoms:
                    lst_incluyed_atoms.append(i)
                    lst.append(i)
        lst.append('')
        lst.append("; Include Query topology")
        for class_query in lst_class_queries:
            lst.append("#include \"" + class_query.prefix_out_q + ".itp\"")
    return lst

def new_topology_queries( prefix_out, cfg, lst_class_queries):
    lst_top = []
    lst_top.append(";This file was generated by generate_topogy.py")
    lst_top.append('')
    lst_top.append('; Include forcefield parameters')
    lst_top.append('#include "{}"'.format(os.path.join(cfg.path_forece_files_choice+'forcefield.itp' )))

    lst_top.append('')
    lst_top = lst_top + include_itps_atoms_query(lst_class_queries)
    lst_top.append('')
    lst_top.append('; Include water topology')
    lst_top.append('#include "{}"'.format(os.path.join(cfg.path_forece_files_choice+cfg.solvent+'.itp') ))
    lst_top.append('')
    lst_top.append('# ifdef POSRES_WATER')
    lst_top.append('; Position restraint for each water oxygen')
    lst_top.append('[position_restraints]')
    lst_top.append(';  i   funct   fcx   fcy   fcz')
    lst_top.append('1    1       1000       1000       1000')
    lst_top.append('#endif')
    lst_top.append('')
    lst_top.append('; Include topology for ions')
    lst_top.append('#include "{}"'.format(os.path.join(cfg.path_forece_files_choice+'ions.itp' )))
    lst_top.append('')
    lst_top.append('[ system ]')
    #lst_top.append('; Name')
    #lst_top.append('No Name')
    lst_top.append('')
    lst_top.append('[ molecules ]')
    lst_top.append('; Compound        #mols')
    for class_query in lst_class_queries:
        lst_top.append('{0:19} {1:20}'.format(class_query.name_query, "1").strip() )

    write_file(lst_top, prefix_out + cfg.prefix_sumulation + ".top")
    print(prefix_out + cfg.prefix_sumulation + ".top")

    # include "/home/horacio/jorge/lanzador/lanzador/externalSw/gromacs/campoFuerza/amber99sb.ff/forcefield.itp"


#def mod_target_top(target_file_top, prefix_out, name_query, lst_atoms_type, prefix_out_q):
def mod_target_top_queries(prefix_out, cfg, lst_class_queries):
    """
        Modifica el fichero tpo de la proteina a??adiendo el ligando los atomos y
        modificacnndo las rutas de los ficheros a absolutas

    """
    lst_top = []
    first = True
    contador = 0
    with open(cfg.target_file_top, "r") as myfile:
        for line in myfile:
            line = line.strip()
            #   Cambiamos las ruts del cmapo de fuerza a las nuestras propias
            if '#include' in  line:
                if 'forcefield.itp' in line or cfg.solvent+'.itp' in line or 'ions.itp' in line:
                    line = line.strip().split(" ")
                    aux = line[1][1:-1] # se eiminan las comillas
                    aux = os.path.basename(aux)
                    line = line[0]+" \""+cfg.path_forece_files_choice + aux + "\""
                    lst_top.append(line)
                    if first: # a??adimos el itp del ligando despues del primer include
                        lst_top = lst_top + include_itps_atoms_query(lst_class_queries)
                        first = False
                elif 'porse.itp' in line:
                    line = line.strip().split(" ")
                    aux = line[1][1:-1] #se eiminan las moillas
                    shutil.copyfile(aux,prefix_out+cfg.prefix_sumulation+"_target_porse.itp" )
                    lst_top.append("#include \""+cfg.path+prefix_out+cfg.prefix_sumulation+"_target_porse.itp\"")
                elif 'Protein_chain' in line:     #se supone que solo puede existir estas lineas en las topologias compuestas de varias cadenas
                    line = "#include \"" + cfg.itp_target_chain[contador]+"\""
                    lst_top.append(line)
                    contador += 1
                elif 'DNA_chain' in line:     #se supone que solo puede existir estas lineas en las topologias compuestas de varias cadenas
                    line = "#include \"" + cfg.itp_target_chain[contador]+"\""
                    lst_top.append(line)
                    contador += 1
            else:
                lst_top.append(line)
    write_file(lst_top, prefix_out+cfg.prefix_sumulation+".top")
    with open(prefix_out+cfg.prefix_sumulation+".top", "a") as myfile:
        for class_query in lst_class_queries:
            myfile.write('{0:19} {1:20}'.format(class_query.name_query, "1").strip()+"\n")



def get_query_charge(itp_query):
    f = open(itp_query)
    contador = 0
    sumaChargeLig = 0
    for i in f:
        if i.strip() == "[ bonds ]":
            contador += 1
        if contador == 1 and not i.startswith(";") and re.sub(' +', ' ', i).strip() != "":
            aux=re.sub(' +',' ',i).strip().split(" ")
            sumaChargeLig += float(aux[6])
        if i.strip() == "[ atoms ]":
            contador += 1
    return round(float(sumaChargeLig), 3)

def create_charge_file(charge_target, lst_class_queries, prefix_out, cfg):
    lst = []
    total_charge =float(charge_target)
    lst.append("Charge Target: "+str(round(float(charge_target), 3)))

    for class_query in lst_class_queries:
        lst.append('Charge '+class_query.name_query+': '+str(class_query.charge) )
        total_charge += class_query.charge
    lst.append("Suma: "+str( total_charge ))
    write_file(lst,prefix_out+cfg.prefix_sumulation+".eng")


def create_conf_file( lst, prefix_out, cfg):
    lst.append("Profile:"+cfg.profile)
    write_file(lst, prefix_out+cfg.prefix_sumulation+".conf")


def get_lst_gro(file_gro):
    """
        return:
          lst atoms gro
          last Line file "last_line_gro"
    """
    lst_gro = []
    first_line = True
    with open(file_gro) as myfile:
        for line in myfile:
            if not first_line:
                aux = re.sub(' +', ' ', line).strip().split(" ")
                if len(aux) == 7 or len(aux) == 6: #7 tienen los liagndos 6 las protienas
                    lst_gro.append([line[0:5].strip(), line[5:10].strip(), line[10:15].strip(), line[15:20].strip(), line[20:28].strip(), line[28:36].strip(), line[36:44].strip()])
                if len(aux) == 5:
                	lst_gro.append([line[0:5].strip(), line[5:10].strip(), line[10:15].strip(), line[15:20].strip(), line[20:28].strip(), line[28:36].strip(), line[36:44].strip()])
            else:
                first_line=False
        last_line_gro = line.strip()

    return last_line_gro, lst_gro


def generate_porse(prefix_out, lst_gro, cfg):
    """
        Genera fichero con restricciones de movimmmiento porse.itp
    """
    lst_porse = []
    lst_porse.append("; In this topology include file, you will find position restraint")
    lst_porse.append("; entries for all the heavy atoms in your original pdb file.")
    lst_porse.append("; This means that all the protons which were added by pdb2gmx are")
    lst_porse.append("; not restrained.")
    lst_porse.append("[ position_restraints ]")
    lst_porse.append("; atom  type      fx      fy      fz")
    for i in range(len(lst_gro)):
        lst_porse.append(cfg.format_itp.format(i+1, "1", "1000", "1000", "1000"))
    write_file(lst_porse, prefix_out+"_porse.itp")
    return prefix_out+"_porse.itp"


def create_ligand_topology(cfg, file_query):
    """
        Convierte el ligando a formato gromacs con ambertools y acepype
        Ojo a acepype no se le puede indicar la ruta de salida, para solucionarlo creamos una carpeta y nos movemos a ella
    """

    cfg.num_queries = rename_res_name_mol2( file_query, cfg.num_queries)
    name_query, ext_query = os.path.splitext(os.path.basename(file_query))

    folder_aux_query = name_query + "/"
    if not os.path.exists(folder_aux_query):
        os.mkdir(folder_aux_query)
    else:
        shutil.rmtree(folder_aux_query)
        os.mkdir(folder_aux_query)
    os.chdir(folder_aux_query)
    shutil.copyfile(join(cfg.path,cfg.dir_queries,name_query+ext_query), name_query+ext_query)
    cmd = 'ln -s {} acepype.py'.format(cfg.acpype)
    execute_cmd(cmd)
    cmd = '{} acepype.py  -i {} -o gmx -c user'.format(cfg.python_run, name_query+ext_query)
    execute_cmd(cmd)


    prefix_in = join(name_query+'.acpype',name_query)
    prefix_out_q = join(cfg.path, cfg.dir_queries,  cfg.name_target + '_' + name_query + cfg.prefix_sumulation + '_query')  # prefix salida ficheros query
    shutil.copyfile(prefix_in+'_GMX.gro', prefix_out_q + '.gro')
    shutil.copyfile(prefix_in+'_GMX.itp', prefix_out_q + '.itp')

    cmd = '{}  {} {}'.format(cfg.python_run, cfg.mod_diedrals_phosphate, prefix_out_q + '.itp')
    execute_cmd(cmd)
    #
    # A??adimos en el itp el ficehro porse
    #
    with open(prefix_out_q+'.itp', "a") as myfile:
        myfile.write("; Include Position restraint file\n")
        myfile.write("#ifdef POSRES\n")
        myfile.write('#include "' + prefix_out_q + '_porse.itp" \n')
        myfile.write("#endif\n")

    os.chdir(cfg.path)
    shutil.rmtree(folder_aux_query)
    class_query = ClassQuery()
    class_query.lst_atoms_type = extract_atoms_type(prefix_out_q + ".itp")
    class_query.last_line_gro, class_query.lst_query_gro = get_lst_gro(prefix_out_q + '.gro')
    class_query.prefix_out_q = prefix_out_q
    class_query.porse = generate_porse(class_query.prefix_out_q, class_query.lst_query_gro , cfg)
    class_query.itp = class_query.prefix_out_q + ".itp"
    class_query.charge = get_query_charge(class_query.itp)
    class_query.name_query = name_query

    #os.remove(join(prefix_out_q + ".gro"))
    return class_query



def create_topology_target(cfg):
    """
        Crea la topologia de la proteina mediante gromcas
    """

    os.chdir(cfg.path_force_field)
    cmd = "echo "+str(cfg.option_force_field)+" | {} pdb2gmx -v -f {} -o {} -i {} -p {} -water {} {}".format(cfg.gromacs_run, cfg.path + cfg.target,
                              cfg.target_file_gro, cfg.target_file_itp, cfg.target_file_top, cfg.solvent, cfg.ignh)

    print(cfg.path_force_field)
    print(cmd)
    for line in execute_cmd(cmd):
        if line.startswith("Total charge"):  # OJO NO PROBADO CON MONOMEROS
            charge = line.split(" ")[2]
        if line.startswith("Total charge in system"):
            charge = line.split(" ")[4]
    os.chdir(cfg.path)

    return charge
