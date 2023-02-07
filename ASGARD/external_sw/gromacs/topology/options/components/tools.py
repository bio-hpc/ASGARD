#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import shutil
import subprocess
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def check_conformations_ligand_mol2(file_mol2):
    """
        devuleve las conformaciones de un liagndo mol2
    """
    tokens = [ "@<TRIPOS>ATOM","@<TRIPOS>BOND"]
    tokens_cnt = 0
    with open(file_mol2) as fp:
        for line in fp:
            if line.strip() in tokens :
                tokens_cnt += 1
    return tokens_cnt / 2


def write_file(lista, file):
    """
        Escribe una lista en un ficher
    """
    f = open(file, 'w')
    for i in lista:
        f.write(str(i)+"\n")
    f . close()


def rename_res_name_mol2(file_mol2, num_queries):
    """
        Renombea el nombre del residuo a L00 , L01 ...
        Asi eliminamos el error con el nombre

    """
    lst = []
    mode_file = False
    tokens = ["@<TRIPOS>ATOM", "@<TRIPOS>BOND"]
    tokens_cnt = 0
    with open(file_mol2) as fp:
        for line in fp:

            if line.strip() in tokens:
                tokens_cnt += 1
            if tokens_cnt == 1:
                if line.strip() != "":
                    aux = re.sub(' +', ' ', line).strip().split(" ")
                    if aux[0] not in tokens:

                        if len(aux) < 8:
                            line = line.strip() + "\tLIG"
                            mode_file = True
                        if len(aux) > 8:
                            if num_queries < 10:
                                line = line.replace(str(aux[7]), "L0" + str(num_queries))
                            else:
                                line = line.replace(str(aux[7]), "L" + str(num_queries))
                            mode_file = True
                    lst.append(line[:-1])
            else:
                lst.append(line[:-1])

        if mode_file:
            shutil.copyfile(file_mol2, file_mol2 + ".bk")
            write_file(lst, file_mol2)
        return num_queries + 1


def execute_cmd(cmd):
    try:
        out_cmd = subprocess.check_output(cmd, stderr=subprocess.STDOUT,shell=True).decode("utf8").strip()
    
        if "ERROR" in out_cmd:
            for line in out_cmd.split('\n'):
                if "ERROR" in line.upper() or "FALLO" in line.upper():
                    print(bcolors.FAIL + "ERROR: " + line.strip() + bcolors.ENDC)
                else:
                    print("ERROR: " + line.strip())
    
    except subprocess.CalledProcessError as exc:
        for line in exc.output.decode("utf8").split('\n'):
            if "ERROR" in line.upper() or "FALLO" in line.upper():
                print (bcolors.FAIL +"ERROR: "+line.strip() + bcolors.ENDC )
            else:
                print("ERROR: " + line.strip() )
        exit()
    return out_cmd.split('\n')


def extract_atoms_type(file_gro):
    """
        Devueleve los atoms_type que deben ponerse en el itp de la proteina, esto es encesario cuando se hacen
        DM de varios ligandos iguales con una proteinas

    """
    lst_atoms_tyoe = []
    lst_itp = []
    token_start = "[ atomtypes ]"
    token_end = "[ moleculetype ]"
    contador = 0
    f = open(file_gro)
    for i in f:
        if i.startswith(token_start) or i.startswith(token_end):
            contador += 1
        if contador == 1 and i.strip() != "":
            lst_atoms_tyoe.append(i.strip())
            lst_itp.append(i[:-1])
    f.close()
    comment_atoms_type_itp_ligand(lst_itp, file_gro)
    return lst_atoms_tyoe


def comment_atoms_type_itp_ligand(lst_itp, file_gro):
    """
        Añade caracter de comentario a los atoms type del itp del ligand asi no repetir con la prtoteina

    """
    for i in lst_itp:  # se comentasn las lineas de atom type
        i = i.replace("[ ", "\\[ ")  # para la etiquerta atomstype
        cmd = 'sed -i "s/' + i + '/;' + i + '/" ' + file_gro
        execute_cmd(cmd)


def print_error_parameters(txt):
    print("\n"+txt+" Este script necesita que se le pase una proteina en formato pdb y una carpeta de ligandos en formato mool2")
    print ("para proteina solo añadir opcion -p only")
    print ("para gromacs gmx -g gmx default(gmx_mpi)")
    print("Ojo gromacs debe estar instalado o cargado como modulo")
    print ("Profiles:")
    print("{:25}{:20}".format("\tTARGET_QUERY","Default"))
    print("{:25}".format("\tTARGET"))
    print("{:25}".format("\tTARGET_QUERIES"))
    print("{:25}".format("\tTARGET_ONE_QUERY"))
    print("{:25}".format("\tQUERIES"))
    print("{:25}".format("\tDNA_QUERY"))
    print("{:25}".format("\tBIPHSIC_SYSTEMS"))
    print("Ejemplo:")
    print("python ShuttleMol/external_sw/gromacs/topology/generate_topology.py -t targets/test/1le0.pdb -q queries/test/")
    print("python ShuttleMol/external_sw/gromacs/topology/generate_topology.py -t targets/test/1le0.pdb -q queries/test/ -p  TARGET_QUERIES")
    print("python ShuttleMol/external_sw/gromacs/topology/generate_topology.py -t targets/test/1le0.pdb -p TARGET_ONE_QUERY")
    print("python ShuttleMol/external_sw/gromacs/topology/generate_topology.py -t targets/test/1le0.pdb -p TARGET")
    print("python ShuttleMol/external_sw/gromacs/topology/generate_topology.py -q queries/test/ -p QUERIES")
    print("python ShuttleMol/external_sw/gromacs/topology/generate_topology.py -q queries/test/ -p BIPHSIC_SYSTEMS")

    print("\n")
    exit(2)