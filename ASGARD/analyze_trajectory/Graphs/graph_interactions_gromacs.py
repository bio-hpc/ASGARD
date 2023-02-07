#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    version: 1
#
#   Este script generara la graficas de todos los atomos que interaccionan proteina ligando,
#   1 diagrama cartesiano eje y todos los residuso eje x tiempo
#   2 histograma con los residuos eje x y la suma de energias eje y
#   3 histograma con los residuos eje x y las energias desglosadas
#   4 diagrama poseview
#   5 tabla latex
#
import subprocess
import sys
import os
import numpy as np
from operator import itemgetter # ordenar lista
import re #expresion regular

from GenerateGraph.GenerateGraph import GenerateGraph

generateGraph = GenerateGraph()

solvent = "SOL"
path_asgard= repr(sys.argv[0])												#busco la ruta de lanzador para aniadir los paquetes de python
path_asgard= path_asgard[1:path_asgard.rfind("ASGARD/")+len("/ASGARD")]
poseview_run = path_asgard + "external_sw/poseview/poseview"
#babel_run = path_asgard + "external_sw/babel/babel"
babel_run = "babel"


class Energy():
    def __init__(self, avarage_all, lst_all, last_all, add_coul, avarage_coul , last_coul, lst_coul, add_lj, avarage_lj, last_lj, lst_lj):
        self.avarage_all = avarage_all
        self.lst_all = lst_all
        self.last_all = last_all

        self.avarage_coul = avarage_coul
        self.last_coul = last_coul
        self.add_coul = add_coul
        self.lst_coul = lst_coul

        self.avarage_lj = avarage_lj
        self.last_lj = last_lj
        self.add_lj = add_lj
        self.lst_lj = lst_lj



def generate_poseview(pdb):
    prefix = '{}_{}_{}'.format(os.path.splitext(pdb)[0], mol_target_original_name, mol_query_original_name)
    cmd = 'cat {} | grep -v L[0-9][0-9] |grep -v Li[0-9] |grep -v NA |grep -v CL |grep -v {}  |grep -v TER |grep -v ENDMDL > {}'.format(pdb, solvent, prefix+"_target.pdb")
    subprocess.check_output(cmd, shell=True)
    cmd = 'cat {} | grep {} > {}'.format(pdb, mol_query_name, prefix+"_query.pdb" )
    subprocess.check_output(cmd, shell=True)
    cmd = '{} -ipdb {} -omol2 {} 2> /dev/null'.format(babel_run, prefix+"_target.pdb",  prefix+"_target.mol2")
    subprocess.check_output(cmd, shell=True)
    cmd = '{} -ipdb {} -omol2 {} 2> /dev/null'.format(babel_run,  prefix+"_query.pdb",prefix+"_query.mol2" )
    subprocess.check_output(cmd, shell=True)
    cmd = '{} -p {} -l {} -o {} -t  {}'.format(
        poseview_run,
        prefix+"_target.mol2",
        prefix+"_query.mol2",
        n_g_poseview,
        "\""+mol_target_original_name +" "+mol_query_original_name+  "\""
    )
    try:
      subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
      print("Poseview license failed")
      cmd = ' convert -size 560x560 xc:white -font Palatino-Bold -pointsize 50 -fill red -draw \"text 20,55 {}\" {}'.format('\'Poseview license failed\'',n_g_poseview)
      subprocess.check_output(cmd, shell=True)
    cmd = 'rm {}*'.format(prefix)
    #subprocess.check_output(cmd, shell=True)

def read_xvg_energy_file(file):
    """
    Retorna los valores del fichero xml en una matriz
    time, Coul-SR, LJ-SR, Coul-14, LJ-14 y se le añade un campo que es la suma de los 4 componenetes
    """
    with open(file) as f:
        return [re.sub(' +',' ',line).strip().split(" ") for line in f if not line.startswith("#") and not line.startswith("@")]
def get_fields_xvg(lst):
    """
    :return: media global, suma_coul, media_coul, last_coul, suma lj, media lj, last_lj
    """
    coul = 0
    lj = 0
    lst_coul = []
    lst_lj = []
    lst_all = []
    sum = 0
    cnt_elements = 0
    for i in lst:
        last_coul = float(i[1]) + float(i[3])
        last_lj = float(i[2]) + float(i[4])
        lst_coul.append(last_coul)
        lst_lj.append(last_lj)
        lst_all.append(last_coul + last_lj)
        sum += last_coul + last_lj
        coul += last_coul
        lj += last_lj
        cnt_elements += 1

    return Energy (sum/cnt_elements, lst_all, (last_coul + last_lj), coul, coul/cnt_elements, last_coul, lst_coul, lj, lj/cnt_elements, last_lj, lst_lj)

#
#   Devuelve el nomrbe del residuo
#
avail_resiues= []
def get_residue(file):
    res = os.path.splitext(file)[0]
    num = res[res.rfind('_')+1:]
    aux = res[:res.rfind('_')]
    res = aux[aux.rfind('_')+1:]
    if res + num not in avail_resiues and num +res not in avail_resiues:
        res = res + num
        avail_resiues.append(res)

    else:
        res = None
    return res

if len(sys.argv) != 9:
    print ("Parameters: ")
    print("\t 1º  Xvg Directory")
    print("\t 2º  Output Prefix")
    print("\t 3º  Name Target")
    print("\t 4º  Name Target Original")
    print("\t 5º  Name Query")
    print("\t 6º  Name Query Original")
    print("\t 7º  Pdb del sistema para posview")
    print("\t 8º  media para descartar las energias")
    exit()

folder_xvg      = sys.argv[1]							#directorio donde se encuentran los xvg
results_prefix  = sys.argv[2]									#salida de la grafica
mol_target_name = sys.argv[3]								#fichero de proteina para hacer el poseview
mol_target_original_name = sys.argv[4]
mol_query_name  = sys.argv[5]								#nom del ligando L1, L2...
mol_query_original_name = sys.argv[6]
system_pdb = sys.argv[7]
energies_discard = float(sys.argv[8])
lst_energies = []
aux_target = []



#
#   Titulo 3 primeras grafcicas
#
title = "Gromacs Energies "+mol_target_original_name+" vs "+mol_query_original_name
x_title = "Time (ps)"
y_title = "(kJ/mol)"
#
#   Salida images
#
n_g_global_line_res = results_prefix + "_line_global_energy_res.png"
n_g_global_hist_res = results_prefix + "_hist_global_energy_res.png"
n_g_split_hist_res = results_prefix + "_hist_split_energy_res.png"
n_g_poseview = results_prefix + "_poseview.png"
n_t_latex = os.path.dirname(os.path.dirname(results_prefix))+"/"
n_t_latex += os.path.basename(results_prefix) + "_table_interations.tex"
n_g_join_hist = results_prefix + "_hist_global_split_energy_res.png"
n_g_join_line_pose = results_prefix + "_line_poseview.png"


#
#	Recogemos los datos de la simulacion, no group
#
for file in sorted(os.listdir(folder_xvg)): #devuelve los ficheros ordenados
    if not file.startswith("#") and "no_group" in file:
        name_residue = get_residue(file)
        if name_residue != None:
            interactions_xvg = read_xvg_energy_file(folder_xvg+file)
            aux =  get_fields_xvg(interactions_xvg)

            if abs(aux.avarage_all) > energies_discard:
                if 'Protein' in name_residue:
                    aux_target = (['Protein', aux])
                else:

                    lst_energies.append([name_residue, aux])

lst_step_md = [float(i[0]) for i in interactions_xvg]
lst_energies = sorted(lst_energies, key=itemgetter(0))
#if aux_target:
#    lst_energies.insert(0, aux_target)      # Aqu� a�ade Protein a la gr�fica (Quitar y sustituir unicamente por el valor)
    #lst_energies.append(aux_target)
#
#   Grafica 1 energia por atomos en el tiempo
#
name_residues = [row[0] for row in lst_energies]
y = [row[1].lst_all for row in lst_energies]
generateGraph.line_graph(name_residues, lst_step_md, y, n_g_global_line_res, x_title, y_title, title, "")
#
#   Grafica 2 histograma ccon energias global de los residuos
#
datos = []
datos.append([i[1].avarage_all for i in lst_energies])
datos.append([i[1].last_all for i in lst_energies])
legend = ["Mean","Last Step" ]
generateGraph.generate_multiple_bar(legend, datos,name_residues, n_g_global_hist_res, y_title ,title)
#
#   Grafica 3histograma ccon energias split de los residuos
#
datos = []
datos.append([i[1].avarage_coul for i in lst_energies])
datos.append( [i[1].last_coul for i in lst_energies ])
datos.append([i[1].avarage_lj for i in lst_energies])
datos.append([i[1].last_lj for i in lst_energies])
legend = ["Coul-SR+Coul-14","Last Steep Coul-SR+Coul-14","LJ-SR+LJ-14","Last Steep LJ-SR+LJ-14"]
generateGraph.generate_multiple_bar(legend, datos,name_residues ,n_g_split_hist_res, y_title, title)

#
#	4 diagrama poseview
#
generate_poseview(system_pdb)
#
# Tabla latex con energias
#
f = open(n_t_latex, "w")
f.write('{}'.format('\\begin{tabular}{ l r l r l r l r } \n') )
f.write('\t\multicolumn{6}{c} {Energies '+mol_target_original_name +' '+mol_query_original_name+'} \\\\'+"\n")
f.write('\t\hline \n' )
f.write('\t\tResidue & Energy & LJ &  Energy & Coul & Total  & Energy  & last step\\\\ \n')
f.write('\t\hline \n' )
round_decimals = 1
sum_fileds = [0, 0, 0, 0, 0, 0, 0]
for i in lst_energies:
    f.write( '\t\t{:<10} & {:<8} & $\pm {:<8}$ & {:<8} & $\pm {:<8}$ & {:<8} & $\pm {:<8}$ & {:<8} \\\\ \n'.format(
        i[0],
        round(i[1].avarage_lj, round_decimals),
        round(np.std(i[1].lst_lj), round_decimals),

        round(i[1].avarage_coul, round_decimals),
        round(np.std(i[1].lst_coul), round_decimals),

        round(i[1].avarage_all, round_decimals),
        round(np.std(i[1].lst_all), round_decimals),

        round(i[1].last_all, round_decimals)
        ))
    sum_fileds[0] += i[1].avarage_lj
    sum_fileds[1] += np.std(i[1].lst_lj)
    sum_fileds[2] += i[1].avarage_coul
    sum_fileds[3] += np.std(i[1].lst_coul)
    sum_fileds[4] += i[1].avarage_all
    sum_fileds[5] += np.std(i[1].lst_all)
    sum_fileds[6] += i[1].last_all

f.write( '\t\t{:<10} & {:<8} & $\pm {:<8}$ & {:<8} & $\pm {:<8}$ & {:<8} & $\pm {:<8}$ & {:<8} \\\\ \n'.format(
    "",
    round(sum_fileds[0], round_decimals), round(sum_fileds[1], round_decimals), round(sum_fileds[2], round_decimals), round(sum_fileds[3], round_decimals),round(sum_fileds[4], round_decimals),
    round(sum_fileds[5] , round_decimals),round(sum_fileds[6], round_decimals)
))
f.write('\\end{tabular}'+"\n")

f.close()
