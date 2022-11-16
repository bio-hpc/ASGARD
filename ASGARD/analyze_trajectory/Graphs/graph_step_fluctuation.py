#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Genera grafica de fluctuacion cada x steps (Heat Map)
#
import subprocess
import sys
import numpy
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()

GRAPIHC_CA = 1 #protein


def isfloat(value):
    return str(value).replace('.','',1).isdigit()


if len(sys.argv) != 8:
    print("El script necesita:")
    print("1º una carpeta donde almacenar los tmp")
    print("2º fichero minimizacion.tpr")
    print("3º fichero simulacion.xtc")
    print("4º comando gormacs (dependiendo de la  version cambia la forma de llamarlo, ejemplo g_rmsF_mpi)")
    print("5º paso para generar las graficas")
    print("6º maximo para generar las graficas")
    print("7º Output grafica png")
    exit()
folder_fluctuation = sys.argv[1]
min_tpr = sys.argv[2]
md_xtc = sys.argv[3]
cmd_gromacs = sys.argv[4]
step = int(sys.argv[5])
max_steps = int(sys.argv[6])
out_png = sys.argv[7]
x_ticks = []
y_ticks = []
lst_files = []
lst_steps = []

for i in range(0, max_steps):

    ini = i*step
    fin = (i * step )+step
    lst_steps.append( fin)
    out_file = folder_fluctuation + "/" + str(ini) + "_" + str(fin) + ".xvg"
    cmd = 'echo {} | {} -s {} -f {} -b {} -e {} -o {} -res'.format(GRAPIHC_CA, cmd_gromacs, min_tpr, md_xtc, str(ini), str(fin), out_file)
    subprocess.check_output(cmd, shell=True)
    lst_files.append(out_file)
    x_ticks.append(fin)

#
#  matriz de 0s con el tamaño de filas y columan
#
n_res_aux, rms_f, _, _, _, _ = generateGraph.read_xvg(lst_files[0])
nRes = [int(i) for i in n_res_aux]
datos = numpy.zeros(shape=(len(nRes), max_steps))
#
#   Datos para graficar
#

cnt_row = 0
for i in lst_files:
    n_res_aux, rms_f, _, _, _, _ = generateGraph.read_xvg(i)
    cnt_col = 0
    for j in rms_f:
        if isfloat(j): #si se encuentra un valor no digito se pone el valor de la columan anterior
            datos[cnt_col][cnt_row] = j
        else:
            if cnt_row > 0:
                datos[cnt_col][cnt_row] = datos[cnt_col][cnt_row-1]
            else:
                datos[cnt_col][cnt_row] = 0
        cnt_col += 1
    cnt_row += 1

#
# Generar grafica
#
titulo = "Rmsd Fluctuation steps"
y_label= "ResN"
x_label = "(ps)"

generateGraph.generate_heatmap(datos,lst_steps, out_png, titulo, x_label, y_label, x_ticks, y_ticks, max_steps)
