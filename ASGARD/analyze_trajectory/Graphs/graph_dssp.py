#!/usr/bin/env pythoncat
# -*- coding: utf-8 -*-
#
#	Genera grafica DSSP
#
import sys
import re   # expresion regular
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()
if len(sys.argv) != 3:
    print("Parameters: ")
    print("1º file xvg dssp")
    print("2º out put png")
    exit()
file = sys.argv[1]
out_png = sys.argv[2]
#
#   Leyenda de la proteina
#   
x, y, title, x_title, y_title, subtitle = generateGraph.read_xvg(file)
legend = []
legend.append("Structure")           #s0
legend.append("Coil")                #s1
legend.append("B-Sheet")             #s2
legend.append("Bend")                #s3
legend.append("Turn")                #s4
legend.append("3-Helix")             #s5
legend.append("Chain_Separator")     #s6 

time = []
datos = []
for i in range(6):
    datos.append( [] )
f = open(file)
for i in f:
    if not i.startswith("#") and not i.startswith("@"):
        i = re.sub(' +',' ',i).strip() #eliminamos espacios dobles inicial y final
        aux = i.split(" ")#2 num Atom 3Hydrophobic Hydrophobic Total D Gsolv
        time.append(float(aux[0]))
        datos[0].append(float(aux[1]))
        datos[1].append(float(aux[2]))
        datos[2].append(float(aux[3]))
        datos[3].append(float(aux[4]))
        if len(aux) > 5:
            datos[4].append(float(aux[5]))
        if len(aux) > 6:
            datos[5].append(float(aux[6]))
f.close()
lst_remove = []
for i in range(len(datos)-1,0,-1): #por si no llegan a 6
    if len(datos[i]) == 0:
        datos.remove(datos[i])
generateGraph.line_graph(legend, time, datos, out_png, x_title, y_title, title, "")
