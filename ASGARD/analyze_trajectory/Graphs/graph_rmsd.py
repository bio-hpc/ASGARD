#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    version:    1
#
#   Genera grafica de rmsd o rmsdf para gromcas y su distribucion
#
import sys
import os
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()
ALLOW_EXTENSION = ".xvg"
TOKEN_DISTRIBUTION = [ "rmsd_distribution", "distance_distribution" ]

TOKEN = ["rmsd", 'distance']

if len(sys.argv) < 4:
    print("Parameters")
    print("1 file xvg")
    print("2 file xvg")
    print("3 ...")
    print("4 ...")
    print("n-2 Title")
    print("n-1 Out png")
    print("n distribution true,  flase)")
    print("Se le pasa 1 o varios ficheros xml y al finalizar el titulo y el output")
    exit()
title = sys.argv[len(sys.argv)-2]
out_png = sys.argv[len(sys.argv)-1]
lst_x = []
lst_y = []
legend = []
for i in range(1, len(sys.argv)-2):
    file_xvg = sys.argv[i]
    if os.path.splitext(file_xvg)[1] == ALLOW_EXTENSION:
        if TOKEN[0] in file_xvg:
            aux_title = file_xvg[:file_xvg.find(TOKEN[0])-1]
        if TOKEN[1] in file_xvg:
            aux_title = file_xvg[:file_xvg.find(TOKEN[1]) - 1]
        legend.append( aux_title[aux_title.rindex("_")+1:])
        x, y, _, x_title, y_title, _ = generateGraph.read_xvg(file_xvg)
        lst_x.append(x)
        lst_y.append(y)
    else:
        print("Error file: " + file_xvg)
        exit()
if [ True for i in TOKEN_DISTRIBUTION  if i in sys.argv[1] ]:
    x_title=str("Frequency")
    y_title=str("RMSD (nm)")
    generateGraph.graph_doble_line(legend, lst_y, lst_x, out_png,x_title, y_title, title, "")
else:
    generateGraph.graph_doble_line(legend, lst_x, lst_y, out_png, x_title, y_title, title, "")
#if [ True for i in TOKEN_DISTRIBUTION  if i in sys.argv[1] ]:
#    generateGraph.graph_doble_line(legend, lst_y, lst_x, out_png,x_title, y_title, title, "")
#else:
#    generateGraph.graph_doble_line(legend, lst_x, lst_y, out_png, x_title, y_title, title, "")