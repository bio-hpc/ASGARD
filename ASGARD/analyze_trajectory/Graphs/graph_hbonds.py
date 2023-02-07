#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Genera grafica de los numeros de hbonds durante la simulacion
#
import sys
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()
if len(sys.argv) !=4:
    print("debe introudir")
    print("1 nombre del fichero de la graafica")
    print("2 Titulo")
    print ("file_out.png")
    exit()
file = sys.argv[1]
title = sys.argv[2]
out_png = sys.argv[3]
step, numH, _, x_title, y_title, _ = generateGraph.read_xvg(file)
generateGraph.generate_histogram([], step, numH, out_png, x_title, y_title, title)
