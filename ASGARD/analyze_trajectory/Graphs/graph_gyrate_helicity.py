#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#     version: 2
#  		 Generate Gyrate graph
#
import subprocess
import sys
from GenerateGraph.GenerateGraph import GenerateGraph

generateGraph = GenerateGraph()

if len(sys.argv) != 3:
    print("You must provide")
    print("xvg out ")
    print("png out ")
    sys.exit()

fichero = sys.argv[1]
out_png = sys.argv[2]
x, y, title, x_title, y_title, _ = generateGraph.read_xvg(fichero)

comando = 'cat '+fichero+' |grep s0 |awk -F\\"  \'{print $2}\''

legend = []
legend.append(subprocess.getoutput(comando))
generateGraph.line_graph(legend, x, y, out_png, x_title, y_title, title, "")

if "distance" in str(fichero) and "fluct" in str(fichero):
    x_title = "Frequency"
    y_title = "Distance (nm)"
    generateGraph.line_graph(legend, x, y, out_png, x_title, y_title, title, "")
