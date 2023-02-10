#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    version:    1
#   Generate a standard graph for a xvg file
#
import sys
import os
ALLOW_EXTENSION = ".xvg"

from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()
if len(sys.argv) < 4:
    print("Parameters")
    print("1 file xvg")
    print("2 file xvg")
    print("3 ...")
    print("4 ...")
    print("n-1 Title")
    print("n Out png")
    print("Introduce 1 or more xvg files, and the title and output file")
    exit()
#file_xvg = sys.argv[1]
title = sys.argv[len(sys.argv)-2]
out_png = sys.argv[len(sys.argv)-1]
lst_x = []
lst_y = []
for i in range(1, len(sys.argv)-2):
    file_xvg = sys.argv[i]
    if os.path.splitext(file_xvg)[1] == ALLOW_EXTENSION:
        x, y, _, x_title, y_title, _ = generateGraph.read_xvg(file_xvg)
        lst_x.append(x)
        lst_y.append(y)
    else:
        print ("Error file: " + file_xvg)
        exit()
if len(lst_x) == 1:
    lst_x = lst_x[0]
    lst_y = lst_y[0]
generateGraph.graph_doble_line("", lst_x, lst_y, out_png, x_title, y_title, title, "")



