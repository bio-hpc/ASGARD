#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Generate graph with the number of hydrogen bonds for the simulation
#
import sys
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()
if len(sys.argv) !=4:
    print("debe introudir")
    print("1 filename of the graph")
    print("2 Tittle")
    print ("file_out.png")
    exit()
file = sys.argv[1]
title = sys.argv[2]
out_png = sys.argv[3]
step, numH, _, x_title, y_title, _ = generateGraph.read_xvg(file)
generateGraph.generate_histogram([], step, numH, out_png, x_title, y_title, title)
