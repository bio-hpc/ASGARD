#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#   Generate graph with the number of hydrogen bonds for the simulation
#
import sys
from GenerateGraph.GenerateGraph import GenerateGraph

generateGraph = GenerateGraph()

if len(sys.argv) != 4:
    print("You must provide:")
    print("1. Filename of the graph")
    print("2. Title")
    print("3. Output file (file_out.png)")
    exit()

file = sys.argv[1]
title = sys.argv[2]
out_png = sys.argv[3]

step, numH, _, x_title, y_title, _ = generateGraph.read_xvg(file)
print(len(numH))
# for i in range(0,len(step)):
#     print(step[i],numH[i])
generateGraph.generate_histogram([], step, numH, out_png, x_title, y_title, title)