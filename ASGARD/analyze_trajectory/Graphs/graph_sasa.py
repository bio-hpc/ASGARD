#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Generate SAS graph, for the whole system

import sys
from GenerateGraph.GenerateGraph import GenerateGraph

generateGraph = GenerateGraph()

if len(sys.argv) != 4:
    print("Parameters: ")
    print("file xvg sasa o")
    print("file xvg sasa odp")
    print("out png")
    exit()

file_o = sys.argv[1]
file_odp = sys.argv[2]
out_png = sys.argv[3]

lst = []
steeps = []

x, lst_total, _, _, _, _ = generateGraph.read_xvg(file_o)
lst.append(lst_total)
steeps.append(x)

x, lst_solv, _, _, _, _ = generateGraph.read_xvg(file_odp)
lst.append(lst_solv)
steeps.append(x)

x_title = "Time (ps)"
y_title = "Area"
legend = ["Total", "DG solv"]
title = "Solvent Accessible Surface"
generateGraph.line_graph(legend, steeps, lst, out_png, x_title, y_title, title, "")
