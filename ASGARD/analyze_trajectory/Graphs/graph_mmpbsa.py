#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
#   version: 	 1
#   Generate MM-PBSA graph with the Vdw, elec and total energies
#
import sys
import re
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()

def check_data(data):
    try:
        if data[0].isdigit() or data[0] == "-":
            return float(data)
        else:
            return 0
    except:
        return 0

if len(sys.argv) != 4:
    print("debe introudir")
    print("1 nombre del fichero de la graafica")
    print("titulo")
    print("fie_out")
    exit()
file = sys.argv[1]
title = sys.argv[2]
out_png = sys.argv[3]
_, _, _, x_title, y_title, _ = generateGraph.read_xvg(file)
time = []
lst = []
e_van = []
e_elec = []
e_total = []
f = open(file)
for i in f:
    if not i.startswith("#") and not i.startswith("@"):
        i = re.sub(' +', ' ', i).strip().split(" ")
        time.append(check_data(i[0]))
        e_van.append(check_data(i[5]))
        e_elec.append(check_data(i[6]))
        e_total.append(check_data(i[7]))

f.close()
lst.append(e_van)
lst.append(e_elec)
lst.append(e_total)

legend = ["Vdw", "Elec", "Total"]
generateGraph.line_graph(legend, time, lst, out_png, x_title, y_title, title, "")
