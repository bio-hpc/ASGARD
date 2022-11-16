#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

#
#    Autor:      Jorge De La Peña Garcia (jorge.dlpg@gmail.com)
#    version: 1
#  		 Genera grafica Gyrate
#
import commands
import sys
from GenerateGraph.GenerateGraph import GenerateGraph
generateGraph = GenerateGraph()

if len(sys.argv) != 3:
	print("debe introudir")
	print("xvg out ")
	print("png out ")
	exit()
fichero = sys.argv[1]
out_png = sys.argv[2]
x, y, title, x_title, y_title, _ = generateGraph.read_xvg(fichero)
#
#	Titurlos que pone en el fichero
#
comando = 'cat '+fichero+' |grep s0 |awk -F\\"  \'{print $2}\''
legend = []
legend.append(commands.getoutput(comando))
generateGraph.line_graph(legend, x, y, out_png, x_title, y_title, title, "")
