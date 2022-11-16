stribución de personas por especialidade#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys
import math
import subprocess

#
#    Mide la distnacia maxima del ligando busca x1, y1, z1 maximo y x1, y1, z1 minimo y calcula la diagonal
#    depende de:
#        standarFileCoords.py
#
StandarCoords = "lanzador/scriptsLanzador/standarFileCoords.py"


def distance(_a, _b):
    _x = math.pow((_a[0] - _b[0]), 2)
    _y = math.pow((_a[1] - _b[1]), 2)
    _z = math.pow((_a[2] - _b[2]), 2)
    suma = _x + _y + _z
    return round(math.sqrt(suma), 3)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("debe introducir un fichero de ligando")
        exit()
    fichero = sys.argv[1]

    comando = "python " + StandarCoords + " " + fichero + " | grep -v \# >" + fichero + "TMP"
    lineas = subprocess.check_output(comando, shell=True)

    x = []
    y = []
    z = []

    f = open(fichero+"TMP")
    for lineas in f:
        aux = lineas.split(":")
        x.append(float(aux[0]))
        y.append(float(aux[1]))
        z.append(float(aux[2]))
    f.close()
    os.unlink(fichero+"TMP")

    distanciaTotal = 0
    for i in range(len(x)):
        a = [x[i], y[i], z[i]]
        for j in range(len(x)):
            b = [x[j], y[j], z[j]]
            d = distance(a, b)
            if d > distanciaTotal:
                distanciaTotal = d

    print(distanciaTotal)
