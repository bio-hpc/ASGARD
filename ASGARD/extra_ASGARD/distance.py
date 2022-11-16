#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys
import math
import subprocess

#
#    Comprueba la distancia entre ligando y proteina, quedandose con los residuos menores de longitudMaxima
#    Es util para gromacs
#    depende de:
#         lanzador/scriptsLanzador/standarFileCoords.py


rutaStandarFileCoors = "lanzador/scriptsLanzador/standarFileCoords.py"


def leer_file(filename, array):
    f = open(filename)
    for linea in f:
        if not linea.startswith("#") and linea.strip() != "":
            aux = linea.strip().split(":")
            v = list()
            v.append(float(aux[0]))
            v.append(float(aux[1]))
            v.append(float(aux[2]))
            v.append(aux[3])
            v.append(aux[5])
            v.append(aux[6])
            array.append(v)
    f.close()
    return array


def distancia(array_lig, array_prot):
    # if not arrayProt[4]+"_"+arrayProt[5] in residuosCerca:
    x = math.pow(array_lig[0] - array_prot[0], 2)
    y = math.pow(array_lig[1] - array_prot[1], 2)
    z = math.pow(array_lig[2] - array_prot[2], 2)
    ret = math.sqrt(x + y + z)
    return ret


def create_hash(residuo, distance, lig_atom):
    if residuo not in hashCarboHidrato:
        hashCarboHidrato[residuo] = [distance, lig_atom]
    else:
        if distance < hashCarboHidrato[residuo][0]:
            hashCarboHidrato[residuo] = [distance, lig_atom]


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("error:")
        print("debe introducir proteina y legando")
        exit()
    prot = sys.argv[1]
    lig = sys.argv[2]
    ligCoords = []
    protCoords = []
    hashCarboHidrato = {}

    longitudMaxima = 5

    comando = "python2 " + rutaStandarFileCoors + " " + prot + " >" + prot + ".tmp"
    subprocess.check_output(comando, shell=True)
    comando = "python2 " + rutaStandarFileCoors + " " + lig + " >" + lig + ".tmp"
    subprocess.check_output(comando, shell=True)
    ligCoords = leer_file(lig + ".tmp", ligCoords)
    protCoords = leer_file(prot + ".tmp", protCoords)

    #    Borramos los temporales
    os.unlink(lig + ".tmp")
    os.unlink(prot + ".tmp")

    # for i in ligCoords:
    #    print i

    for i in ligCoords:
        for j in protCoords:
            dist = distancia(i, j)
            if dist < longitudMaxima:
                create_hash(j[4] + "_" + j[5], dist, i[3])

    for res, numRepe in hashCarboHidrato.items():
        print(res + " " + str(numRepe[1]) + " " + str(numRepe[0]))
