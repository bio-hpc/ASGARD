#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
    OJO: Script totalmente experimental para incorporar BD flexible.
"""

import os
import sys
import math
import argparse

import numpy as np

package_path = os.path.join(os.path.dirname(__file__), './../../externalSw/python2-packages')
sys.path.append(os.path.realpath(package_path))

import pybel

# Bloque "prestado" de Get_histogram.Tools ########################################################################


def read_molecule(fname):
    # pybel requiere especificar el tipo de archivo que se va a leer
    # así que tomamos la extensión para averiguarlo
    file_name, file_ext = os.path.splitext(fname)
    mols = [mol for mol in pybel.readfile(file_ext.strip('.'), fname)]
    assert len(mols) == 1
    return mols[0]


def get_sq_dist(point1, point2):
    return sum((point1 - point2)**2)


def get_dist(point1, point2):
    return math.sqrt(get_sq_dist(point1, point2))


###################################################################################################################

"""
  Estrategia para determinar qué residuos hay que hacer flexibles:

  - Calcular distancia desde el Carbono alfa a todos los demás residuos.
    La distancia se define respecto al átomo no-hidrógeno más cercano.

"""


def get_residue_ac(_residue):
    """
    Encontrar el carbono alfa de un residuo
    """

    for _atom in pybel.ob.OBResidueAtomIter(_residue):
        if _residue.GetAtomID(_atom).strip() == 'CA':
            return _atom

    return None


def get_residue_by_id(_id, _mol):
    for _residue in pybel.ob.OBResidueIter(_mol):
        if _residue.GetNum() == _id:
            return _residue
    return None


def get_residue_neighbors(_residue, _dist, _mol):
    """
    Marcar los residuos que están a la distancia indicada o menos de cada carbono alfa.

    OJO: la matriz resultante no es simétrica:
    1a. coordenada es origen (CA), la 2a. coord es destino (cualquier átomo del residuo)
    """

    _ca = get_residue_ac(_residue)

    if not hasattr(_ca, 'coords'):
        _ca.coords = np.array((_ca.GetX(), _ca.GetY(), _ca.GetZ()))

    _neighbors = []
    for _other in pybel.ob.OBResidueIter(_mol):
        for _atom in pybel.ob.OBResidueAtomIter(_other):

            if _atom.IsHydrogen():
                continue

            if not hasattr(_atom, 'coords'):
                _atom.coords = np.array((_atom.GetX(), _atom.GetY(), _atom.GetZ()))

            d = get_dist(_atom.coords, _ca.coords)

            if d <= _dist:
                _neighbors.append(_other)
                break

    return _neighbors


if __name__ == '__main__':

    ap = argparse.ArgumentParser(description='Generar listado de residuos a flexibilizar para hacer BD flexible.')
    ap.add_argument('prot', help='Archivo de la proteína a examinar.')
    ap.add_argument('res_num', type=int, help='Residuo origen para el que obtener indicador de flexibilidad.')
    ap.add_argument('-d', dest='dist', type=float, default=5., help='Radio en el que buscar residuos cercanos.')

    args = ap.parse_args()

    pybel.ob.obErrorLog.SetOutputLevel(pybel.ob.obError)

    pdb_code, _ = os.path.splitext(os.path.basename(args.prot))

    mol = read_molecule(args.prot)
    mol = mol.OBMol

    # Obtener el residuo
    residue = get_residue_by_id(args.res_num, mol)

    # Debemos obtener un residuo
    assert residue

    # Encontrar vecinos
    neighbors = get_residue_neighbors(residue, args.dist, mol)

    # Debe haber al menos un residuo (el de entrada)
    assert neighbors

    # Crear la lista
    print(
        ','.join(
            '{}:{}:{}{}'.format(pdb_code, neighbor.GetChain().strip(),
                                neighbor.GetName().strip(), neighbor.GetNum())
            for neighbor in neighbors
        )
    )

