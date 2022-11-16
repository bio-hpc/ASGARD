"""
This module is written for reading the atom and bond information from mol2
file.
"""
from __future__ import absolute_import
from pymsmtmol.mol import Atom, Residue, Molecule
from pymsmtmol.element import ionnamel, Metalpdb
import sys
import linecache

def get_atominfo(fname):

    #Detect the line numbers of each part information
    fp = open(fname, 'r')
    lnum = 1
    for line in fp:
      if (line == "@<TRIPOS>ATOM\n"):
        atbgin = lnum + 1
      elif (line == "@<TRIPOS>BOND\n"):
        atend = lnum
      lnum = lnum + 1
    fp.close()

    Atoms = {}
    Residues = {}

    atids = []
    resids = []
    resnamedict = {}
    conterdict = {}

    for i in range(atbgin, atend):
      atid, atname, crdx, crdy, crdz, atomtype, resid, resname, charge = \
      linecache.getline(fname, i).split()[:9]

      #for atom part
      gtype = "ATOM"
      atid = int(atid)
      atids.append(atid)
      crd = (float(crdx),float(crdy),float(crdz))
      charge = float(charge)
      resid = int(resid)

      if (resname, atname) in Metalpdb.keys():
        element = Metalpdb[(resname, atname)]
      else:
        element = atname[0]

      Atoms[atid] = Atom(gtype, atid, atname, element, atomtype, crd, charge, resid, resname)

      #for the residue part
      if resid not in resids:
        resids.append(resid)
      if resid not in resnamedict.keys():
        resnamedict[resid] = resname

    #clean the memory
    linecache.clearcache()

    resids.sort()

    for i in resids:
      preconter = []
      for j in atids:
        if (Atoms[j].resid == i) and (j not in preconter):
          preconter.append(j)
      preconter.sort()
      conterdict[i] = preconter

    for i in resids:
      resname = resnamedict[i]
      resconter = conterdict[i]
      Residues[i] = Residue(i, resname, resconter)

    del resnamedict
    del conterdict

    mol = Molecule(Atoms, Residues)

    return mol, atids, resids

def get_bondinfo(fname):
    fp = open(fname, 'r')
    lnum = 1
    for line in fp:
      if (line == "@<TRIPOS>BOND\n"):
        bdbgin = lnum + 1
      elif (line == "@<TRIPOS>SUBSTRUCTURE\n"):
        bdend = lnum
      lnum = lnum + 1
    fp.close()

    blist = []

    for i in range(bdbgin, bdend):
      a, b, c, d = linecache.getline(fname, i).split()[:4]

      if (int(b) > int(c)):
        blist.append((int(c), int(b), d))
      else:
        blist.append((int(b), int(c), d))

    linecache.clearcache()
    return blist

def get_pure_type(onelist):
    newlist = []
    for i in onelist:
      if (not i in newlist) & (not i[::-1] in newlist):
        newlist.append(i)
    return newlist

def get_pure_num(onelist):
    newlist = []
    for i in onelist:
      if not i in newlist:
        newlist.append(i)
    return newlist    
