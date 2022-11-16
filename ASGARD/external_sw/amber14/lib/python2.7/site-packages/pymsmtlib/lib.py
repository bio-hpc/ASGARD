"""
This module is used for getting the parameter information from mol2 and
parm*.dat file.
"""
from pymsmtmol.readmol2 import get_atominfo
import os
import numpy
from scipy.optimize import curve_fit
import linecache

#-----------------------------------------------------------------------------

#Get the charge library file and the atom type definition for the amino acid
#atoms

#-----------------------------------------------------------------------------

def get_lib_dict(parms):

    amberhome = os.getenv('AMBERHOME')

    add = amberhome + '/AmberTools/src/pymsmt/pymsmtlib/'

    if parms in ['ff94', 'ff99', 'ff99SB']:
      mol, atids, resids = get_atominfo(add + 'parm94.mol2')
    elif (parms == 'ff03'):
      mol, atids, resids = get_atominfo(add + 'parm03.mol2')
    elif (parms == 'ff03.r1'):
      mol, atids, resids = get_atominfo(add + 'parm03_r1.mol2')
    elif (parms == 'ff10'):
      mol, atids, resids = get_atominfo(add + 'parm10.mol2')
    elif parms in ['ff12SB', 'ff14SB']:
      mol, atids, resids = get_atominfo(add + 'parm12.mol2')
    else:
      mol, atids, resids = get_atominfo(parms)

    libdict = {} #resname + atname : atom type, atom charge
    chargedict = {} #resname : charge

    for i in resids:
      charge = 0.0
      for j in mol.residues[i].resconter:
        #key is residue name and atom name
        key =  mol.residues[i].resname + '-' +  mol.atoms[j].atname
        if len(mol.atoms[j].atomtype) == 1:
          mol.atoms[j].atomtype = mol.atoms[j].atomtype + ' '
        #value is atom type and atom charge
        val =  (mol.atoms[j].atomtype, mol.atoms[j].charge)
        libdict[key] = val
        #cummulative charge of the residue
        charge = charge + mol.atoms[j].charge
      chargedict[mol.residues[i].resname] = charge

      #Alias HN as H
      atnames = [mol.atoms[k].atname for k in mol.residues[i].resconter]
      if set(['H', 'N', 'C', 'O']) < set(atnames):
        libdict[mol.residues[i].resname + '-HN'] = \
        libdict[mol.residues[i].resname + '-H']
    return libdict, chargedict

#-----------------------------------------------------------------------------
#Read the frcmod file
#-----------------------------------------------------------------------------

class Parms:
    def __init__(self, mass, bond, ang, dih, imp, nb):
        self.mass = mass
        self.bond = bond
        self.ang = ang
        self.dih = dih
        self.imp = imp
        self.nb = nb

    def combine(self, Parms2):

        #Mass
        self.mass.update(Parms2.mass)

        #Bond
        for i in Parms2.bond.keys():
          if (i in self.bond.keys()) or (i[::-1] in self.bond.keys()):
            self.bond[i] = Parms2.bond[i]
          else:
            self.bond[i] = Parms2.bond[i]

        #Angle
        for i in Parms2.ang.keys():
          if (i in self.ang.keys()) or (i[::-1] in self.ang.keys()):
            self.ang[i] = Parms2.ang[i]
          else:
            self.ang[i] = Parms2.ang[i]

        #Dih
        for i in Parms2.dih.keys():
          if (i in self.dih.keys()) or (i[::-1] in self.dih.keys()):
            self.dih[i] = Parms2.dih[i]
          else:
            self.dih[i] = Parms2.dih[i]

        #Imp
        for i in Parms2.imp.keys():
          if (i in self.imp.keys()) or ((i[0], i[3], i[2], i[1]) in self.imp.keys()) \
            or ((i[1], i[0], i[2], i[3]) in self.imp.keys()) \
            or ((i[1], i[3], i[2], i[0]) in self.imp.keys()) \
            or ((i[3], i[0], i[2], i[1]) in self.imp.keys()) \
            or ((i[3], i[1], i[2], i[0]) in self.imp.keys()):
            self.imp[i] = Parms2.imp[i]
          else:
            self.imp[i] = Parms2.imp[i]

        #Nb
        self.nb.update(Parms2.nb)

#-----------------Read parameters from dat or frcmod file----------------------

def readmass(massparms, line):
    if len(line) > 2:
      attyp = line[0:2]
      if attyp not in ['\n', '']:
        mass = line[2:].strip('\n')
        massparms[attyp] = mass
    return massparms

def readbond(bondparms, line):
    if len(line) > 5:
      if line[2] == '-':
        at1 = line[0:2]
        at2 = line[3:5]
        bondparm = line[5:].strip('\n')
        bondparms[(at1, at2)] = bondparm
    return bondparms

def readang(angparms, line):
    if len(line) > 8:
     if line[2] == '-' and line[5] == '-':
       at1 = line[0:2]
       at2 = line[3:5]
       at3 = line[6:8]
       angparm = line[8:].strip('\n')
       angparms[(at1, at2, at3)] = angparm
    return angparms

def readdih(dihparms, line):
    if len(line) > 11:
      if line[2] == '-' and line[5] == '-' and line[8] == '-':
        at1 = line[0:2]
        at2 = line[3:5]
        at3 = line[6:8]
        at4 = line[9:11]
        dihtyp = (at1, at2, at3, at4)
        nvnp = line[11:43]
        pero = line[43:52].strip('\n')
        annot = line[52:].strip('\n')
        dihparm = [nvnp, pero, annot]
        if dihtyp in dihparms.keys():
          has_pero = dihparms[dihtyp][1::3]
          has_pero = [abs(int(i.strip().strip('.'))) for i in has_pero]
          if abs(int(dihparm[1].strip().strip('.'))) not in has_pero:
            dihparm = dihparms[dihtyp] + dihparm
        elif dihtyp[::-1] in dihparms.keys():
          has_pero = dihparms[dihtyp[::-1]][1::3]
          has_pero = [abs(int(i.strip().strip('.'))) for i in has_pero]
          if abs(int(dihparm[1].strip().strip('.'))) not in has_pero:
            dihparm = dihparms[dihtyp[::-1]] + dihparm
        dihparms[dihtyp] = dihparm
    return dihparms

def readgaffdih(dihparms, line):
    if len(line) > 11:
      if line[2] == '-' and line[5] == '-' and line[8] == '-':
        at1 = line[0:2]
        at2 = line[3:5]
        at3 = line[6:8]
        at4 = line[9:11]
        dihtyp = (at1, at2, at3, at4)
        nvnp = line[11:43]
        pero = line[43:52].strip('\n')
        annot = line[52:].strip('\n')
        dihparm = [nvnp, pero, annot]
        if dihtyp in dihparms.keys():
          has_pero = dihparms[dihtyp][1::3]
          has_pero = [abs(int(i.strip().strip('0').strip('-').strip('.'))) for i in has_pero]
          if abs(int(dihparm[1].strip().strip('0').strip('-').strip('.'))) not in has_pero:
            dihparm = dihparms[dihtyp] + dihparm
        elif dihtyp[::-1] in dihparms.keys():
          has_pero = dihparms[dihtyp[::-1]][1::3]
          has_pero = [abs(int(i.strip().strip('0').strip('-').strip('.'))) for i in has_pero]
          if abs(int(dihparm[1].strip().strip('0').strip('-').strip('.'))) not in has_pero:
            dihparm = dihparms[dihtyp[::-1]] + dihparm
        dihparms[dihtyp] = dihparm
    return dihparms

def readimp(impparms, line):
    if len(line) > 11:
      if line[2] == '-' and line[5] == '-' and line[8] == '-':
        impparm = line[11:].strip('\n')
        at1 = line[0:2]
        at2 = line[3:5]
        at3 = line[6:8]
        at4 = line[9:11]
        impparms[(at1, at2, at3, at4)] = impparm
    return impparms

def readeqnb(eqdict, line):
    eqatms = line.split()
    for i in eqatms:
      if len(i) == 1:
        i = i + ' '
    eqdict[eqatms[0]] = eqatms[1:]
    return eqdict

def readnb(nbparms, line):
    if len(line) > 2:
      at1 = line.split()[0]
      if at1 not in ['\n', '']:
        if len(at1) == 1:
          at1 = at1 + ' '
        nbparm = line[4:].strip('\n')
        nbparms[at1] = nbparm
    return nbparms

def read_frcmod_file(frcmodf):

    #Get range of each parameter part in the frcmodf
    rfrcmodf = open(frcmodf, 'r')
    cardlist = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB']
    lnlist1 = []
    lnlist2 = []
    ln = 1
    for line in rfrcmodf:
      for card in cardlist:
        if line[0:len(card)] == card:
          lnlist1.append(card)
          lnlist2.append(ln)
      ln = ln + 1
    tln = ln - 1
    rfrcmodf.close()

    lndict = {}
    for i in range(0, len(lnlist1)-1):
      lndict[lnlist1[i]] = (lnlist2[i]+1, lnlist2[i+1])

    lndict[lnlist1[-1]] = (lnlist2[-1]+1, tln)

    #Read the parameter into dicts
    massparms = {}
    bondparms = {}
    angparms = {}
    dihparms = {}
    impparms = {}
    nbparms = {}

    for i in lndict.keys():
      if i == "MASS":
        for j in range(lndict[i][0],lndict[i][1]):
          line = linecache.getline(frcmodf, j)
          massparms = readmass(massparms, line)
      elif i == "BOND":
        for j in range(lndict[i][0], lndict[i][1]):
          line = linecache.getline(frcmodf, j)
          bondparms = readbond(bondparms, line)
      elif i == "ANGL":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(frcmodf, j)
           angparms = readang(angparms, line)
      elif i == "DIHE":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(frcmodf, j)
           dihparms = readdih(dihparms, line)
      elif i == "IMPR":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(frcmodf, j)
           impparms = readimp(impparms, line)
      elif i == "NONB":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(frcmodf, j)
           nbparms = readnb(nbparms, line)

    linecache.clearcache()
    parmdict = Parms(massparms, bondparms, angparms, dihparms, impparms, nbparms)

    return parmdict

def get_parm_dict(ffchoice, gaff, frcmodfs):

    #-------------------------------------------------------------------------
    #1. Read the parm*.dat file
    #-------------------------------------------------------------------------
    add = os.getenv('AMBERHOME') + '/dat/leap/parm/'

    if ffchoice == 'ff94': #parm94
      parmf = add + 'parm94.dat'
      lndict = {'MASS': (2, 57), 'BOND': (59, 142), 'ANGL': (143, 334),
                'DIHE': (335, 416), 'IMPR': (417, 448), 'EQUA': (451, 453),
                'NONB': (455, 489)}
    elif ffchoice in ['ff99', 'ff99SB', 'ff03', 'ff03.r1']: 
      parmf = add + 'parm99.dat'
      lndict = {'MASS': (2, 66), 'BOND': (68, 184), 'ANGL': (185, 466),
                'DIHE': (467, 631), 'IMPR': (632, 670), 'EQUA': (673, 675),
                'NONB': (677, 719)}
    elif ffchoice in ['ff10', 'ff12SB', 'ff14SB']:
      parmf = add + 'parm10.dat'
      lndict = {'MASS': (2, 65), 'BOND': (67, 218), 'ANGL': (219, 619),
                'DIHE': (620, 895), 'IMPR': (896, 955), 'EQUA': (958, 960),
                'NONB': (962, 1001)}

    #define the parameter dicts
    massparms = {}
    bondparms = {}
    angparms = {}
    dihparms = {}
    impparms = {}
    eqdict = {}
    nbparms = {}

    for i in lndict.keys():
      if i == "MASS":
        for j in range(lndict[i][0],lndict[i][1]):
          line = linecache.getline(parmf, j)
          massparms = readmass(massparms, line)
      elif i == "BOND":
        for j in range(lndict[i][0], lndict[i][1]):
          line = linecache.getline(parmf, j)
          bondparms = readbond(bondparms, line)
      elif i == "ANGL":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(parmf, j)
           angparms = readang(angparms, line)
      elif i == "DIHE":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(parmf, j)
           dihparms = readdih(dihparms, line)
      elif i == "IMPR":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(parmf, j)
           impparms = readimp(impparms, line)
      elif i == "EQUA":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(parmf, j)
           eqdict = readeqnb(eqdict, line)
      elif i == "NONB":
         for j in range(lndict[i][0], lndict[i][1]):
           line = linecache.getline(parmf, j)
           nbparms = readnb(nbparms, line)

    linecache.clearcache()

    #Deal with the equil atoms
    for i in eqdict.keys():
      for j in eqdict[i]:
        if len(i) == 1:
          nbparms[j] = nbparms[i + ' ']
        else:
          nbparms[j] = nbparms[i]

    parmdict = Parms(massparms, bondparms, angparms, dihparms, impparms, nbparms)

    #-------------------------------------------------------------------------
    #2. Read the frcmod file for each force field
    #-------------------------------------------------------------------------
    if ffchoice in ['ff03', 'ff03.r1', 'ff99SB', 'ff12SB', 'ff14SB']:
      if ffchoice in ['ff03', 'ff03.r1']: #Year: 2003
        parmf1 = add + 'frcmod.ff03'
      elif ffchoice == 'ff99SB': #Year: 2006
        parmf1 = add + 'frcmod.ff99SB'
      elif ffchoice == 'ff12SB': #Year: 2012
        parmf1 = add + 'frcmod.ff12SB'
      elif ffchoice == 'ff14SB': #Year: 2014
        parmf1 = add + 'frcmod.ff14SB'
      parmdict1 = read_frcmod_file(parmf1)
      parmdict.combine(parmdict1)

    #-------------------------------------------------------------------------
    #3. GAFF
    #-------------------------------------------------------------------------
    if gaff == 1:

      parmf2 = add + 'gaff.dat'

      lndict2 = {'MASS': (2, 73), 'BOND': (75, 882), 'ANGL': (883, 5131),
                 'DIHE': (5132, 5848), 'IMPR': (5849, 5887),
                 'NONB': (5892, 5963)}

      massparms2 = {}
      bondparms2 = {}
      angparms2 = {}
      dihparms2 = {}
      impparms2 = {}
      nbparms2 = {}

      for i in lndict2.keys():
        if i == "MASS":
          for j in range(lndict2[i][0], lndict2[i][1]):
            line = linecache.getline(parmf2, j)
            massparms2 = readmass(massparms2, line)
        elif i == "BOND":
          for j in range(lndict2[i][0], lndict2[i][1]):
            line = linecache.getline(parmf2, j)
            bondparms2 = readbond(bondparms2, line)
        elif i == "ANGL":
           for j in range(lndict2[i][0], lndict2[i][1]):
             line = linecache.getline(parmf2, j)
             angparms2 = readang(angparms2, line)
        elif i == "DIHE":
           for j in range(lndict2[i][0], lndict2[i][1]):
             line = linecache.getline(parmf2, j)
             dihparms2 = readgaffdih(dihparms2, line)
        elif i == "IMPR":
           for j in range(lndict2[i][0], lndict2[i][1]):
             line = linecache.getline(parmf2, j)
             impparms2 = readimp(impparms2, line)
        elif i == "NONB":
           for j in range(lndict2[i][0], lndict2[i][1]):
             line = linecache.getline(parmf2, j)
             nbparms2 = readnb(nbparms2, line)

      linecache.clearcache()
      parmdict2 = Parms(massparms2, bondparms2, angparms2, dihparms2,
                        impparms2, nbparms2)
      parmdict.combine(parmdict2)

    #-------------------------------------------------------------------------
    #4. Additional frcmod file
    #-------------------------------------------------------------------------

    for i in frcmodfs:
      parmdict3 = read_frcmod_file(i)
      parmdict.combine(parmdict3)

    return parmdict

def expf(x, a, b, c):
    return a * numpy.exp(-b * x) + c

def getfc(fname, dis):
    amberhome = os.getenv('AMBERHOME')
    add = amberhome + '/AmberTools/src/pymsmt/pymsmtlib/' + fname

    lengthl = []
    fcl = []
    fcf = open(add, 'r')
    for line in fcf:
      length, fc = line.split()[:2]
      length = float(length)
      lengthl.append(length)
      fc = float(fc)
      fcl.append(fc)
    fcf.close()

    initguess = [5E6, 5.5, 0.0]
    lengthl = numpy.array(lengthl)
    fcl = numpy.array(fcl)

    optparas, convar = curve_fit(expf, lengthl, fcl, p0=initguess, maxfev=10000)
    a, b, c = optparas

    val = expf(dis, a, b, c)
    val = round(val, 1)
    return val

