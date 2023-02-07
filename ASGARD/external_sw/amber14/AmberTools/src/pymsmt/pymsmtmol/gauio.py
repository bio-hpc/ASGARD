"""
Module for writting a Gaussian file and read the coordinates and force
constants from Gaussian output file.
"""
import numpy

def write_gauatm(gauatm, fname):
    wf = open(fname, 'a')
    print >> wf, "%-6s   %8.3f%8.3f%8.3f" %(gauatm.element, \
        gauatm.crdx, gauatm.crdy, gauatm.crdz)
    wf.close()

def write_gauatm_opth(gauatm, fname):
    wf = open(fname, 'a')
    if gauatm.element == "H":
      print >> wf, "%-6s  0 %8.3f%8.3f%8.3f" %(gauatm.element, \
          gauatm.crdx, gauatm.crdy, gauatm.crdz)
    else:
      print >> wf, "%-6s -1 %8.3f%8.3f%8.3f" %(gauatm.element, \
          gauatm.crdx, gauatm.crdy, gauatm.crdz)
    wf.close()

def get_crds_from_fchk(fname, bstring, estring):

    crds = []

    bl = len(bstring)
    el = len(estring)

    fp = open(fname, 'r')
    i = 1
    for line in fp:
      if line[0:bl] == bstring:
        beginl = i
      elif line[0:el] == estring:
        endl = i
      i = i + 1
    fp.close()

    fp1 = open(fname, 'r')
    i = 1
    for line in fp1:
      if (i > beginl) and (i < endl):
        crd = line.split()
        for j in crd:
          if j != ' ':
            crds.append(float(j))
      i = i + 1
    fp1.close()

    return crds

def get_matrix_from_fchk(fname, bstring, estring, msize):

    crds = []

    bl = len(bstring)
    el = len(estring)

    fp = open(fname, 'r')
    i = 1
    for line in fp:
      if line[0:bl] == bstring:
        beginl = i
      elif line[0:el] == estring:
        endl = i
      i = i + 1
    fp.close()

    fcmatrix = numpy.array([[float(0) for x in range(msize)] for x in range(msize)])

    fp1 = open(fname, 'r')
    i = 1
    for line in fp1:
      if (i > beginl) and (i < endl):
        crd = line.split()
        for j in crd:
          if j != ' ':
            crds.append(float(j))
      i = i + 1
    fp1.close()

    i = 0
    for j in range(0, msize):
      for k in range(0, j+1):
        fcmatrix[j][k] = crds[i]
        fcmatrix[k][j] = crds[i]
        i = i + 1

    return fcmatrix

def get_fc_from_log(logfname):

    stringle = 'calculate D2E/DX2 analytically'
    stringfc = ' Internal force constants:'

    sturefs = []
    vals = []

    ##Get the values for each bond, angle and dihedral
    fp = open(logfname, 'r')
    for line in fp:
      if stringle in line:
        line = line.strip('\n')
        line = line.strip('!')
        line = line.lstrip(' !') 
        line = line.split()
        val = float(line[2])
        vals.append(val)
        typ = line[0][0]
        line[1] = line[1].strip(typ)
        line[1] = line[1].lstrip('(')
        line[1] = line[1].strip(')')
        ats = line[1].split(',')
        ats = [int(i) for i in ats]
        sturefs.append(tuple(ats))
    fp.close()

    maxnum = len(sturefs)

    ##Get the force constant begin line
    fp = open(logfname, 'r')
    lnum = 1
    for line in fp:
      if stringfc in line:
        blnum = lnum + 2
      lnum = lnum + 1
    fp.close()

    blnum1 = blnum

    ##Get the line number list to read the force constant
    numl = []

    if (maxnum%5 != 0):
      for i in range(0, int(maxnum/5)+1):
        gap = maxnum - len(numl) + 1
        for j in range(0, 5):
          if len(numl) < maxnum:
            numl.append(blnum1 + j)
        blnum1 = blnum1 + gap
    else:
      for i in range(0, int(maxnum/5)):
        gap = maxnum - len(numl) + 1
        for j in range(0, 5):
          if len(numl) < maxnum:
            numl.append(blnum1 + j)
        blnum1 = blnum1 + gap

    fcs = []
    fp = open(logfname, 'r')
    lnum = 1
    for line in fp:
      if lnum in numl:
        line = line.strip('\n')
        line = line.split()
        numv = len(line)
        fc = float(line[-1].replace('D', 'e'))
        fcs.append(fc)
      lnum = lnum + 1
    fp.close()

    ##Return three lists: identifications, values, force constants
    return sturefs, vals, fcs

def get_crds_from_log(logfname, g0x):

    if g0x == 'g03':
      fp = open(logfname)
      ln = 1
      for line in fp:
        if 'Redundant internal coordinates' in line:
          bln = ln + 3
        elif 'Recover connectivity data from disk' in line:
          eln = ln - 1
        ln = ln + 1
      fp.close()
    elif g0x == 'g09':
      fp = open(logfname)
      ln = 1
      for line in fp:
        if 'Redundant internal coordinates' in line:
          bln = ln + 1
        elif 'Recover connectivity data from disk' in line:
          eln = ln - 1
        ln = ln + 1
      fp.close()

    crds = []
    ln = 1
    fp1 = open(logfname)
    for line in fp1:
      if (ln >= bln) and (ln <= eln):
        line = line.strip('\n')
        line = line.split(',')
        line = line[-3:]
        line = [float(i) for i in line]
        crds += line
      ln = ln + 1
    fp1.close()

    return crds

