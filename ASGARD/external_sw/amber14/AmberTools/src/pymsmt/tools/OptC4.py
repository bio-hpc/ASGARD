#!/usr/bin/env python
# Filename: OptC4.py
"""
This is the OptC4.py program, written by Pengfei Li in Merz Research Group,
Michigan State University. All Rights Reserved. It needs the OpenMM and
ParmEd installed and also needs the SciPy module. It was designed to optimize
the C4 terms for the metal complex in the protein system. The program is not
gurantee to work due to bug may exist. Suggestions and bug reports are welcome
to send to Pengfei Li (Email address: ldsoar1990@gmail.com).

Please cite the following paper if you use the software to perform the
modeling:

The 12-6-4 LJ-type potential:
** P. Li, K. M. Merz, JCTC, 2014, 10, 289-297
"""
from __future__ import division

# OpenMM Imports
import simtk.openmm as mm  #about force field
import simtk.openmm.app as app #about algorithm

# ParmEd imports
from chemistry import unit as u
from chemistry.amber import AmberParm, AmberMask
from chemistry.openmm import StateDataReporter, NetCDFReporter, RestartReporter

# pyMSMT Imports
from pymsmtmol.cal import calc_bond
from pymsmtmol.element import Atnum, CoRadiiDict
from interface.AmberParm import read_amber_prm

# Other Imports
from optparse import OptionParser
import os
import sys

#-----------------------------------------------------------------------------#
# Functions
#-----------------------------------------------------------------------------#

def get_typ_dict(typinds, typs):

    typdict = {}
    for i in range(0, len(typinds)):
      if typinds[i] not in typdict.keys():
        typdict[typinds[i]] = [typs[i]]
      elif typs[i] in typdict[typinds[i]]:
        continue
      else:
        typdict[typinds[i]].append(typs[i])
    return typdict

def get_rmsd(initparas):

    global idxs

    print initparas

    #Modify the C4 terms in the prmtop file
    for i in range(0, len(idxs)):
      prmtop.parm_data['LENNARD_JONES_CCOEF'][idxs[i]] = initparas[i]

    #Overwrite the prmtop file
    prmtop.write_parm('OptC4.top')

    #Perform the OpenMM optimization
    #Use AmberParm function to transfer the topology and
    #coordinate file to the object OpenMM can use
    Ambermol = AmberParm('OptC4.top', options.cfile)

    # Create the OpenMM system
    print 'Creating OpenMM System'
    if options.simupha == 'gas':
      system = Ambermol.createSystem(nonbondedMethod=app.NoCutoff)
    elif options.simupha == 'liquid':
      system = Ambermol.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=8.0*u.angstroms,
                                     constraints=app.HBonds,)

    # Create the integrator to do Langevin dynamics
    # Temperature of heat bath, Friction coefficient, Time step
    integrator = mm.LangevinIntegrator(300*u.kelvin, 1.0/u.picoseconds,
                                       1.0*u.femtoseconds,)

    # Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not
    # specify the platform to use the default (fastest) platform
    # Create the Simulation object
    if options.platf == 'ref':
       platform = mm.Platform.getPlatformByName('Reference')
       sim = app.Simulation(Ambermol.topology, system, integrator, platform)
    if options.platf == 'cpu':
       platform = mm.Platform.getPlatformByName('CPU')
       sim = app.Simulation(Ambermol.topology, system, integrator, platform)
    elif options.platf == 'cuda':
       platform = mm.Platform.getPlatformByName('CUDA')
       prop = dict(CudaPrecision='mixed')
       sim = app.Simulation(Ambermol.topology, system, integrator, platform,
                            prop)
    elif options.platf == 'opencl':
       platform = mm.Platform.getPlatformByName('OpenCL')
       prop = dict(CudaPrecision='mixed')
       sim = app.Simulation(Ambermol.topology, system, integrator, platform,
                            prop)

    # Set the particle positions
    sim.context.setPositions(Ambermol.positions)

    # Output the rst file
    restrt = RestartReporter(options.rfile, 100, write_velocities=False)
    sim.reporters.append(restrt)

    # Minimize the energy
    print 'Minimizing energy'
    sim.minimizeEnergy(maxIterations=options.maxsteps)

    # Overwrite the final file
    state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    restrt.report(sim, state)

    # Perform the RMSD calcualtion, using ptraj in AmberTools
    os.system("cpptraj -p OptC4.top -i OptC4_ptraj.in > OptC4_ptraj.out")

    ptrajof = open('OptC4_rmsd.txt', 'r')
    ln = 1
    for line in ptrajof:
      if ln == 3:
        rmsd = float(line[12:21])
      ln += 1
    ptrajof.close()

    print 'RMSD is: ', rmsd
    return rmsd

#-----------------------------------------------------------------------------#
# Main Program
#-----------------------------------------------------------------------------#

parser = OptionParser("usage: -m amber_mask -p topology_file "
                      "-c coordinate_file -r restart_file \n"
                      "       [--maxsteps maxsteps] "
                      "[--phase simulation_phase] \n"
                      "       [--size optimization_step_size] "
                      "[--method optimization_method] \n"
                      "       [--platform device_platform] "
                      "[--model metal_complex_model]")

parser.set_defaults(simupha='gas', maxsteps=1000, stepsize=10.0, minm='bfgs',
                    platf='cpu', model=1)

parser.add_option("-m", dest="ion_mask", type='string', help="Amber mask of "
                  "the center metal ion")
parser.add_option("-p", dest="pfile", type='string', help="Topology file")
parser.add_option("-c", dest="cfile", type='string', help="Coordinate file")
parser.add_option("-r", dest="rfile", type='string', help="Restart file")
parser.add_option("--maxsteps", dest="maxsteps", type='int', \
                  help="Maximum minimization steps performed by OpenMM "
                       "in each parameter optimization cycle. "
                       "[Default: 1000]")
parser.add_option("--phase", dest="simupha", type='string', \
                  help="Simulation phase, either gas or liquid. "
                       "[Default: gas]")
parser.add_option("--size", dest="stepsize", type='float', \
                  help="Step size chosen by the user for the C4 value "
                       "during parameter searching. [Default: 10.0]")
parser.add_option("--method", dest="minm", type='string', \
                  help="Optimization method of the C4 terms. The options "
                       "are: powell, cg, bfgs, ncg, l_bfgs_b, tnc or slsqp. "
                       "[Default: bfgs] Please check the website: "
                       "http://docs.scipy.org/doc/scipy/reference/optimize.htm"
                       " for more information if interested.")
parser.add_option("--platform", dest="platf", type='string', \
                  help="Platform used. The options are: reference, cpu, gpu "
                       "or opencl. [Default: cpu] Here we use the OpenMM "
                       "software to perform the structure minimization. "
                       "Please check OpenMM user guide for more information "
                       "if interested.")
parser.add_option("--model", dest="model", type='int', \
                  help="The metal ion complex model chosen to calculate the "
                  "RMSD value. RMSD value is the criterion for the "
                  "optimization (with a smaller RMSD value, better the "
                  "parameters). The options are: 1 or 2. "
                  "1 means a small model (only contains the metal ion and "
                  "binding heavy atoms) while 2 means a big (contains the "
                  "metal ion and heavy atoms in the ligating residues). "
                  "[Default: 1]")
(options, args) = parser.parse_args()

#lowercase the input variables
options.simupha=options.simupha.lower()
options.minm=options.minm.lower()
options.platf=options.platf.lower()

if options.minm == 'powell':
  from scipy.optimize import fmin_powell as fmin
elif options.minm == 'cg':
  from scipy.optimize import fmin_cg as fmin
elif options.minm == 'bfgs':
  from scipy.optimize import fmin_bfgs as fmin
elif options.minm == 'ncg':
  from scipy.optimize import fmin_ncg as fmin
elif options.minm == 'l_bfgs_b':
  from scipy.optimize import fmin_l_bfgs_b as fmin
elif options.minm == 'tnc':
  from scipy.optimize import fmin_tnc as fmin
elif options.minm == 'cobyla':
  from scipy.optimize import fmin_cobyla as fmin
elif options.minm == 'slsqp':
  from scipy.optimize import fmin_slsqp as fmin

#Get the metal center information from prmtop and coordinate files
#Get the metal center
prmtop, mol, atids, resids = read_amber_prm(options.pfile, options.cfile)
mask = AmberMask(prmtop, options.ion_mask)

c4terms = prmtop.parm_data['LENNARD_JONES_CCOEF']

#Get the metal ion ids
metids = []
mettyps = []
mettypinds = []
for i in mask.Selected():
  #Atom IDs
  j = i + 1
  metids.append(j)
  #Amber Atom Type
  atyp = prmtop.parm_data['AMBER_ATOM_TYPE'][i]
  mettyps.append(atyp)
  #Atom Type Index
  mettypind = prmtop.parm_data['ATOM_TYPE_INDEX'][i]
  mettypinds.append(mettypind)

mcids = []
mctyps = []
mctypinds = []

if options.model == 1: #Small model
  for i in metids:
    crdi = mol.atoms[i].crd
    atmi = mol.atoms[i].element
    radiusi = CoRadiiDict[atmi]
    for j in atids:
      if j != i:
        crdj = mol.atoms[j].crd
        atmj = mol.atoms[j].element
        dis = calc_bond(crdi, crdj)
        radiusj = CoRadiiDict[atmj]
        radiusij = radiusi + radiusj
        if (dis <= radiusij + 0.4) and (dis >= 0.1) and (atmj != 'H'):
          k = j - 1
          #Atom IDs
          mcids.append(j)
          #Amber Atom Type
          atyp = prmtop.parm_data['AMBER_ATOM_TYPE'][k]
          mctyps.append(atyp)
          #Atom Type Index
          mctypind = prmtop.parm_data['ATOM_TYPE_INDEX'][k]
          mctypinds.append(mctypind)
elif options.model == 2: #Big model
  mcresids = []
  for i in metids:
    crdi = mol.atoms[i].crd
    atmi = mol.atoms[i].element
    radiusi = CoRadiiDict[atmi]
    for j in atids:
      if j != i:
        crdj = mol.atoms[j].crd
        atmj = mol.atoms[j].element
        dis = calc_bond(crdi, crdj)
        radiusj = CoRadiiDict[atmj]
        radiusij = radiusi + radiusj
        if (dis <= radiusij + 0.4) and (dis >= 0.1) and (atmj != 'H'):
          if mol.atoms[j].resid not in mcresids:
            mcresids.append(mol.atoms[j].resid)

  for i in mcresids:
    for j in mol.residues[i].resconter:
      if mol.atoms[j].element != 'H':
        k = j - 1
        mcids.append(j)
        #Amber Atom Type
        atyp = prmtop.parm_data['AMBER_ATOM_TYPE'][k]
        mctyps.append(atyp)
        #Atom Type Index
        mctypind = prmtop.parm_data['ATOM_TYPE_INDEX'][k]
        mctypinds.append(mctypind)

#Get the atom type dictionary for people to see
mettypdict = get_typ_dict(mettypinds, mettyps)
mctypdict = get_typ_dict(mctypinds, mctyps)

#Delete the repeat elements sort the list
mettypinds = list(set(mettypinds))
mettypinds.sort()
mctypinds = list(set(mctypinds))
mctypinds.sort()

#-----------------------------------------------------------------------------#
#Get the Amber mask of the metal center complex and print it into ptraj.in file
#-----------------------------------------------------------------------------#
maskns = []
for i in mcids:
  maskn = str(mol.atoms[i].resid) + '@' + mol.atoms[i].atname
  maskns.append(maskn)
for i in metids:
  maskn = str(mol.atoms[i].resid) + '@' + mol.atoms[i].atname
  maskns.append(maskn)

maskslet = ':'
for i in range(0, len(maskns)-1):
  maskslet += maskns[i] + ','
maskslet += maskns[-1]

ptrajif = open('OptC4_ptraj.in', 'w')
print >> ptrajif, 'trajin %s' %options.cfile
print >> ptrajif, 'trajin %s' %options.rfile
print >> ptrajif, 'rms %s first out OptC4_rmsd.txt' %maskslet
print >> ptrajif, 'go'
ptrajif.close()

#Detect the C4 terms which relates to the metal center complex
ntyps = prmtop.pointers['NTYPES']
idxs = []
iddict = {}
for i in mettypinds:
  j = i - 1
  for k in mctypinds:
    l = k - 1
    if k < i:
      idx = prmtop.parm_data['NONBONDED_PARM_INDEX'][ntyps*j+l] - 1
    else:
      idx = prmtop.parm_data['NONBONDED_PARM_INDEX'][ntyps*l+j] - 1
    idxs.append(idx)
    iddict[idx] = (i, k)

idxs.sort()
initparas = [c4terms[i] for i in idxs]

#Doing optimization of the parameters, initial was the normal C4 term
xopt = fmin(get_rmsd, initparas, epsilon=options.stepsize)

print xopt

print "The final RMSD is: "
frmsd = get_rmsd(xopt)

for i in range(0, len(idxs)):
  print 'The optimal value of atomtypes ' + str(iddict[idxs[i]]) + \
        ' is ' + str(xopt[i])

print 'The following is the dictionary of the atom types: '
print mettypdict
print mctypdict

