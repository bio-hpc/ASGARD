"""
This module contains the function to model leap input files for use.
"""
from __future__ import absolute_import
from pymsmtmol.readpdb import get_atominfo_fpdb, writepdb
from pymsmtmol.readmol2 import get_atominfo
from pymsmtmol.getlist import get_mc_blist
from pymsmtmol.element import resnamel, IonHFEparal, IonCMparal, IonIODparal
from mcpb.rename_residues import rename_res, get_diS_bond
from pymsmtexp import *
import warnings
import os

##############################################################################
# Related functions
##############################################################################
def gene_ion_libfile(resname, atname, element, charge):

    atomtyp = element + str(charge) + '+'
    ionfname = '%s.cmd' %resname
    ionf = open(ionfname, 'w')
    print >> ionf, 'i = createAtom   %s  %s  %s' %(atname, atomtyp,
                                                   str(charge))
    print >> ionf, 'set i    element %s' %element
    print >> ionf, 'set i    position { 0 0 0 }'
    print >> ionf, 'r = createResidue %s' %resname
    print >> ionf, 'add r i'
    print >> ionf, '%s = createUnit %s' %(resname, resname)
    print >> ionf, 'add %s r' %resname
    print >> ionf, 'saveOff %s ./%s' %(resname, resname + '.lib')
    print >> ionf, 'quit'
    ionf.close()

    os.system('tleap -s -f %s.cmd > %s.log' %(resname, resname))

def get_frcmod_fname(element, charge, watermodel, paraset):
    ##Get the frcmod file which will be loaded
    #Check whether there is parameter for the ion

    atomtyp0 = element + str(charge)

    if paraset in ['iod', '12_6_4']:
      if atomtyp0 not in IonIODparal:
        raise pymsmtError('There is no %s parameter set for %s ion with '
                          '%+d charge for %s water model in the '
                          'database.' \
                          %(paraset.upper(), element, charge, watermodel))
    if paraset in ['hfe']:
      if atomtyp0 not in IonHFEparal:
        raise pymsmtError('There is no %s parameter set for %s ion with '
                          '%+d charge for %s water model in the '
                          'database.' \
                          %(paraset.upper(), element, charge, watermodel))
    if paraset in ['cm']:
      if charge in [1, -1]:
        warnings.warn('There is no CM parameter set for %+d ions. '
                      'Here use HFE parameter set instead, as '
                      'author suggested.' %charge, pymsmtWarning)
        paraset = 'hfe' #For +1 ion, hfe is the default
      elif charge in [3, 4]:
        warnings.warn('There is no CM parameter set for %+d ions. '
                      'Here use IOD parameter set instead, as author '
                      'suggested.' %charge, pymsmtWarning)
        paraset = 'iod' #For +1 ion, iod is the default
      elif charge == 2:
        if atomtyp0 not in IonCMparal:
          raise pymsmtError('There is no %s parameter set for %s ion with '
                            '%+d charge for %s water model in the '
                            'database.' \
                            %(paraset, element, charge, watermodel))

    #For +1 metal ions
    if charge in [1, -1]:
      frcmodf = 'frcmod.ions1lsm_'
      if paraset == 'hfe':
          frcmodf = frcmodf + paraset + '_' + watermodel
      elif paraset == 'iod':
        frcmodf = frcmodf + paraset
      elif paraset == '12_6_4':
          frcmodf = frcmodf + '1264_' + watermodel
    #For +2 metal ions, cm is the default
    elif charge == 2:
      frcmodf = 'frcmod.ionslrcm_'
      if paraset in ['hfe', 'cm']:
          frcmodf = frcmodf + paraset + '_' + watermodel
      elif paraset == 'iod':
        frcmodf = frcmodf + paraset
      elif paraset == '12_6_4':
        frcmodf = 'frcmod.ionslm_'
        frcmodf = frcmodf + '1264_' + watermodel
    #For +3 and +4 metal ions
    elif charge in [3, 4]:
      frcmodf = 'frcmod.ions34lsm_'
      if paraset == 'hfe':
          frcmodf = frcmodf + paraset + '_' + watermodel
      elif paraset == 'iod':
        frcmodf = frcmodf + 'iod'
      elif paraset == '12_6_4':
        frcmodf = frcmodf + '1264_' + watermodel

    return frcmodf

##############################################################################
#The following function has ability to generate three different leaprc files
#with model variable to switch them:
#0) For the bonded model with refitting the charge
#1) For the normal nonbonded model with refitting the charge
#2) For the normal nonbonded model without refitting the charge
##############################################################################

def gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids,\
                ionmol2fs, ioninf, mcresname, naamol2fs, ff_choice, frcmodfs,
                finfcdf, ileapf, model, watermodel='tip3p', paraset='cm'):

    print "******************************************************************"
    print "*                                                                *"
    print "*=================Generating input file for leap=================*"
    print "*                                                                *"
    print "******************************************************************"

    #---------------------Generate the new pdb file--------------------------

    mol0, atids0, resids0 = get_atominfo_fpdb(orpdbf)
    mol, atids, resids = get_atominfo_fpdb(orpdbf)

    #rename the residue names into AMBER style, e.g. HIS --> HID, HIP, HIE
    mol = rename_res(mol, atids)

    #get the disulfur bond information
    disul = get_diS_bond(mol, atids)

    #resname the old residue names to new ones if it is not metal ion
    if model in [1, 2]:
      metcenres1 = [] #Old residue ids
      fp = open(stfpf, 'r')
      for line in fp:
        if line[0:4] != "LINK":
          line = line.split('-')
          if int(line[0]) not in metcenres1:
            metcenres1.append(int(line[0]))
      fp.close()
      metcenres2 = mcresname #New residue names
      resndict = {}
      resns = []
      for i in range(0, len(metcenres1)):
          resndict[metcenres1[i]] = metcenres2[i]
          resns.append(metcenres2[i])

      for i in resndict.keys():
        mol.residues[i].resname = resndict[i]

    writepdb(mol0, mol, atids, fipdbf)

    #----------------------------get the atom names which changed atom type
    if model in [1, 2]:
      atomdefs = []
      fp0 = open(stfpf, 'r')
      for line in fp0:
        if line[0:4] != "LINK":
          line = line.strip('\n')
          line = line.split()
          if line[2] != line[4]:
            atnewtype = line[4]
            element = mol.atoms[int(line[1])].element
            atomdefs.append((atnewtype, element))
      fp0.close()
    #---------------------------Get the bond information
    if model == 1:
      mol2, atids2, resids2 = get_atominfo_fpdb(stpdbf)
      blist = get_mc_blist(mol2, atids2, ionids, stfpf)
      blist1 = [(i[0], i[1]) for i in blist]

    #----------------------Generate the lib file and get the frcmod file name
    if model == 2:
      frcmodfs = []
      for ionmol2f in ionmol2fs:
        ionmol, ionatids, ionresids = get_atominfo(ionmol2f)
        for i in ionatids:
          element = ionmol.atoms[i].element
          chg = int(ionmol.atoms[i].charge)
          frcmodf = get_frcmod_fname(element, chg, watermodel, paraset)
          if frcmodf not in frcmodfs:
            frcmodfs.append(frcmodf)
    elif model == 3 and ioninf != []:
      #get the metal information
      metresns = ioninf[0::4]
      metatns = ioninf[1::4]
      metelmts = ioninf[2::4]
      metelmts = [i[0] + i[1:].lower() for i in metelmts]
      metchgs = ioninf[3::4]
      metchgs = [int(i) for i in metchgs]
      #check the charge of the metal ions
      for metchg in metchgs:
        if metchg < -1 or metchg > 4:
          raise pymsmtError('Could not deal with atomic ion which has charge '
                            'less than -1 or bigger than +4.')
      frcmodfs = []
      for i in range(0, len(metresns)):
        if metchgs[i] > 1: #if it is -1 or +1 ions, no need to create the lib file
          gene_ion_libfile(metresns[i], metatns[i], metelmts[i], metchgs[i])
          frcmodf = get_frcmod_fname(metelmts[i], metchgs[i], watermodel, paraset)
          if frcmodf not in frcmodfs:
            frcmodfs.append(frcmodf)

    #-----------------------Generate the leap input file-----------------------

    print 'Generating the leap input file...'

    ##Generate the tleap.in file
    lp = open(ileapf, 'w')
    if ff_choice in ['ff94', 'ff99', 'ff99SB', 'ff03', 'ff10']:
      print >> lp, 'source oldff/leaprc.%s' %ff_choice
    elif ff_choice in ['ff03.r1', 'ff12SB', 'ff14SB']:
      print >> lp, 'source leaprc.%s' %ff_choice

    print >> lp, 'source leaprc.gaff'
    #Add atom types, for bonded model
    if model in [1, 2]:
      if atomdefs != []:
        print >> lp, 'addAtomTypes {'
        for i in atomdefs:
          print >> lp, '        { "%s"  "%s" "sp3" }' %(i[0], i[1])
        print >> lp, '}'

    #load lib and frcmod files for monovalent ions (for salt)
    print >> lp, 'loadoff atomic_ions.lib'
    print >> lp, 'loadamberparams frcmod.ions1lsm_hfe_%s' %watermodel

    #Load mol2 file for the refitting charge residues
    if model in [1, 2]:
      for i in resns:
        print >> lp, '%s = loadmol2 %s.mol2' %(i, i)
    elif model == 3:
      for i in naamol2fs:
        print >> lp, '%s = loadmol2 %s.mol2' %(i, i)

    #Load frcmod files for non-standard residues and metal site
    for i in frcmodfs:
      print >> lp, 'loadamberparams %s' %i

    if model == 1:
      print >> lp, 'loadamberparams %s' %finfcdf
    elif model == 2:
      for frcmodf in frcmodfs:
        print >> lp, 'loadamberparams %s' %frcmodf
    elif model == 3:
      for metresn in metresns:
        print >> lp, 'loadoff %s.lib' %metresn
      for frcmodf in frcmodfs:
        print >> lp, 'loadamberparams %s' %frcmodf

    #load pdb file
    print >> lp, 'mol = loadpdb %s' %fipdbf

    ##The Disulfur Bond information
    if disul != []:
      for i in disul:
        at1 = i[0]
        at2 = i[1]
        resid1 = mol.atoms[at1].resid
        resid2 = mol.atoms[at2].resid
        atname1 = mol.atoms[at1].atname
        atname2 = mol.atoms[at2].atname
        print >> lp, 'bond', 'mol.' + str(resid1) + '.' + atname1, 'mol.' + \
                 str(resid2) + '.' + atname2

    ##Bond including the metal ions
    if model == 1:
      for bond in blist1:
        if list(set(bond) & set(ionids)) != []:
          at1 = bond[0]
          at2 = bond[1]
          resid1 = mol.atoms[at1].resid
          resid2 = mol.atoms[at2].resid
          atname1 = mol.atoms[at1].atname
          atname2 = mol.atoms[at2].atname
          print >> lp, 'bond', 'mol.' + str(resid1) + '.' + atname1, 'mol.' + \
                   str(resid2) + '.' + atname2

    ##Nonstandard residues with nearby residues
    if model in [1, 2]:
      for i in metcenres1:
        resname = mol0.residues[i].resname
        print 'Renamed residues includes: ' + str(i) + '-' + resname
        if resname in resnamel:
          print >> lp, 'bond', 'mol.' + str(i-1) + '.C', 'mol.' + str(i) + '.N'
          print >> lp, 'bond', 'mol.' + str(i) + '.C', 'mol.' + str(i+1) + '.N'

    #Save dry structure
    print >> lp, 'savepdb mol %s_dry.pdb' %gname
    print >> lp, 'saveamberparm mol %s_dry.prmtop %s_dry.inpcrd' \
                 %(gname, gname)
    #Solvatebox
    if watermodel == 'tip3p':
      print >> lp, 'solvatebox mol TIP3PBOX 10.0'
    elif watermodel == 'spce':
      print >> lp, 'solvatebox mol SPCBOX 10.0'
      print >> lp, 'loadamberparams frcmod.spce'
    elif watermodel == 'tip4pew':
      print >> lp, 'solvatebox mol TIP4PEWBOX 10.0'
      print >> lp, 'loadamberparams frcmod.tip4pew'
    #Add counter ions
    print >> lp, 'addions mol Na+ 0'
    print >> lp, 'addions mol Cl- 0'
    #Save solvated structure
    print >> lp, 'savepdb mol %s_solv.pdb' %gname
    print >> lp, 'saveamberparm mol %s_solv.prmtop %s_solv.inpcrd' \
                 %(gname, gname)
    print >> lp, 'quit'
    print >> lp, ' '

    lp.close()

    print 'Finish generating the leap input file.'

    #tleap
    #os.system('tleap -s -f %s > %s' %(ileapf, oleapf))

