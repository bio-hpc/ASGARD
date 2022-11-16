"""
This module was written for generating the resp fitting input file and doing
the RESP charge fitting.
"""
from __future__ import absolute_import
from pymsmtmol.element import resnamel, Atnum
from pymsmtmol.readpdb import get_atominfo_fpdb
from pymsmtmol.getlist import get_mc_blist
from pymsmtlib.lib import get_lib_dict
from pymsmtexp import *
import os

def read_resp_file(fname):
    chgs = []
    fp = open(fname, 'r')
    for line in fp:
      line = line.strip("\n")
      tmp = line.split(" ")
      for i in tmp:
        if (i != ''):
          chgs.append(float(i))
    fp.close()
    return chgs

def print_mol2f(resid, resname1, resname2, resconter, mol, iddict1, sddict, \
                stdict, blist_each):

    #iddict1: atom id
    #mol: atname, crd
    #sddict: atomtype
    #stdict: atom charge
    #blist_each: bond information

    mol2f = open(resname2 + '.mol2', 'w')

    print '***Generating the ' + resname2 + '.mol2 file...'

    ##1. molecule information
    print >> mol2f, "@<TRIPOS>MOLECULE"
    print >> mol2f, resname2
    print >> mol2f, '%5d%6d%6d%6d%6d' %(len(resconter), len(blist_each),
                                       1, 0, 0) #atom number and bond number
    print >> mol2f, 'SMALL'
    print >> mol2f, 'RESP Charge'
    print >> mol2f, ' '
    print >> mol2f, ' '

    ##2. atom information
    print >> mol2f, '@<TRIPOS>ATOM'
    for atm in resconter:
      #new atom id
      atid = mol.atoms[atm].atid
      natid = iddict1[atid]
      atname = mol.atoms[atm].atname
      crd = mol.atoms[atm].crd
      key = str(resid) + '-' +  resname1 + '-' + atname
      atomtype = sddict[key]
      chg = stdict[key]
      print >> mol2f, '%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f'\
                      %(natid, atname, crd[0], crd[1], crd[2], atomtype,
                       1, resname2, chg)

    ##3. bond information
    print >> mol2f, '@<TRIPOS>BOND'
    for bonds in range(0, len(blist_each)):
      print >> mol2f, '%6d%5d%5d%2d' %(bonds+1, blist_each[bonds][0],
                      blist_each[bonds][1], 1) #all as single bond

    ##4. substructure information
    print >> mol2f, '@<TRIPOS>SUBSTRUCTURE'
    print >> mol2f, '     1', resname2, '        1 TEMP' + \
                    '              0 ****  ****    0 ROOT'
    mol2f.close()

def gene_resp_input_file(lgpdbf, ionids, stfpf, ffchoice, mol2fs,
                         chgmod, ionchgr):

    libdict, chargedict = get_lib_dict(ffchoice)

    for mol2f in mol2fs:
      libdict1, chargedict1 = get_lib_dict(mol2f)
      libdict.update(libdict1)
      chargedict.update(chargedict1)

    #-------------------------------------------------------------------------
    ##############RESP1.IN file###############################################
    #-------------------------------------------------------------------------

    print "***Generating the 1st stage resp charge fitting input file..."

    #print the 1st part, the title
    fresp1 = open('resp1.in', 'w')
    print >> fresp1, "Resp charges for organic molecule"
    print >> fresp1, " "
    print >> fresp1, " &cntrl"
    print >> fresp1, " "
    print >> fresp1, " nmol = 1,"
    print >> fresp1, " ihfree = 1,"
    print >> fresp1, " ioutopt = 1,"
    print >> fresp1, " "
    print >> fresp1, " &end"
    print >> fresp1, "    1.0"
    print >> fresp1, "Resp charges for organic molecule"

    mol, atids, resids = get_atominfo_fpdb(lgpdbf)

    blist = get_mc_blist(mol, atids, ionids, stfpf)

    #Get the total charge of the system
    totchg = 0.0

    for i in resids:
      resname = mol.residues[i].resname
      reschg = chargedict[resname]
      totchg = totchg + reschg

    print >> fresp1, "%5d" %int(totchg),
    print >> fresp1, "%4d" %len(atids)

    #print the 2nd part, the free and fozen atoms

    #new atids
    natids = [i for i in range(1, len(atids)+1)] 

    iddict = {}

    for i in range(0, len(atids)):
      iddict[atids[i]] = natids[i]

    for i in atids:
      element = mol.atoms[i].element
      elenum = Atnum[element]
      iddict[i] = (iddict[i], elenum)

    #for ACE, NME and GLY residues
    acedict = {}
    nmedict = {}
    glydict = {}
    ace1st = 0
    nme1st = 0
    gly1st = 0

    for i in resids:
      if mol.residues[i].resname == 'ACE':
        if ace1st < 1:
        #get the ids in the first ACE group
          for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HH31', 'CH3', 'C', 'O']):
              iddict[j] = (iddict[j][0], iddict[j][1], 0)
              acedict[atname] = iddict[j][0]
            elif (atname in ['HH32', 'HH33']):
              iddict[j] = (iddict[j][0], iddict[j][1], acedict['HH31'])
        else:
          for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HH31', 'CH3', 'C', 'O']):
              iddict[j] = (iddict[j][0], iddict[j][1], acedict[atname])
            elif (atname in ['HH32', 'HH33']):
              iddict[j] = (iddict[j][0], iddict[j][1], acedict['HH31'])
        ace1st = ace1st + 1
      elif mol.residues[i].resname == 'NME':
        if nme1st < 1:
        #get the ids in the first ACE group
          for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'CH3', 'HH31', 'H']):
              iddict[j] = (iddict[j][0], iddict[j][1], 0)
              nmedict[atname] = iddict[j][0]
            elif (atname in ['HH32', 'HH33']):
              iddict[j] = (iddict[j][0], iddict[j][1], nmedict['HH31'])
        else:
          for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'CH3', 'HH31', 'H']):
              iddict[j] = (iddict[j][0], iddict[j][1], nmedict[atname])
            elif (atname in ['HH32', 'HH33']):
              iddict[j] = (iddict[j][0], iddict[j][1], nmedict['HH31'])
        nme1st = nme1st + 1
      #for GLY residues
      elif mol.residues[i].resname == 'GLY':
        if gly1st < 1:
        #get the ids of the first GLY group
          for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'H', 'C', 'O', 'CA', 'HA2']):
              iddict[j] = (iddict[j][0], iddict[j][1], 0)
              glydict[atname] = iddict[j][0]
            elif (atname in ['HA3']):
              iddict[j] = (iddict[j][0], iddict[j][1], glydict['HA2'])
        else:
          for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'H', 'C', 'O', 'CA', 'HA2']):
              iddict[j] = (iddict[j][0], iddict[j][1], glydict[atname])
            elif (atname in ['HA3']):
              iddict[j] = (iddict[j][0], iddict[j][1], glydict['HA2'])
        gly1st = gly1st + 1
      #for other residues to get the hydrogen equal information 
      else:
        #for each residue to get the bonds inside the residue
        blistinres = []
        atinres = mol.residues[i].resconter
        for bondinfo in blist:
          at1 = bondinfo[0]
          at2 = bondinfo[1]
          bondinfo = [at1, at2]
          if set(atinres).intersection(set(bondinfo)) == set(bondinfo):
            blistinres.append(bondinfo)
 
        #for the bonds inside a residue, to see whether there are two atoms
        #connect to the same atoms
        for j in range(0, len(blistinres)):
          blistj = blistinres[j]
          for k in range(j+1, len(blistinres)):
            blistk = blistinres[k]
 
            if list(set(blistj).intersection(set(blistk))) != []:
            #if there are atoms connect to the same atom, get the two atoms
              unionjk = set(blistj) | set(blistk)
              intersecjk = set(blistj) & set(blistk)
              diffset = sorted(list(unionjk - intersecjk))
 
              at1 = diffset[0]
              at2 = diffset[1]
              atc = list(intersecjk)[0]
 
              if (iddict[at1][1] == 1) and (iddict[at2][1] == 1) and (iddict[atc][1] == 6): 
              #if it is a CH2 group or CH3 group
                if (len(iddict[at2]) == 2):
                #if there is no restriction index exist for atom 2
                  iddict[at2] = (iddict[at2][0], iddict[at2][1], iddict[at1][0])
                  #if there is no restriction index exist for atom 1
                  if len(iddict[at1]) == 2:
                    iddict[at1] = (iddict[at1][0], iddict[at2][1], 0)
                  if len(iddict[atc]) == 2:
                    iddict[atc] = (iddict[atc][0], iddict[atc][1], 0)

    #other atoms are frozen
    for i in atids:
      if (len(iddict[i]) == 2):
        iddict[i] = (iddict[i][0], iddict[i][1], -99)

    for i in atids:
      if iddict[i][2] == -99:
        print >> fresp1, "%5d" %iddict[i][1],
        print >> fresp1, "%4s" %'0'
      else:
        print >> fresp1, "%5d" %iddict[i][1],
        print >> fresp1, "%4s" %iddict[i][2]

    #print the 3rd part, the atoms related to ACE and NME
    for i in resids:
      resname = mol.residues[i].resname
      #print resname
      if resname in ['ACE', 'NME', 'GLY']:
        print >> fresp1, "%5d" %len(mol.residues[i].resconter),
        print >> fresp1, "%9s" %"0.00000"
        print >> fresp1, "",
        for j in mol.residues[i].resconter:
          print >> fresp1, "%4d%5d" %(1, iddict[j][0]),
        print >> fresp1, "\n",

    #add the 4th part, backbone restriction
    if chgmod == 0:
      pass
    elif (chgmod in [1, 2, 3]):
      for i in resids:
        resname = mol.residues[i].resname
        if (resname in resnamel) and (resname not in ['ACE', 'GLY', 'NME']):
          for j in range(0, len(mol.residues[i].resconter)):
            #get new atom id
            atj = mol.residues[i].resconter[j]
            natj = iddict[atj][0]
            #get its charge
            atnamej = resname + '-' + mol.atoms[atj].atname
            chg = libdict[atnamej][1]
            #atom name
            atnamejs = mol.atoms[atj].atname
            #print the charge of restrict atom
            if chgmod == 1:
              if (atnamejs in ['CA', 'N', 'C', 'O']):
                print >> fresp1, "%5d%10.5f" %(1, chg)
                print >> fresp1, "%5d%5d" %(1, natj)
            elif chgmod == 2:
              if (atnamejs in ['CA', 'H', 'HA', 'N', 'C', 'O', 'HN']):
                print >> fresp1, "%5d%10.5f" %(1, chg)
                print >> fresp1, "%5d%5d" %(1, natj)
            elif chgmod == 3:
              if (atnamejs in ['CA', 'H', 'HA', 'N', 'C', 'O', 'HN', 'CB']):
                print >> fresp1, "%5d%10.5f" %(1, chg)
                print >> fresp1, "%5d%5d" %(1, natj)
    else:
      raise pymsmtError('Please choose chgmod among 0, 1, 2 and 3.')

    #add the 5th part, the restriction of the metal ion
    if (ionchgr == 1):
      for i in ionids:
        resname = mol.atoms[i].resname
        ionchg = chargedict[resname]
        print >> fresp1, "%5d%10.5f" %(int(1), ionchg)
        print >> fresp1, "%5d%5d" %(int(1), iddict[i][0])

    print >> fresp1, "\n"
    print >> fresp1, "\n"
    fresp1.close()

    #-------------------------------------------------------------------------
    ####################RESP2.IN file#########################################
    #-------------------------------------------------------------------------
    print "***Generating the 2nd stage resp charge fitting input file..."

    #print the 1st part, the title
    fresp2 = open('resp2.in', 'w')
    print >> fresp2, "Resp charges for organic molecule"
    print >> fresp2, " "
    print >> fresp2, " &cntrl"
    print >> fresp2, " "
    print >> fresp2, " nmol = 1,"
    print >> fresp2, " ihfree = 1,"
    print >> fresp2, " ioutopt = 1,"
    print >> fresp2, " iqopt = 2,"
    print >> fresp2, " qwt = 0.001,"
    print >> fresp2, " "
    print >> fresp2, " &end"
    print >> fresp2, "    1.0"
    print >> fresp2, "Resp charges for organic molecule"
    print >> fresp2, "%5d" %int(totchg),
    print >> fresp2, "%4d" %len(atids)

    ##forzen the group which is CH2 or CH3 group
    for i in resids:
      if mol.residues[i].resname == 'ACE':
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          if (atname in ['C', 'O']):
            iddict[j] = (iddict[j][0], iddict[j][1], -99)
      elif mol.residues[i].resname == 'NME':
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          if (atname in ['N', 'H']):
            iddict[j] = (iddict[j][0], iddict[j][1], -99)
      elif mol.residues[i].resname == 'GLY':
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          if (atname in ['N', 'H', 'C', 'O']):
            iddict[j] = (iddict[j][0], iddict[j][1], -99)

    #print the 2nd part, the free or frozen information
    for i in atids:
      print >> fresp2, "%5d" %iddict[i][1],
      print >> fresp2, "%4s" %iddict[i][2]

    #print the 3rd part, the group restriction
    for i in resids:
      resname = mol.residues[i].resname
      #print resname
      if resname in ['ACE', 'NME', 'GLY']:
        print >> fresp2, "%5d" %len(mol.residues[i].resconter),
        print >> fresp2, "%9s" %"0.00000"
        print >> fresp2, "",
        for j in mol.residues[i].resconter:
          print >> fresp2, "%4d%5d" %(1, iddict[j][0]),
        print >> fresp2, "\n",

    #add the 4th part, the backbone restriction
    if chgmod == 0:
      pass
    elif (chgmod in [1, 2, 3]):
      for i in resids:
        resname = mol.residues[i].resname
        if (resname in resnamel) and (resname != 'GLY'):
        #if resname not in ['ACE', 'NME']:
          for j in range(0, len(mol.residues[i].resconter)):
            #get new atom id
            atj = mol.residues[i].resconter[j]
            natj = iddict[atj][0]
            #get charge of the atom
            atnamej = resname + '-' + mol.atoms[atj].atname
            chg = libdict[atnamej][1]
            #get the atom name
            atnamejs = mol.atoms[atj].atname
            if chgmod == 1:
              if atnamejs in ['CA', 'N', 'C', 'O']:
                print >> fresp2, "%5d%10.5f" %(1, chg)
                print >> fresp2, "%5d%5d" %(1, natj)
            elif chgmod == 2:
              if atnamejs in ['CA', 'H', 'HA', 'N', 'C', 'O', 'HN']:
                print >> fresp2, "%5d%10.5f" %(1, chg)
                print >> fresp2, "%5d%5d" %(1, natj)
            elif chgmod == 3:
              if atnamejs in ['CA', 'H', 'HA', 'N', 'C', 'O', 'HN', 'CB']:
                print >> fresp2, "%5d%10.5f" %(1, chg)
                print >> fresp2, "%5d%5d" %(1, natj)
    else:
      raise pymsmtError('Please choose chgmod among 0, 1, 2 and 3.')

    #add the 5th part, the restriction of metal ion
    if (ionchgr == 1):
      for i in ionids:
        resname = mol.atoms[i].resname
        ionchg = chargedict[resname]
        print >> fresp2, "%5d%10.5f" %(int(1), ionchg)
        print >> fresp2, "%5d%5d" %(int(1), iddict[i][0])

    print >> fresp2, "\n"
    print >> fresp2, "\n"
    fresp2.close()

def resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids,\
                 ffchoice, mol2fs, metcenres2, chgmod, ionchgr):

    print "******************************************************************"
    print "*                                                                *"
    print "*======================RESP Charge fitting=======================*"
    print "*                                                                *"
    print "******************************************************************"

    gene_resp_input_file(lgpdbf, ionids, stfpf, ffchoice, mol2fs,
                         chgmod, ionchgr)

    #-------------------------------------------------------------------------
    ####################RESP charge fitting###################################
    #-------------------------------------------------------------------------

    print '***Doing the RESP charge fiting...'

    espf = mklogf.strip('.log') + '.esp'
    os.system("espgen -i %s -o %s" %(mklogf, espf))
    os.system("resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg \
               -e %s -s resp1_calc.esp" %espf)
    os.system("resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg \
              -t resp2.chg -e %s -s resp2_calc.esp" %espf)

    #-------------------------------------------------------------------------
    ####################Collecting the atom type and charge data##############
    #-------------------------------------------------------------------------

    #------------Atom type----------
    sddict = {} #get the atom type information from standard model
    scf = open(stfpf, 'r')
    for line in scf:
      if line[0:4] != "LINK":
        line = line.strip('\n')
        line = line.split(' ')
        line = [i for i in line if i != '']
        if len(line[-1]) == 1:
          line[-1] = line[-1] + ' '
        sddict[line[0]] = line[-1]
    scf.close()

    #------------Charge-------------
    chgs = read_resp_file('resp2.chg')

    metcenres1 = [] #original name of the metal center residue
    stlist = [] #get the atom name list from the standard model
    stf = open(stfpf, 'r')
    for line in stf:
      if line[0:4] != "LINK":
        line = line.strip('\n')
        line = line.split()
        stlist.append(line[0])
        line = line[0].split('-')
        lresname = line[0] + '-' + line[1]
        if lresname not in metcenres1:
          metcenres1.append(lresname)
    stf.close()

    llist = [] #get the atom name list from the large model
    lf = open(lgfpf, 'r')
    for line in lf:
      line = line.strip('\n')
      llist.append(line)
    lf.close()

    ldict = {} #get charge of the large model, one-to-one relationship
    for i in range(0, len(llist)):
      ldict[llist[i]] = chgs[i]

    stdict = {} #get the charge of the standard model
    for i in ldict.keys():
      if i in stlist:
        stdict[i] = ldict[i]

    #-------------------------------------------------------------------------
    ####################Checking Models#########################
    #-------------------------------------------------------------------------

    print "=========================Checking models=========================="
    print '***Check the large model...'
    if len(chgs) != len(llist):
      raise pymsmtError('Error: the charges and atom numbers are mismatch '
                        'for the large model!')
    else:
      print 'Good. The charges and atom numbers are match for the ' + \
            'large model.'
      print 'Good. There are ' + str(len(llist)) + ' atoms in the ' + \
            'large model.'

    print '***Check the standard model...'
    if len(stlist) != len(stdict):
      raise pymsmtError('Error: the charges and atom numbers are mismatch '
                        'for the standard model!')
    else:
      print 'Good. The charges and atom numbers are match for the ' + \
            'standard model.'
      print 'Good. There are ' + str(len(stlist)) + ' atoms in the ' + \
            'standard model.'

    print '***Check the residue names provided...'
    if len(metcenres1) != len(metcenres2):
      print 'You gave the residue names: ', str(metcenres2)
      print 'Database had them: ', str(metcenres1)
      raise pymsmtError('Error: The number of the residue names given is '
                        'mismatch the database!')
    else:
      print 'Good. The number of the residue names given matches ' + \
            'the database.'

    #-------------------------------------------------------------------------
    ####################Building mol2 files for modeling######################
    #-------------------------------------------------------------------------
 
    #Load the force field
    libdict, chargedict = get_lib_dict(ffchoice)

    ##get the bondlist
    mol, atids, resids = get_atominfo_fpdb(stpdbf) #from standard pdb

    blist = get_mc_blist(mol, atids, ionids, stfpf)

    blist2 = [(i[0], i[1]) for i in blist]

    print "=======================Building mol2 files========================"
    #for each residue, print out the mol2 file
    for i in range(0, len(resids)):

      resconter = mol.residues[resids[i]].resconter #atom ids in the residue
      resname1 =  mol.residues[resids[i]].resname #load residue name
      resname2 = metcenres2[i] #new residue name

      #Backbone atoms are use the AMBER backbone atom type
      if resname1 in resnamel:
        for bbatm in ['N', 'H', 'CA', 'HA','C', 'O']:
          key1 = str(resids[i]) + '-' + resname1 + '-' + bbatm
          key2 = resname1 + '-' + bbatm
          sddict[key1] = libdict[key2][0]

      #New id dict
      iddict1 = {}
      for j in range(0, len(resconter)):
        iddict1[resconter[j]] = j + 1

      #Bond list for each residue with new atom ids
      blist_each = []
      for k in blist2:
        if set(resconter) & set(k) == set(k):
          blist_each.append((iddict1[k[0]], iddict1[k[1]]))

      print_mol2f(resids[i], resname1, resname2, resconter, mol, iddict1, \
                  sddict, stdict, blist_each)


