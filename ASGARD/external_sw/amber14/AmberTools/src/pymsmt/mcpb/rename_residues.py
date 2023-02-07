"""
This module is written for detecting the disulfide bond and renaming the
residues.
"""
from pymsmtmol.cal import calc_bond
from pymsmtmol.element import resnamel

def get_diS_bond(mol, atids):

    #Residue IDs
    resids = []
    for i in atids:
      if mol.atoms[i].resid not in resids:
        resids.append(mol.atoms[i].resid)

    disul = []
    #Correct the names of the HIS, ASP, GLU, LYS, CYS
    for i in resids:
      if mol.residues[i].resname == 'CYM':
        #for atom in CYM residue
        for j in mol.residues[i].resconter:
          if mol.atoms[j].atname == 'SG':
            sgcrd = mol.atoms[j].crd
            #for every atom
            for k in atids:
              if mol.atoms[k].atname == 'SG' and k != j:
                atkcrd = mol.atoms[k].crd
                dis = calc_bond(sgcrd, atkcrd)
                if dis <= 2.50:
                  if j < k:
                    disuls.append((j, k))
                  else:
                    disuls.append((k, j))

    return disul

def rename_res(mol, atids):

    #Residue IDs
    resids = []
    for i in atids:
      if mol.atoms[i].resid not in resids:
        resids.append(mol.atoms[i].resid)

    #Correct the names of the HIS, ASP, GLU, LYS, CYS
    for i in resids:
      if mol.residues[i].resname == 'HIS':
        hasatoms = []
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          hasatoms.append(atname)
        if ('HD1' in hasatoms) and ('HE2' in hasatoms):
          mol.residues[i].resname = 'HIP'
        elif ('HD1' in hasatoms):
          mol.residues[i].resname = 'HID'
        elif ('HE2' in hasatoms):
          mol.residues[i].resname = 'HIE'
      elif mol.residues[i].resname == 'ASP':
        hasatoms = []
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          hasatoms.append(atname)
        if ('HD1' in hasatoms) or ('HD2' in hasatoms):
          mol.residues[i].resname = 'ASH'
      elif mol.residues[i].resname == 'GLU':
        hasatoms = []
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          hasatoms.append(atname)
        if ('HE1' in hasatoms) or ('HE2' in hasatoms):
          mol.residues[i].resname = 'GLH'
      elif mol.residues[i].resname == 'LYS':
        hasatoms = []
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          hasatoms.append(atname)
        if ('HZ1' not in hasatoms):
          mol.residues[i].resname = 'LYN'
      elif mol.residues[i].resname == 'CYS':
        hasatoms = []
        for j in mol.residues[i].resconter:
          atname = mol.atoms[j].atname
          hasatoms.append(atname)
        if ('HG' not in hasatoms): ##There are two different situations
          #for atom in CYS residue
          for j in mol.residues[i].resconter:
            if mol.atoms[j].atname == 'SG':
              sgcrd = mol.atoms[j].crd
              #for every atom
              for k in atids:
                if mol.atoms[k].atname == 'SG' and k != j:
                  atkcrd = mol.atoms[k].crd
                  dis = calc_bond(sgcrd, atkcrd)
                  if dis <= 2.50:
                    mol.residues[i].resname = 'CYS'
                  else:
                    mol.residues[i].resname = 'CYM'

    #rename the HN atom to H atom in amino acid residues
    for i in resids:
      resname = mol.residues[i].resname
      if resname in resnamel:
        for j in mol.residues[i].resconter:
          if mol.atoms[j].atname == 'HN':
            mol.atoms[j].atname = 'H'

    return mol

