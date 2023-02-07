"""
module for define the normally used data about elements, amino acids et al.
"""
from chemistry.periodic_table import AtomicNum
#-----------------------------------------------------------------------------

#Ion names which has the VDW parameters

#-----------------------------------------------------------------------------

movions = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Tl', 'Cu', 'Ag', 'Mo']
divions = ['Be', 'Cu', 'Ni', 'Pt', 'Zn', 'Co', 'Pd', 'Ag', 'Cr',
            'Fe', 'Mg', 'V', 'Mn', 'Hg', 'Cd', 'Yb', 'Ca', 'Sn',
            'Pb', 'Eu', 'Sr', 'Sm', 'Ba', 'Ra', 'Mo']
trivions = ['Al', 'Fe', 'Cr', 'In', 'Tl', 'Y', 'La', 'Ce', 'Pr',
          'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Er', 'Tm', 'Lu', 'As',
          'Ru', 'Mo']
tetvions = ['Zr', 'Ce', 'U', 'Pu', 'Th', 'Mo']
petvions = ['Mo']
sxavions = ['Mo']

ionnamel = movions + divions + trivions + tetvions + petvions + sxavions
ionnamel2 = [i[0] + i[1].upper() for i in ionnamel if len(i) > 1]
ionnamel = ionnamel + ionnamel2

#-----------------------------------------------------------------------------

#Residue names and their letters

#-----------------------------------------------------------------------------

resdict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B',
           'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G',
           'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
           'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'SEC': 'U',
           'TRP': 'W', 'TYR': 'Y', 'XAA': 'X', 'VAL': 'V'}

#-----------------------------------------------------------------------------

#Residue names avaiable in the AMBER FFs

#-----------------------------------------------------------------------------

resnamel = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
            'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
            'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL']

resnamecl = ['C' + i for i in resnamel]
resnamenl = ['N' + i for i in resnamel]
resnamel = resnamel + resnamecl + resnamenl

reschg = {'ALA': 0, 'ARG': 1, 'ASH': 0, 'ASN': 0, 'ASP': -1, 'CYM': -1,
          'CYS': 0, 'CYX': 0, 'GLH': 0, 'GLN': 0, 'GLU': -1, 'GLY':  0,
          'HIS': 0, 'HID': 0, 'HIE': 0, 'HIP': 1, 'ILE':  0, 'LEU':  0,
          'LYN': 0, 'LYS': 1, 'MET': 0, 'PHE': 0, 'PRO': 0,  'SER':  0,
          'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0}

reschgnt = dict([('N' + k, v+1) for k, v in reschg.iteritems()])
reschgct = dict([('C' + k, v-1) for k, v in reschg.iteritems()])

ResChgDict = reschg
ResChgDict.update(reschgnt)
ResChgDict.update(reschgct)

#-----------------------------------------------------------------------------

#Mass, making it compatible with AMBER FFs

#-----------------------------------------------------------------------------

OrganicMass = {'H': 1.008, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'S': 32.06,
               'P': 30.97, 'F': 19.00, 'Cl': 35.45, 'Br': 79.90, 'I': 126.9
              }
 
Metal1mass = {'Li': 6.94, 'Na': 22.99, 'K': 39.10, 'Rb': 85.47, 'Cs': 132.91,
              'Mo': 95.96
             }

Metal2mass = {'Be': 9.01,'Cu': 63.55,'Ni': 58.69,'Pt': 195.08,'Zn': 65.4,
              'Co': 58.93,'Pd': 106.42,'Ag': 107.87,'Cr': 52.00,'Fe': 55.85,
              'Mg': 24.305,'V':  50.94,'Mn': 54.94,'Hg': 200.59,'Cd': 112.41,
              'Yb': 173.05,'Ca': 40.08,'Sn': 118.71,'Pb': 207.2,'Eu': 151.96,
              'Sr': 87.62,'Sm': 150.36,'Ba': 137.33,'Ra': 226.03, 'Mo': 95.96
             }

Metal3mass = {'Al': 26.98,'Fe': 55.85,'Cr': 52.00,'In': 114.82,'Tl': 204.38,
              'Y': 88.91,'La': 138.91,'Ce': 140.12,'Pr': 140.91,'Nd': 144.24,
              'Sm': 150.36,'Eu': 151.96,'Gd': 157.25,'Tb': 158.93,'Dy': 162.5,
              'Er': 167.26,'Tm': 168.93,'Lu': 174.97, 'As': 74.92, 'Ru': 101.07,
              'Mo': 95.96
             }

Metal4mass = {'Hf': 178.49,'Zr': 91.22,'Ce': 140.12,'U': 238.03,'Pu': 244.06,
              'Th': 232.04, 'Mo': 95.96
             }

Mass = OrganicMass
Mass.update(Metal1mass)
Mass.update(Metal2mass)
Mass.update(Metal3mass)
Mass.update(Metal4mass)

"""
Metal1pdb = {'F-':  ('F', 'F'), 'Br-': ('BR', 'BR'), 'Cl-': ('CL', 'CL'), 
             'I-':  ('IOD', 'I'), 'Li+': ('LI', 'LI'), 'Na+': ('NA', 'NA'),
             'Rb+': ('RB', 'RB'), 'Tl+': ('TL', 'TL'), 'Cs+': ('CS', 'CS'),
             'K+':  ('K', 'K'), 'Cu+': ('CU1', 'CU'), 'Ag+': ('AG', 'AG'),
             'Au+': ('AU', 'AU')}

Metal2pdb = {'Cu2+': ('CU', 'CU'), 'Ni2+': ('NI', 'NI'), 'Pt2+': ('PT', 'PT'),
             'Zn2+': ('ZN', 'ZN'), 'Co2+': ('CO', 'CO'), 'Pd2+': ('PD', 'PD'),
             'Fe2+': ('FE2', 'FE2'), 'Mg2+': ('MG', 'MG'), 'Mn2+': ('MN', 'MN'),
             'Hg2+': ('HG', 'HG'), 'Cd2+': ('CD', 'CD'), 'Yb2+': ('YB2', 'YB2'),
             'Ca2+': ('CA', 'CA'), 'Pb2+': ('PD', 'PD'), 'Eu2+': ('EU', 'EU'),
             'Sr2+': ('SR', 'SR'), 'Ba2+': ('BA', 'BA')}

Metal3pdb = {'Al3+': ('AL', 'AL'), 'Fe3+': ('FE', 'FE'), 'Cr3+': ('CR', 'CR'),
             'In3+': ('IN', 'IN'), 'Y3+': ('Y', 'Y'), 'La3+': ('LA', 'LA'),
             'Ce3+': ('CE', 'CE'), 'Pr3+': ('PR', 'PR'), 'Sm3+': ('SM', 'SM'),
             'Eu3+': ('EU3', 'EU'), 'Gd3+': ('GD3', 'GD'), 'Tb3+': ('TB', 'TB'),
             'Lu3+': ('LU', 'LU'), 'V3+': ('V', 'V'), 'As3+': ('ARS', 'AS'),
             'Ru3+': ('RU', 'RU')}

Metal4pdb = {'Ir4+': ('IR', 'IR'), 'Mo4+': ('MO', 'MO')}
"""

#-----------------------------------------------------------------------------

#PDB names, the usually resname and atname for the metal ions

#-----------------------------------------------------------------------------



Metal1pdb = {('F', 'F'): 'F', ('BR', 'BR'): 'Br', ('CL', 'CL'): 'Cl', 
             ('IOD', 'I'): 'I', ('LI', 'LI'): 'Li', ('NA', 'NA'): 'Na',
             ('RB', 'RB'): 'Rb', ('TL', 'TL'): 'Tl', ('CS', 'CS'): 'Cs',
             ('K', 'K'): 'K', ('CU1', 'CU'): 'Cu', ('AG', 'AG'): 'Ag',
             ('AU', 'AU'): 'Au'}

Metal2pdb = {('CU', 'CU'): 'Cu', ('NI', 'NI'): 'Ni', ('PT', 'PT'): 'Pt',
             ('ZN', 'ZN'): 'Zn', ('CO', 'CO'): 'Co', ('PD', 'PD'): 'Pd',
             ('FE2', 'FE2'): 'Fe', ('MG', 'MG'): 'Mg', ('MN', 'MN'): 'Mn',
             ('HG', 'HG'): 'Hg', ('CD', 'CD'): 'Cd', ('YB2', 'YB2'): 'Yb',
             ('CA', 'CA'): 'Ca', ('PB', 'PB'): 'Pb', ('EU', 'EU'): 'Eu',
             ('SR', 'SR'): 'Sr', ('BA', 'BA'): 'Ba'}

Metal3pdb = {('AL', 'AL'): 'Al', ('FE', 'FE'): 'Fe', ('CR', 'CR'): 'Cr',
             ('IN', 'IN'): 'In', ('Y', 'Y'): 'Y', ('LA', 'LA'): 'La',
             ('CE', 'CE'): 'Ce', ('PR', 'PR'): 'Pr', ('SM', 'SM'): 'Sm',
             ('EU3', 'EU'): 'Eu', ('GD3', 'GD'): 'Gd', ('TB', 'TB'): 'Tb',
             ('LU', 'LU'): 'Lu', ('V', 'V'): 'V', ('ARS', 'AS'): 'As',
             ('RU', 'RU'): 'Ru'}

Metal4pdb = {('IR', 'IR'): 'Ir', ('MO', 'MO'): 'Mo'}

Metalpdb = Metal1pdb
Metalpdb.update(Metal2pdb)
Metalpdb.update(Metal3pdb)
Metalpdb.update(Metal4pdb)

#-----------------------------------------------------------------------------
#Covalent radii are used for determine the bonds between two atoms.
#The tolerance is 0.4 Angstrom. (<= Coradius1 + Coradius2 + 0.4)
#From Elaine C. Meng and Richard A. Lewis, Journal of Computational Chemistry,
#1991, 12(7), 891-898
#-----------------------------------------------------------------------------

CoRadiiDict = {'H': 0.23, 'He': 1.50, 'Li': 0.68, 'Be': 0.35, 'B': 0.83,
               'C': 0.68, 'N': 0.68, 'O': 0.68, 'F': 0.64, 'Ne': 1.50,
               'Na': 0.97, 'Mg': 1.10, 'Al': 1.35, 'Si': 1.20, 'P': 1.05,
               'S': 1.02, 'Cl': 0.99, 'Ar':1.51, 'K': 1.33, 'Ca': 0.99,
               'Sc': 1.44, 'Ti': 1.47, 'V': 1.33, 'Cr': 1.35, 'Mn': 1.35,
               'Fe': 1.34, 'Co': 1.33, 'Ni': 1.50, 'Cu': 1.52, 'Zn': 1.45,
               'Ga': 1.22, 'Ge': 1.17, 'As': 1.21, 'Se': 1.22, 'Br': 1.21,
               'Kr': 1.50, 'Rb': 1.47, 'Sr': 1.12, 'Y': 1.78, 'Zr': 1.56,
               'Nb': 1.48, 'Mo': 1.47, 'Tc': 1.35, 'Ru': 1.40, 'Rh': 1.45,
               'Pd': 1.50, 'Ag': 1.59, 'Cd': 1.69, 'In': 1.63, 'Sn': 1.46,
               'Sb': 1.46, 'Sb': 1.46, 'Te': 1.47, 'I': 1.40, 'Cs': 1.67,
               'Ba': 1.34, 'La': 1.87, 'Ce': 1.83, 'Pr': 1.82, 'Nd': 1.81, 
               'Pm': 1.80, 'Sm': 1.80, 'Eu': 1.99, 'Gd': 1.79, 'Tb': 1.76,
               'Dy': 1.75, 'Ho': 1.74, 'Er': 1.73, 'Tm': 1.72, 'Yb': 1.94,
               'Lu': 1.72, 'Hf': 1.57, 'Ta': 1.43, 'W': 1.37, 'Re': 1.35,
               'Os': 1.37, 'Ir': 1.32, 'Pt': 1.50, 'Au': 1.50, 'Hg': 1.70,
               'Tl': 1.55, 'Pb': 1.54, 'Bi': 1.54, 'Po': 1.68, 'Ra': 1.90,
               'Ac': 1.88, 'Th': 1.79, 'Pa': 1.61, 'U': 1.58, 'Np': 1.55,
               'Pu': 1.53, 'Am': 1.51} #91 elements


#-----------------------------------------------------------------------------

#VDW radii
#From http://periodictable.com/Properties/A/VanDerWaalsRadius.an.html 

#-----------------------------------------------------------------------------

#Neutral
VdwRadiiDict = {'H': 1.20, 'He': 1.40, 'Li': 1.82, 'C': 1.70, 'N': 1.55,
                'O': 1.52, 'F': 1.47, 'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73,
                'Si': 2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Ar': 1.88,
                'K': 2.75, 'Ni': 1.63, 'Cu': 1.40, 'Zn': 1.39, 'Ga': 1.87,
                'As': 1.85, 'Se': 1.90, 'Br': 1.85, 'Kr': 2.02, 'Pd': 1.63,
                'Ag': 1.72, 'Cd': 1.58, 'In': 1.93, 'Sn': 2.17, 'Te': 2.06,
                'I': 1.98, 'Xe': 2.16, 'Pt': 1.75, 'Au': 1.66, 'Hg': 1.55,
                'Tl': 1.96, 'Pb': 2.02, 'U': 1.86} #38 elements

#LJ parameters for the ions in bonded model
def get_ionljparadict(watermodel):

    #Monovalent ions from IOD parameter set
    monoljpara =  {
                'Li1' : (1.315, 0.00594975, 'IOD set from Li et al. JCTC, Accepted'),
                'Na1' : (1.465, 0.02909167, 'IOD set from Li et al. JCTC, Accepted'),
                'K1'  : (1.745, 0.17018074, 'IOD set from Li et al. JCTC, Accepted'),
                'Rb1' : (1.820, 0.22962229, 'IOD set from Li et al. JCTC, Accepted'),
                'Cs1' : (2.000, 0.38943250, 'IOD set from Li et al. JCTC, Accepted'),
                'Tl1' : (1.870, 0.27244486, 'IOD set from Li et al. JCTC, Accepted'),
                'Cu1' : (1.214, 0.00139196, 'IOD set from Li et al. JCTC, Accepted'),
                'Ag1' : (1.500, 0.03899838, 'IOD set from Li et al. JCTC, Accepted'),
                'H1'  : (0.925, 0.00000147, 'IOD set from Li et al. JCTC, Accepted'),
                'Mo1' : (1.20, 0.0125, 'As Zn(II) From Merz et al. JACS, 113, 8262'),
                'F-1' : (1.739, 0.16573832, 'IOD set from Li et al. JCTC, Accepted'),
                'Cl-1': (2.162, 0.53154665, 'IOD set from Li et al. JCTC, Accepted'),
                'Br-1': (2.331, 0.65952968, 'IOD set from Li et al. JCTC, Accepted'),
                'I-1' : (2.590, 0.80293907, 'IOD set from Li et al. JCTC, Accepted'),
                   }

    IonLJParaDict = monoljpara

    #Divalent ions from IOD parameter set
    diljpara = {
              'Be2': (1.168, 0.00063064, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Cu2': (1.409, 0.01721000, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Ni2': (1.373, 0.01179373, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Zn2': (1.395, 0.01491700, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Co2': (1.404, 0.01636246, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Cr2': (1.388, 0.01386171, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Fe2': (1.409, 0.01721000, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Mg2': (1.395, 0.01491700, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'V2':  (1.476, 0.03198620, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Mn2': (1.467, 0.02960343, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Hg2': (1.575, 0.06751391, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Cd2': (1.506, 0.04090549, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Ca2': (1.608, 0.08337961, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Sn2': (1.738, 0.16500296, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Sr2': (1.753, 0.17618319, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Ba2': (1.913, 0.31060194, 'IOD set from Li et al. JCTC, 2013, 9, 2733'),
              'Mo2': (1.20, 0.0125, 'As Zn(II) IOD set from Merz et al. JACS, 113, 8262')
                }

    IonLJParaDict.update(diljpara)

    #Divalent ions from CM parameter set
    if watermodel == 'tip3p':
      diljpara2 = {
                 'Pt2': (1.266, 0.00307642, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Pd2': (1.303, 0.00509941, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Ag2': (1.336, 0.00770969, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Yb2': (1.642, 0.10185975, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Pb2': (1.745, 0.17018074, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Eu2': (1.802, 0.21475916, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Sm2': (1.819, 0.22878796, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Ra2': (2.019, 0.40664608, 'CM set for TIP3P water from Li et al. JCTC, 2013, 9, 2733'), 
                   }
    elif watermodel == 'spce':
      diljpara2 = {
                 'Pt2': (1.272, 0.00334975, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Pd2': (1.305, 0.00523385, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Ag2': (1.337, 0.00780282, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Yb2': (1.634, 0.09731901, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Pb2': (1.731, 0.15989650, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Eu2': (1.786, 0.20184160, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Sm2': (1.800, 0.21312875, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Ra2': (1.980, 0.37126402, 'CM set for SPC/E water from Li et al. JCTC, 2013, 9, 2733'), 
                   }
    elif watermodel == 'tip4pew':
      diljpara2 = {
                 'Pt2': (1.251, 0.00247282, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Pd2': (1.288, 0.00417787, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Ag2': (1.323, 0.00657749, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Yb2': (1.654, 0.10888937, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Pb2': (1.758, 0.17997960, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Eu2': (1.823, 0.23213110, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Sm2': (1.838, 0.24480038, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                 'Ra2': (2.050, 0.43454345, 'CM set for TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'), 
                   }

    IonLJParaDict.update(diljpara2)

    #Trivalent and tetravalent ions from IOD parameter set
    if watermodel == 'tip3p':
      higljpara = {
                 'Al3': (1.297, 0.00471279, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Fe3': (1.386, 0.01357097, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Cr3': (1.344, 0.00848000, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'In3': (1.461, 0.02808726, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Tl3': (1.513, 0.04321029, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Y3':  (1.602, 0.08034231, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'La3': (1.718, 0.15060822, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Ce3': (1.741, 0.16721338, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Pr3': (1.733, 0.16134811, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Nd3': (1.681, 0.12564307, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Sm3': (1.659, 0.11189491, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Eu3': (1.666, 0.11617738, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Gd3': (1.623, 0.09126804, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Tb3': (1.630, 0.09509276, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Dy3': (1.609, 0.08389240, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Er3': (1.602, 0.08034231, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Tm3': (1.602, 0.08034231, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Lu3': (1.588, 0.07351892, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Hf4': (1.499, 0.03868661, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Zr4': (1.519, 0.04525501, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Ce4': (1.684, 0.12758274, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'U4':  (1.684, 0.12758274, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Pu4': (1.662, 0.11371963, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                 'Th4': (1.708, 0.14364160, 'IOD set for TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   }
    elif watermodel == 'spce':
      higljpara = {
                 'Al3': (1.296, 0.00465074, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Fe3': (1.386, 0.01357097, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Cr3': (1.343, 0.00838052, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'In3': (1.461, 0.02808726, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Tl3': (1.513, 0.04321029, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Y3':  (1.602, 0.08034231, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'La3': (1.718, 0.15060822, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Ce3': (1.741, 0.16721338, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Pr3': (1.734, 0.16207614, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Nd3': (1.681, 0.12564307, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Sm3': (1.659, 0.11189491, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Eu3': (1.666, 0.11617738, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Gd3': (1.623, 0.09126804, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Tb3': (1.630, 0.09509276, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Dy3': (1.609, 0.08389240, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Er3': (1.602, 0.08034231, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Tm3': (1.602, 0.08034231, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Lu3': (1.588, 0.07351892, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Hf4': (1.501, 0.03931188, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Zr4': (1.521, 0.04595090, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Ce4': (1.689, 0.13084945, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'U4':  (1.689, 0.13084945, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Pu4': (1.666, 0.11617738, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                 'Th4': (1.713, 0.14710519, 'IOD set for SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   }
    elif watermodel == 'tip4pew':
      higljpara = {
                 'Al3': (1.285, 0.00401101, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Fe3': (1.375, 0.01205473, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Cr3': (1.333, 0.00743559, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'In3': (1.450, 0.02545423, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Tl3': (1.502, 0.03962711, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Y3':  (1.590, 0.07447106, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'La3': (1.707, 0.14295367, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Ce3': (1.729, 0.15845086, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Pr3': (1.722, 0.15343866, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Nd3': (1.669, 0.11803919, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Sm3': (1.647, 0.10475707, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Eu3': (1.655, 0.10948690, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Gd3': (1.612, 0.08544204, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Tb3': (1.619, 0.08912336, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Dy3': (1.597, 0.07786298, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Er3': (1.590, 0.07447106, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Tm3': (1.590, 0.07447106, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Lu3': (1.577, 0.06841702, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Hf4': (1.483, 0.03393126, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Zr4': (1.503, 0.03994409, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Ce4': (1.667, 0.11679623, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'U4':  (1.667, 0.11679623, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Pu4': (1.645, 0.10359269, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                 'Th4': (1.690, 0.13150785, 'IOD set for TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   }

    IonLJParaDict.update(higljpara)

    higljpara2 = {
                'As3': (1.1778, 0.00075222, 'IOD set for TIP3P water from Li. Personal Communication.'),
                'Ru3': (1.20, 0.0125, 'As Zn(II) From Merz et al. JACS, 113, 8262'),
                'Mo3': (1.20, 0.0125, 'As Zn(II) From Merz et al. JACS, 113, 8262'),
                'Mo4': (1.20, 0.0125, 'As Zn(II) From Merz et al. JACS, 113, 8262')
                  }

    IonLJParaDict.update(higljpara2)

    return IonLJParaDict


#LJ parameter for nonbonded model
IonLJparal = ['Li1', 'Na1', 'K1',  'Rb1', 'Cs1', 'Tl1', 'Cu1', 'Ag1', 'H1', 'Mo1',
              'F-1', 'Cl-1','Br-1','I-1', 'Be2', 'Cu2', 'Ni2', 'Zn2', 'Co2','Cr2',
              'Fe2', 'Mg2', 'V2', 'Mn2', 'Hg2', 'Cd2', 'Ca2', 'Sn2', 'Sr2', 'Ba2',
              'Pt2', 'Pd2', 'Ag2', 'Yb2', 'Pb2', 'Eu2', 'Sm2', 'Ra2', 'Al3','Fe3',
              'Cr3', 'In3', 'Tl3', 'Y3', 'La3', 'Ce3', 'Pr3', 'Nd3', 'Sm3', 'Eu3',
              'Gd3', 'Tb3', 'Dy3', 'Er3', 'Tm3', 'Lu3', 'Hf4', 'Zr4', 'Ce4', 'U4',
              'Pu4', 'Th4']

IonCMparal = ['Be2', 'Cu2', 'Ni2', 'Zn2', 'Co2', 'Cr2', 'Fe2', 'Mg2', 'V2', 'Mn2',
              'Hg2', 'Cd2', 'Ca2', 'Sn2', 'Sr2', 'Ba2', 'Pt2', 'Pd2', 'Ag2', 'Yb2',
              'Pb2', 'Eu2', 'Sm2', 'Ra2']

IonHFEparal = IonLJparal
IonHFEparal.remove('H1')

IonIODparal = list(set(IonLJparal) - \
              set(['Pt2', 'Pd2', 'Ag2', 'Yb2', 'Pb2', 'Eu2', 'Sm2', 'Ra2']))

#-----------------------------------------------------------------------------

#Atom numbers

#-----------------------------------------------------------------------------
Atnum = AtomicNum

AtnumRev = dict([ (v, k) for k, v in Atnum.iteritems()])

bdld = {'CH': 1.090, 'NH': 1.010}
