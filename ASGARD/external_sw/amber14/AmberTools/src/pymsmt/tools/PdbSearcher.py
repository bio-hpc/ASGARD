#!/usr/bin/env python
# Filename: PdbSearcher.py
"""
This is the PdbSeacher.py program written by Pengfei Li in Merz Research Group
in Michigan State University. All Rights Reserved. It is a re-written python
version of Pdbseacher in MTK++. It is designed to find the metal center in the
PDB files and collecting the information and generate the metal center complex
pdb files for each center (with metal ion and ligating residues). The program
are not gurantee to work due to bug may exist. Suggestions and bug reports are
welcome to send to Pengfei Li (Email address: ldsoar1990@gmail.com).

Please cite the following paper if you use the software to perform the
modeling:

The orignal Pdbseacher software is come from:
** M. B. Peters, Y. Yang, B. Wang, L. Fusti-Molnar, M. N. Weaver, K. M. Merz,
   JCTC, 2010, 6, 2935-2947
"""
from pymsmtmol.readpdb import get_atominfo_fpdb, writepdbatm
from pymsmtmol.element import Metalpdb, CoRadiiDict, resdict
from pymsmtmol.mol import pdbatm
from pymsmtmol.cal import calc_bond, det_geo
from optparse import OptionParser
from mcpb.title import print_title
import os

#==============================================================================
# Print the title
#==============================================================================

print_title('PdbSearcher.py')

#==============================================================================
# Setting the options
#==============================================================================

parser = OptionParser("usage: -i/--ion ionname -l/--list input_file \n"
                      "       -e/--env environment_file \n"
                      "       -s/--sum summary_file \n"
                      "       [-c/--cut cutoff]")
parser.add_option("-i", "--ion", type='string', dest="ionname",
                  help="Element symbol of ion, e.g. Zn")
parser.add_option("-l", "--list", type='string', dest="inputf",
                  help="List file name, list file contains one PDB file name "
                       "per line")
parser.add_option("-e", "--env", type='string', dest='envrmtf',
                  help="Environment file name. An environment file is used to "
                       "store the metal center environment information such "
                       "as ligating atoms, distance, geometry etc. For each "
                       "bond, there is a record.")
parser.add_option("-s", "--sum", type='string', dest='sumf',
                  help="Summary file name. A summary file is used to store "
                       "the metal center summary information such as metal "
                       "center geometry, ligating residues etc. For each "
                       "metal center there is a record.")
parser.add_option("-c", "--cut", type='float', dest='cutoff',
                  help="Optional. The cut off value used to detect the bond "
                       "between metal ion and ligating atoms. The unit is "
                       "Angstroms. If there is no value specified, the "
                       "default algorithm will be used. The default algorithm "
                       "recognizes the bond when its distance is no less than "
                       "0.1 (smaller than 0.1 usually indicates a low quality "
                       "structure) and no bigger than the covalent radius sum "
                       "of the two atoms with a tolerance of 0.4.")
(options, args) = parser.parse_args()

#==============================================================================
# Read in th pdb file names from input file
#==============================================================================
#pdb file name list
pdbfnl = []

#read the input file list
fp = open(options.inputf, 'r')
for line in fp:
    line = line.strip('\n').strip()
    pdbfnl.append(line)
fp.close()

#==============================================================================
# Print the title of each file
#==============================================================================

ionname = options.ionname

#Transfer the metal ion name
if len(ionname) == 2:
    ionname = ionname[0] + ionname[1:].lower()

print "The ionname you chosen is : " + ionname

if options.cutoff != None:
    print "The cutoff is: " + str(options.cutoff) + ' Angstrom.'
else:
    print "Using the default method to determine the bond exists."

#summary file
sf = open(options.sumf, 'w')
print >> sf, 'PDB_ID,', 'EXP_TECH,', 'RESOLUTION,', 'ATOM_NUMBER,', \
             'ION_NUMBER,', 'RES_ID,', 'RES_NAME,', 'ATOM_ID,', 'ATOM_NAME,', \
             'COORD_SPHERE,', 'GEOMETRY,','GEO_RMS'

#environment file
ef = open(options.envrmtf, 'w')
print >> ef, 'PDB,', 'ION_RESID,', 'ION_RESNAME,', 'ION_ATOM_ID,', \
             'ION_ATOM_NAME,', 'RESID,', 'RESNAME,', 'ATOM_ID,', 'ATOM_NAME,',\
             'DISTANCE,', 'GEOMETRY,', 'GEO_RMS,', 'COORDINATE_SPHERE,', \
             'EXP_TECH,', 'RESOLUTION'

#==============================================================================
# Do analysis for each pdb file
#==============================================================================

for fname in pdbfnl:
    print "Performing the " + fname + " file"

    #get the metal list
    mol, atids, resids = get_atominfo_fpdb(fname)

    #Get the resolution and method
    fp1 = open(fname, 'r')
    for line in fp1:
        if 'RESOLUTION.' in line:
            line = line.split()
            try:
                reso = float(line[-1])
            except:
                try:
                    reso = float(line[-2])
                except:
                    reso = 'UNKNOWN'
        elif 'EXPERIMENT TYPE' in line:
            line = line.split()
            exptyp = line[-1]
            if line[-1] == 'DIFFRACTION' and line[-2] == 'X-RAY':
               exptyp = 'X-RAY'
    fp1.close()

    #Get the metal ion which is the ion user want to process
    metallist = []
    for i in atids:
        resname = mol.residues[mol.atoms[i].resid].resname
        atname = mol.atoms[i].atname
        if (resname, atname) in Metalpdb.keys():
            if Metalpdb[(resname, atname)] == ionname:
                metallist.append(i)

    #for each metal ion in the metal list, print the metal center
    for i in metallist:

        mccrds = [] #The crds of metal site
        crdi = mol.atoms[i].crd
        elmti = mol.atoms[i].element
        residi = mol.atoms[i].resid
        atnamei = mol.atoms[i].atname
        resnamei = mol.residues[residi].resname
        radiusi = CoRadiiDict[elmti]
        mcresids = [] #MetalCenter residue IDs

        #Get the residues which is the metal site
        for j in atids:
            if j != i:
                atnamej = mol.atoms[j].atname
                crdj = mol.atoms[j].crd
                residj = mol.atoms[j].resid
                resnamej = mol.residues[residj].resname
                elmtj = mol.atoms[j].element
                radiusj = CoRadiiDict[elmtj]
                radiusij = radiusi + radiusj + 0.40
                disij = calc_bond(crdi, crdj)

                if options.cutoff == None:
                    if (disij >= 0.1) and (disij <= radiusij) \
                       and (elmtj != 'H'):
                        mccrds.append(crdi)
                        mccrds.append(crdj)
                        if (residj not in mcresids):
                            mcresids.append(residj)
                else:
                    if (disij >= 0.1) and (disij <= options.cutoff) \
                       and (elmtj != 'H'):
                        mccrds.append(crdi)
                        mccrds.append(crdj)
                        if (residj not in mcresids):
                            mcresids.append(residj)

        #Getting the ligating reidue letters
        reslets = ''
        for j in mcresids:
            resname = mol.residues[j].resname
            if resname in resdict.keys():
                reslet = resdict[resname]
            else:
                reslet = 'X'
            reslets = reslets + reslet
        nospace = ''
        reslets = nospace.join(sorted(reslets))
        print len(mcresids), reslets

        #Get the geometry and geometry rms
        geo, georms = det_geo(mccrds)

        #add the metal ions into the mcresids
        if mol.atoms[i].resid not in mcresids:
            mcresids.append(mol.atoms[i].resid)

        mcpdbfn = fname.strip('.pdb') + '_res_' + str(i) + '_MetalCenter.pdb'
        if os.path.isfile(mcpdbfn):
            print "Overwritting the metal center pdb file: " + mcpdbfn
            os.system("rm %s" %mcpdbfn)

        #print the residue which is in the cut off into the pdb file
        for j in mcresids:
            for k in mol.residues[j].resconter:
                tiker = mol.atoms[k].gtype
                atid = mol.atoms[k].atid
                atname = mol.atoms[k].atname
                resname = mol.atoms[k].resname
                chainid = 'A'
                resid = mol.atoms[k].resid
                crdx = round(mol.atoms[k].crd[0], 3)
                crdy = round(mol.atoms[k].crd[1], 3)
                crdz = round(mol.atoms[k].crd[2], 3)
                occp = 1.00
                tempfac = 0.00
                atmj = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmj, mcpdbfn)

        #Print the environment
        for j in atids:
            atnamej = mol.atoms[j].atname
            crdj = mol.atoms[j].crd
            residj = mol.atoms[j].resid
            resnamej = mol.residues[residj].resname
            elmtj = mol.atoms[j].element
            radiusj = CoRadiiDict[elmtj]
            radiusij = radiusi + radiusj + 0.40
            disij = calc_bond(crdi, crdj)

            if options.cutoff == None:
                if (disij >= 0.1) and (disij <= radiusij) \
                   and (elmtj != 'H'):
                    #print the environment file, for each bond in the metal site
                    print >> ef, fname.strip('.pdb'),',', residi,',', resnamei, \
                             ',', i, ',', atnamei,',', residj,',', resnamej,\
                             ',', j,',', atnamej, ',', round(disij, 3),',', geo,\
                             ',', round(georms, 3),',', reslets,',', exptyp,\
                             ',', reso
            else:
                if (disij >= 0.1) and (disij <= options.cutoff) \
                   and (elmtj != 'H'):
                    #print the environment file, for each bond in the metal site
                    print >> ef, fname.strip('.pdb'),',', residi,',', resnamei, \
                             ',', i, ',', atnamei,',', residj,',', resnamej,\
                             ',', j,',', atnamej, ',', round(disij, 3),',', geo,\
                             ',', round(georms, 3),',', reslets,',', exptyp,\
                             ',', reso

        #print the summary file, for each metal site
        print >> sf, fname.strip('.pdb'),',', exptyp,',', reso, \
                 ',', len(atids), ',', len(metallist),',', residi,\
                 ',', resnamei,',', i,',', atnamei,',', reslets,',', geo,\
                 ',', round(georms, 3)

sf.close()
ef.close()

