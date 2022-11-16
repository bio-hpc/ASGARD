#!/usr/bin/env python
# Filename: MCPB.py
"""
This is the MCPB.py program written by Pengfei Li in Merz Research Group,
Michigan State University. All Rights Reserved. It is a re-written python
program of MCPB in MTK++. It is written to assist the metal center modeling
in mixed systems (especially the protein system). It optimize the workflow
of the MCPB and has better supports of different ions and force feilds. It
supports modeling of more than 50 ions from +1 to +4 oxidation states with
modeling by bonded and nonbonded model. It supports a series of AMBER force
fields (ff94, ff99, ff99SB, ff03, ff03.r1, ff10, ff12SB, ff14SB). The program
is not gurantee to work due to bugs may exist. Suggestions and bug reports
are welcome to send to Pengfei Li (Email address: ldsoar1990@gmail.com).

Please cite the following papers if you use the software to perform the
modeling:

The parameterization scheme is come from:
** M. B. Peters, Y. Yang, B. Wang, L. Fusti-Molnar, M. N. Weaver, K. M. Merz,
   JCTC, 2010, 6, 2935-2947

The Seminario method is from:
** J. M. Seminario IJQC, 1996, 30, 1271-1277

The Emprical method is from:
** P. Li, K. M. Merz, In preparation.

The RESP fitting radii, VDW parameters and 12-6-4 parameter sets of +1 metal
ions and halide ions are come from:
** P. Li, L. F. Song, K. M. Merz, JCTC, Accepted

The RESP fitting radii and VDW parameters of +2 metal ions are come from:
** P. Li, B. P. Roberts, D. K. Chakravorty, K. M. Merz, JCTC, 2013, 9,
   2733-2748

The 12-6-4 parameter set of +2 metal ions are from:
** P. Li, K. M. Merz, JCTC, 2014, 10, 289-297

The RESP fitting radii, VDW parameters and 12-6-4 parameter sets of +3 and
+4 metal ions are come from:
** P. Li, L. F. Song, K. M. Merz, JPCB, 2015, 119, 883-895
"""
#==============================================================================
# Load the MCPB module
#==============================================================================
from mcpb.gene_model_files import get_ms_resnames, gene_model_files
from mcpb.resp_fitting import resp_fitting
from mcpb.gene_pre_frcmod_file import gene_pre_frcmod_file
from mcpb.gene_final_frcmod_file import (gene_by_empirical_way,
            gene_by_QM_fitting_sem, gene_by_QM_fitting_zmatrix)
from mcpb.amber_modeling import gene_leaprc
from mcpb.title import print_title
from pymsmtmol.element import resnamel
from pymsmtexp import *
import warnings
import os
from optparse import OptionParser

parser = OptionParser("usage: -i input_file -s/--step step_number")
parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-s", "--step", dest="step", type='string',
                  help="Step number")
(options, args) = parser.parse_args()

#==============================================================================
# Get the input variables
#==============================================================================
# Print the title of the program
print_title('MCPB.py')
options.step = options.step.lower()

# Default values
cutoff = 2.8
ionchgfix = 0
naamol2fs = []
ff_choice = 'ff14SB'
gaff = 1
frcmodfs = []
gname = 'MOL'
g0x = 'g03'
ioninfo = []
sqmopt = 0
watermodel = 'tip3p'
paraset = 'cm'

if options.step not in ['1', '1n', '1m', '1a', '2', '2e', '2s', '2z',
                        '3', '3a', '3b', '3c', '3d', '4', '4b', '4n1',
                       '4n2']:
    raise pymsmtError('Invalid step number chosen. please choose among the '
                      'following values: 1, 1n, 1m, 1a, 2, 2e, 2s, 2z, 3, '
                      '3a, 3b, 3c, 3d, 4, 4b, 4n1, 4n2')

inputf = open(options.inputfile, 'r')
for line in inputf:
    line = line.split()
    if '\n' in line:
        line.remove('\n')
    if ',' in line:
        line.remove(',')
    if '' in line:
        line.remove('')
    if ' ' in line:
        line.remove(' ')
    if ':' in line:
        line.remove(':')
    #comment
    if line[0][0] == '#':
        continue
    #orpdbf
    if line[0].lower() == 'original_pdb':
        if len(line) == 2:
            orpdbf = line[1]
            if os.path.exists(orpdbf):
                continue
            else:
                raise pymsmtError('File %s does not exists.' %orpdbf)
        else:
            raise pymsmtError('%d pdb files provided for original_pdb, only '
                              'need one.' %(len(line)-1))
    #fname
    elif line[0].lower() == 'group_name' :
        if len(line) == 2:
            gname = line[1]
        elif len(line) == 1:
            warnings.warn('None group_name provided in the input file, '
                          'using %s as default.' %gname, pymsmtWarning)
        else:
            raise pymsmtError('More than one group_name provided, need one.')
    #cutoff
    elif line[0].lower() == 'cut_off':
        if len(line) == 2:
            try:
                cutoff = float(line[1])
            except:
                raise pymsmtError('Please provide an float number for the '
                                  'cut_off parameter.')
        elif len(line) == 1:
            warnings.warn('No cut_off parameter provided, Default value '
                          '%5.1f is used.' %cutoff, pymsmtWarning)
        else:
            raise pymsmtError('More than one cut_off values are provided, '
                              'need one.')
    #ionids
    elif line[0].lower() == 'ion_ids':
        if len(line) >= 2:
            try:
                ionids = line[1:]
                ionids = [int(i) for i in ionids]
            except:
                raise pymsmtError('ion_ids need to be integer number(s).')
        else:
            raise pymsmtError('ion_ids need to be provided.')
    #ioninfo
    elif line[0].lower() == 'ion_info':
        if (len(line)-1)%4 == 0:
            try:
                line[4::4] = [int(i) for i in line[4::4]]
                ioninfo = line[1:]
            except:
                raise pymsmtError('The charge of the ion in the ion_info '
                                  'should be integer number.')
        else:
            raise pymsmtError('Wrong ion_info format should be: residue name, '
                              'atom name, element, charge.')
    #ionmol2fs
    elif line[0].lower() == 'ion_mol2files':
        if len(line) >= 2:
            ionmol2fs = line[1:]
            for i in ionmol2fs:
                if os.path.exists(i):
                    continue
                else:
                    raise pymsmtError('File %s does not exists.' %i)
        else:
            raise pymsmtError('Need to provide the mol2 file(s) for '
                              'ion_mol2files.')
    #ionchgfix
    elif line[0].lower() == 'ionchg_fixation':
        if len(line) == 2:
            try:
                ionchgfix = int(line[1])
                if ionchgfix not in [0, 1]:
                    raise pymsmtError('Please provide 0 or 1 for the '
                                      'ionchg_fixation option.')
            except:
                raise pymsmtError('Please provide a integral number for the '
                                  'ionchg_fixation parameter.')
        elif len(line) == 1:
            warnings.warn('No ionchg_fixation parameter provided. '
                          'Default value %d will be used.'
                          %ionchgfix, pymsmtWarning)
        else:
            raise pymsmtError('More than one ionchg_fixation options '
                              'provided.')
    #naamol2fs
    elif line[0].lower() == 'naa_mol2files':
        if len(line) >= 2:
            naamol2fs = line[1:]
            for i in naamol2fs:
                if os.path.exists(i):
                    continue
                else:
                    raise pymsmtError('File %s does not exists.' %i)
        else:
            warnings.warn('No mol2 file is provided for '
                          'naa_mol2files.', pymsmtWarning)
    #gau_version
    elif line[0].lower() == 'gau_version':
        if len(line) == 2:
            g0x = line[1].lower()
            if g0x not in ['g03', 'g09']:
                raise pymsmtError('Please use either g03 or g09, other '
                                  'versions are not gurantee to support.')
        elif len(line) == 1:
            warnings.warn('No g0x parameter provided. Default value '
                          '%s is used.' %g0x, pymsmtWarning)
        else:
            raise pymsmtError('More than one gau_version parameters provided, '
                              'need one.')
    #sqmopt
    elif line[0].lower() == 'sqm_opt':
        if len(line) == 2:
            try:
                sqmopt = int(line[1])
                if sqmopt not in [0, 1, 2, 3]:
                    raise pymsmtError('sqm_opt varible needs to be 0, 1, 2 or '
                                      '3, 0 means not using, 1 means only do '
                                      'optimization for sidechain model. 2 '
                                      'means only do optimization for large '
                                      'model. 3 means do optimization for '
                                      'both models.')
            except:
                raise pymsmtError('sqm_opt value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No sqm_opt parameter provided. '
                                'Default value %d is used.' %sqmopt, pymsmtWarning)
        else:
            raise pymsmtError('More than one sqm_opt parameter provided, '
                              'need one.')
    #ff
    elif line[0].lower() == 'force_field':
        if len(line) == 2:
            ff_choice = line[1]
            if ff_choice not in ['ff94', 'ff99', 'ff99SB', 'ff03', 'ff03.r1',
                'ff10', 'ff12SB', 'ff14SB']:
                raise pymsmtError('Not support %s force field in current '
                                  'version. Only support ff94, ff99, ff99SB, '
                                  'ff03, ff03.r1, ff10, ff12SB, ff14SB.')
        elif len(line) == 1:
            warnings.warn('No force_field parameter provided. '
                          'Default value %s is used.'
                          %ff_choice, pymsmtWarning)
        else:
            raise pymsmtError('More than one force_field parameter provided, '
                            'need one: ff94, ff99, ff99SB, ff03, ff03.r1, '
                            'ff10, ff12SB or ff14SB.')
    #gaff
    elif line[0].lower() == 'gaff':
        if len(line) == 2:
            try:
                gaff = int(line[1])
                if gaff not in [0, 1]:
                    raise pymsmtError('gaff varible needs to be 0 or 1, '
                                      '0 means not using, 1 means using.')
            except:
                raise pymsmtError('gaff value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No gaff parameter provided. '
                          'Default value %d is used.'
                          %gaff, pymsmtWarning)
        else:
            raise pymsmtError('More than one gaff parameter provided, '
                              'need one.')
    #frcmodfs
    elif line[0].lower() == 'frcmod_files':
        if len(line) >= 2:
            frcmodfs = line[1:]
            for i in frcmodfs:
                if os.path.exists(i):
                    continue
                else:
                    raise pymsmtError('File %s does not exists.' %i)
        else:
            warnings.warn('No frcmod files is provided for '
                          'frcmod_files.', pymsmtWarning)
    #watermodel
    elif line[0].lower() == 'water_model':
        if len(line) == 2:
            watermodel = line[1].lower()
            if watermodel not in ['tip3p', 'spce', 'tip4pew']:
                raise pymsmtError('Not support %s water model. Only support '
                                  'TIP3P, SPCE, TIP4PEW' %watermodel.upper())
        elif len(line) == 1:
            warnings.warn('No water_model parameter provided. Default '
                          'value %s is used.'
                          %watermodel.upper(), pymsmtWarning)
        else:
            raise pymsmtError('More than one water_model parameter provided, '
                              'only need one: TIP3P, SPCE or TIP4PEW.')
    #paraset
    elif line[0].lower() == 'ion_paraset':
        if len(line) == 2:
            paraset = line[1].lower()
            if paraset not in ['hfe', 'cm', 'iod', '12_6_4']:
                raise pymsmtError('Do not have %s ion parameter set. Only '
                                  'have HFE, CM, IOD, 12_6_4 parameter set'
                                  %paraset.upper())
        elif len(line) == 1:
            warnings.warn('No ion_paraset parameter provided. '
                          'Default value %s is used.'
                          %paraset.upper(), pymsmtWarning)
        else:
            raise pymsmtError('More than one ion_paraset parameter provided, '
                              'only need one: HFE, CM, IOD or 12_6_4.')
inputf.close()

print "The input file you are using is : %s" %options.inputfile
print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

#Print the input variables
print "The following is the input variable you have:"

try:
    print 'The variable original_pdb is : ', orpdbf
except:
    raise pymsmtError('original_pdb needs to be provided.')

try:
    print 'The variable ion_ids is : ', ionids
except:
    raise pymsmtError('ion_ids needs to be provided.')

try:
    print 'The variable ion_mol2files is : ', ionmol2fs
except:
    raise pymsmtError('ion_mol2files needs to be provided.')

print 'The variable group_name is : ', gname
print 'The variable cut_off is : ', cutoff
print 'The variable ionchg_fixation is : ', ionchgfix
print 'The variable gau_version is : ', g0x
print 'The variable sqm_opt is : ', sqmopt
print 'The variable force_field is : ', ff_choice
print 'The variable gaff is : ', gaff
print 'The variable frcmodfs is : ', frcmodfs
print 'The variable naa_mol2files is : ', naamol2fs
print 'The variable water_model is : ', watermodel.upper()
print 'The variable ion_paraset is : ', paraset.upper(), "(Only for nonbonded model)"

if options.step in ['4n2']:
    if ioninfo == []:
        raise pymsmtError('The variable ion_info need to be provided in step '
                          '%s.' %options.step)
else:
    print 'The variable ion_info is : ', ioninfo

#==============================================================================
# Related define
#==============================================================================
#Get the renamed residue name
mcresname0, mcresname = get_ms_resnames(orpdbf, ionids, cutoff)
for i in mcresname0:
    if (i not in resnamel) and (i+'.mol2' not in naamol2fs):
        raise pymsmtError('%s is required in naa_mol2files but not '
                          'provided.' %i)

#The mol2 file used in pre-generated model
premol2fs = ionmol2fs + naamol2fs

##pdb files
scpdbf = gname + '_sidechain.pdb'
scpdbf2 = gname + '_sqm.pdb'
stpdbf = gname + '_standard.pdb'
lgpdbf = gname + '_large.pdb'
fipdbf = gname + '_mcpbpy.pdb'

##finger print files
stfpf = gname + '_standard.fingerprint'
lgfpf = gname + '_large.fingerprint'

##frcmod files
prefcdf = gname + '_pre.frcmod'
finfcdf = gname + '.frcmod'

##log file
fclogf = gname + '_sidechain_fc.log'
mklogf = gname + '_large_mk.log'

##checkpoint file
fcfchkf = gname + '_sidechain_opt.fchk'

##tleap input file
ileapf = gname + '_tleap.in'

#==============================================================================
# Step 1 General_modeling
#==============================================================================
#1. Generate the modeling files:
#Pdb files for sidechain, standard and large model
#Gaussian input files for sidechain, large model
#Fingerprint files for standard and large model
#Three options:
#1n) Don't rename any of the atom types in the sidechain fingerprint file
#1m) Just rename the metal ion to the AMBER ion atom type style
#1a) Default. Automatically rename the atom type of the atoms in the metal
#    complex.
if (options.step == '1'):
    gene_model_files(orpdbf, ionids, gname, ff_choice, premol2fs, cutoff,
                     watermodel, 2, sqmopt)
elif (options.step == '1n'):
    gene_model_files(orpdbf, ionids, gname, ff_choice, premol2fs, cutoff,
                     watermodel, 0, sqmopt)
elif (options.step == '1m'):
    gene_model_files(orpdbf, ionids, gname, ff_choice, premol2fs, cutoff,
                     watermodel, 1, sqmopt)
elif (options.step == '1a'):
    gene_model_files(orpdbf, ionids, gname, ff_choice, premol2fs, cutoff,
                     watermodel, 2, sqmopt)
#==============================================================================
# Step 2 Frcmod file generation
#==============================================================================
#Mass, dihedral, improper, VDW and metal ion non-related bond and metal
#ion non-related angle parameters are generated first. While the metal ion
#related bond and angle parameters are generated later while they could
#generated by using different methods.
#2e) Empirical method developed by Pengfei Li and co-workers in Merz group
#2s) Default. Seminario method developed by Seminario in 1990s
#2z) Z-matrix method
elif (options.step == '2'): #Default
    gene_pre_frcmod_file(ionids, premol2fs, stpdbf, stfpf, prefcdf, ff_choice,
                         gaff, frcmodfs, watermodel)
    gene_by_QM_fitting_sem(scpdbf, ionids, stfpf, prefcdf, finfcdf, fcfchkf,
                         g0x)
elif (options.step == '2e'):
    gene_pre_frcmod_file(ionids, premol2fs, stpdbf, stfpf, prefcdf, ff_choice,
                         gaff, frcmodfs, watermodel)
    gene_by_empirical_way(scpdbf, ionids, stfpf, prefcdf, finfcdf)
elif (options.step == '2s'):
    gene_pre_frcmod_file(ionids, premol2fs, stpdbf, stfpf, prefcdf, ff_choice,
                         gaff, frcmodfs, watermodel)
    gene_by_QM_fitting_sem(scpdbf, ionids, stfpf, prefcdf, finfcdf, fcfchkf,
                         g0x)
elif (options.step == '2z'):
    gene_pre_frcmod_file(ionids, premol2fs, stpdbf, stfpf, prefcdf, ff_choice,
                         gaff, frcmodfs, watermodel)
    gene_by_QM_fitting_zmatrix(scpdbf, ionids, stfpf, prefcdf, finfcdf, fclogf)
#==============================================================================
# Step 3 Doing the RESP charge fitting and generate the mol2 files
#==============================================================================
#3. Generate mol2 files with the charge parameters after resp charge fitting
#3a) All all the charges of the ligating residues could change
#3b) Default. Restrains the charges of backbone heavy atoms according to the
#    force field chosen
#3c) Restrains the charges of backbone atoms according to the force field
#    chosen
#3d) Restrains the charges of backbone atoms and CB atom in the sidechain
#    according to force field chosen
elif (options.step == '3'): #Default
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
                 premol2fs, mcresname, 1, ionchgfix)
elif (options.step == '3a'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
                 premol2fs, mcresname, 0, ionchgfix)
elif (options.step == '3b'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
                 premol2fs, mcresname, 1, ionchgfix)
elif (options.step == '3c'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
                 premol2fs, mcresname, 2, ionchgfix)
elif (options.step == '3d'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
                 premol2fs, mcresname, 3, ionchgfix)
#==============================================================================
# Step 4 Prepare the modeling file for leap
#==============================================================================
#4. Prepare the final modeling file
#4b) Default. Bonded model
#4n1) Nonbonded model with refitting the charge in the protein complex
#4n2) Normal Nonbonded model (12-6 nonbonded model) without re-fitting charges
elif (options.step == '4'): #Default
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
                ioninfo, mcresname, naamol2fs, ff_choice, frcmodfs, finfcdf,
                ileapf, 1, watermodel, paraset)
elif (options.step == '4b'): #bonded model
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
                ioninfo, mcresname, naamol2fs, ff_choice, frcmodfs, finfcdf,
                ileapf, 1, watermodel, paraset)
elif (options.step == '4n1'): #nonbonded model with refitting the charge
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
                ioninfo, mcresname, naamol2fs, ff_choice, frcmodfs, finfcdf,
                ileapf, 2, watermodel, paraset)
elif (options.step == '4n2'): #normal nonbonded model
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
                ioninfo, mcresname, naamol2fs, ff_choice, frcmodfs, finfcdf,
                ileapf, 3, watermodel, paraset)

quit()

