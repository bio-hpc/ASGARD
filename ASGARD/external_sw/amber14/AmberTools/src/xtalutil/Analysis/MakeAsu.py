#! /usr/bin/env python
import os
from numpy import *
from ReadAmberFiles import *
import argparse 

#======================================================================#
#                                                                      #
# Prepare asu prmtop, rst7, pdb from SC_REFERENCE rst7.
#                                                                      #
#======================================================================#


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--SC_Topology", help="Amber supercell topology")
parser.add_argument("-rst7", "--SC_Reference", help="Amber supercell rst7 with experimental coords")
parser.add_argument("-resid", "--ASU_Residues", help="Number of residues in asu")
args = parser.parse_args()

f=open('ctraj.make_asu.in','w')
f.write('parm %s\n' %args.SC_Topology)
f.write('trajin %s\n' %args.SC_Reference)
f.write('strip \':%d-999999\' outprefix asu\n' %(int(args.ASU_Residues)+1))
f.write('trajout asu.rst7 restart')
f.close()
os.system('cpptraj <ctraj.make_asu.in >>XtalAnalyze.log')
os.system('mv asu.%s asu.prmtop' %(os.path.basename(args.SC_Topology)))
os.system('ambpdb -p asu.prmtop <asu.rst7 >asu.pdb 2>XtalAnalyze.log')
