#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    version: 1
#
# To generate poseview diagram for a protein-ligand complex
#
import subprocess
import sys
import os
import numpy as np
from operator import itemgetter # sort list
import re #regular expression

from GenerateGraph.GenerateGraph import GenerateGraph

generateGraph = GenerateGraph()


solvent = "SOL"
path_asgard= repr(sys.argv[0])										
path_asgard= path_asgard[1:path_asgard.rfind("ASGARD/")+len("/ASGARD")]
poseview_run = path_asgard + "external_sw/poseview/poseview"
babel_run = path_asgard + "external_sw/babel/babel"




class Energy():
    def __init__(self, avarage_all, lst_all, last_all, add_coul, avarage_coul , last_coul, lst_coul, add_lj, avarage_lj, last_lj, lst_lj):
        self.avarage_all = avarage_all
        self.lst_all = lst_all
        self.last_all = last_all

        self.avarage_coul = avarage_coul
        self.last_coul = last_coul
        self.add_coul = add_coul
        self.lst_coul = lst_coul

        self.avarage_lj = avarage_lj
        self.last_lj = last_lj
        self.add_lj = add_lj
        self.lst_lj = lst_lj

def generate_poseview(pdb):
    prefix = '{}_{}_{}'.format(os.path.splitext(pdb)[0], mol_target_original_name, mol_query_original_name)
    cmd = 'cat {} | grep -v L[0-9][0-9] |grep -v Li[0-9] |grep -v NA |grep -v CL |grep -v {}  |grep -v TER |grep -v ENDMDL|grep -v UNL|grep -v UNK > {}'.format(pdb, solvent, prefix+"_target.pdb")
    subprocess.check_output(cmd, shell=True)
    cmd = 'cat {} | grep {} > {}'.format(pdb, mol_query_name, prefix+"_query.pdb" )
    subprocess.check_output(cmd, shell=True)
    cmd = '{} -ipdb {} -omol2 {} 2> /dev/null'.format(babel_run, prefix+"_target.pdb",  prefix+"_target.mol2")
    subprocess.check_output(cmd, shell=True)
    cmd = '{} -ipdb {} -omol2 {} 2> /dev/null'.format(babel_run,  prefix+"_query.pdb",prefix+"_query.mol2" )
    subprocess.check_output(cmd, shell=True)
    cmd = '{} -p {} -l {} -o {} -t  {}'.format(
        poseview_run,
        prefix+"_target.mol2",
        prefix+"_query.mol2",
        n_g_poseview,
        "\""+mol_target_original_name +" "+mol_query_original_name+  "\""
    )
    subprocess.check_output(cmd, shell=True)
    cmd = 'rm {}*'.format(prefix)

results_prefix  = sys.argv[2]									#graph_output

system_pdb = sys.argv[1]

n_g_poseview = results_prefix + "_poseview.png"

#	4 poseview diagram

generate_poseview(system_pdb)

