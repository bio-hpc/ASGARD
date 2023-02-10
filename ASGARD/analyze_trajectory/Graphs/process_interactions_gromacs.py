#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
import subprocess
import sys
import os
import re
os.environ['GMX_MAXBACKUP'] = "-1"
NUM_GROUP_ENERGIES = 4
SPACES = " "
UNDERSCORE = "_"


def add_lst(res1, res2, lst_res, j, aux):
    if res1 != "rest" and res2 != "rest":
        if res1 != res2 and res1 != name_solvent and res2 != name_solvent:
            if res1 == mol_query_name or res2 == mol_query_name:
                if res1 != res2:
                    lst_res.append(aux[j] + " " + aux[j + 1])


def get_intertactions_residues(lst_in, mode):
    """
       mode = true protein ligand (total)
       mode = false residue ligand
    """
    token = "-------------------------------------------------"
    lst_res = []
    cnt_token = 0
    for line in lst_in:
        if token in line:
            cnt_token += 1
        if cnt_token == 1:
            if mol_query_name in line:
                aux = re.sub(' +', ' ', line.strip()).split(" ")
                for j in range(0, len(aux), 2):
                    if ":" in aux[j + 1]:
                        aux_residue = aux[j + 1].split(":")[1]
                        res1 = aux_residue.split("-")[0]
                        res2 = aux_residue.split("-")[1]
                        #if mode and (res1 == "Protein" or res1 == mol_query_name) and (res2 == "Protein" or res2 == mol_query_name):
                        if mode and (res1 == mol_target_name or res1 == mol_query_name) and (res2 == mol_target_name or res2 == mol_query_name):
                            add_lst(res1, res2, lst_res, j, aux)
                        elif not mode:
                            add_lst(res1, res2, lst_res, j, aux)
    return lst_res


def create_energies_xvg( edr_rerun, tpr_min, g_energy, folder_hbonds, mode):
    """
        Genera ficheros xvg de energias grupo y global
    """
    cmd = 'echo 0| {} -f {} -s {}'.format(g_energy, edr_rerun, tpr_min)
    output = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True).stdout.read()
    lst = get_intertactions_residues(output.split("\n"), mode)
    for i in range(0, len(lst), 4):
        aux = lst[i].split(" ")[0], lst[i + 1].split(" ")[0],  lst[i + 2].split(" ")[0], lst[i + 3].split(" ")[0]

        title_no_group = '{}_no_group_{}.xvg'.format(UNDERSCORE.join(aux), lst[i].split(":")[1].replace("-", "_"))
        cmd = 'echo "{} 0" | {} -f {} -o {}'.format(SPACES.join(aux), g_energy, edr_rerun, folder_hbonds + title_no_group)
        subprocess.check_call(cmd, shell=True)

#
#   Only using GetResults
#
if len(sys.argv) != 18:
    print("Parameters: ")
    print("\t1º\t self.cfg.gromacs")
    print("\t2º\t self.cfg.graph")
    print("\t3º\t self.cfg.mpi")
    print("\t4º\t self.cfg.prefix_rerun_md")
    print("\t5º\t self.cfg.prefix_rerun_query")
    print("\t6º\t self.cfg.tpr_min")
    print("\t7º\t folder_xvg")
    print("\t8º\t self.cfg.name_solvent")
    print("\t9º\t self.cfg.graphEnergiesSum")
    print("\t10º\t self.cfg.prefix_results")
    print("\t11º\t self.cfg.mol_query")
    print("\t12º\t self.cfg.mol_query_name_original")
    print("\t13º\t self.cfg.mol_target")
    print("\t14º\t self.cfg.mol_target_name_original")
    print("\t15º\t self.cfg.grapGEnergy")
    print("\t16º\t self.cfg.pdb")
    print("\t17º\t self.cfg.ENERGIES_DISCARD")
    print("\t"+str(len(sys.argv)))
    exit()
gromacs             = sys.argv[1]
graph               = sys.argv[2]
mpi                 = sys.argv[3]
prefix_rerun_md     = sys.argv[4]
prefix_rerun_query  = sys.argv[5]
tpr_min             = sys.argv[6]
folder_xvg          = sys.argv[7]
name_solvent        = sys.argv[8]
graph_interactions_gromacs = sys.argv[9]
prefix_results      = sys.argv[10]
mol_query_name      = sys.argv[11]
mol_query_original_name = sys.argv[12]
mol_target_name     = sys.argv[13]
mol_target_original_name = sys.argv[14]
g_energy            = sys.argv[15]
system_pdb          = sys.argv[16]
energies_discard    = sys.argv[17]

folder_xvg += "/"
os.chdir(folder_xvg)
python_run = "python"


#   1º Generate the interactions files between residue and target
create_energies_xvg(prefix_rerun_query + ".edr", tpr_min, g_energy, folder_xvg, False)
#   2º Generate the interactions files between protein and query
create_energies_xvg(prefix_rerun_md + ".edr", tpr_min, g_energy, folder_xvg, True)

cmd = '{} {}\\\n\t {}\\\n\t {}\\\n\t {}\\\n\t {}\\\n\t {}\\\n\t {}\\\n\t {}\\\n\t  {} '.format(
    python_run, graph_interactions_gromacs, folder_xvg, prefix_results, mol_target_name,
    mol_target_original_name, mol_query_name, mol_query_original_name, system_pdb, energies_discard)
#print (cmd)
#exit()
subprocess.check_output(cmd, shell=True)
