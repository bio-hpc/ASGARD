#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Recibe:
#       proteina, ligando
#   Genera:
#       Topologia para gormacs, cubo de aguas y hace una dm de 1 paso
#
import sys
import os
import subprocess
import re
import shutil


def is_tool(name):
    """Check whether `name` is on PATH."""
    from distutils.spawn import find_executable
    return find_executable(name) is not None


if is_tool("gmx_mpi"):
    GMX = "gmx_mpi"
elif is_tool("gmx"):
    GMX = "gmx"
else:
    print ("Error: Not found gmx_mpi or gmx")
    print ("")
    exit()
PYTHON_RUN = "python"
GENERATE_TOPOL_SCRIPT = " ASGARD/external_sw/gromacs/topology/generate_topology.py "
GENERATE_TOPOL = PYTHON_RUN + GENERATE_TOPOL_SCRIPT + '-t {} -q {} -p TARGET_ONE_QUERY -g '+GMX
PYTHON_EXE = "python"
SCRIPT_CONVERT_TO = "ASGARD/extra_shuttlemol/used_by_ASGARD/convert_to.py"
CONVERT_TO = '{} {} {} {} '
#TYPE_GRID = 'dodecahedron'
TYPE_GRID = 'cubic'
PADDING = 1
F_TPR_CONFIG = 'ASGARD/external_sw/gromacs/config_files/tpr.sh'
F_MD_CONFIG = 'ASGARD/external_sw/gromacs/config_files/simulation.sh'
MAX_WARNINGS = -1
CORES = ""  # "-nt 4" #sependiendo de si gromacs esta compilado con mpi o no este parametro debe estar vacio o con el numero de hilos
SOLVATION = "SOL"
STEPS_MIN = 100
STEPS_MD = 1
os.environ['GMX_MAXBACKUP'] = "-1"


def convert_molecule(file_1, file_pdb):
    cmd = CONVERT_TO.format(PYTHON_EXE, SCRIPT_CONVERT_TO, file_1, file_pdb)
    execute_cmd(cmd)


def execute_cmd(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.PIPE, shell=True, executable="/bin/bash")
    out, err = p.communicate()
    if "error" in str(out).lower():
        print(out)
    if "error" in str(err).lower():
        print(err)
    return out


if len(sys.argv) != 3:
    print("\nERROR: Parameters:")
    print("1º Receptor ")
    print("1º Ligand ")
    exit()
f_bd_rec = sys.argv[1]
f_bd_lig = sys.argv[2]
b_bd_rec, e_bd_rec = os.path.splitext(f_bd_rec)
b_bd_lig, e_bd_lig = os.path.splitext(f_bd_lig)

new_target = os.path.dirname(b_bd_lig)+"/"+os.path.basename(b_bd_lig).split("_")[0]+"_"+os.path.basename(b_bd_rec)+e_bd_rec
#
#      Modificacion para el script generico tanto AD com BD
#
# folder_exp =  "BD_"
# new_target = new_target[new_target.find(folder_exp)+1:] # si se ha introducido el path absoulto se divide por /BD
folder_exp = "ASGARD_analysis/"
new_target = new_target[new_target.rfind(folder_exp)+len("shuttlemol/"):] # si se ha introducido el path absoulto se divide por /BD
shutil.copyfile(f_bd_rec, new_target)
target = os.path.basename(os.path.splitext(new_target)[0])
#
#   Convertimos a formato
#
if e_bd_rec != ".pdb":  
    convert_molecule(new_target, os.path.splitext(new_target)[0] + ".pdb")
    f_bd_rec = os.path.splitext(new_target)[0] + ".pdb"
else:
    f_bd_rec = new_target
if e_bd_lig != ".mol2":
    convert_molecule(f_bd_lig, b_bd_lig + ".mol2")
    f_bd_lig = b_bd_lig + ".mol2"

#
#   Generamos topologia
#

cmd = GENERATE_TOPOL.format(f_bd_rec, f_bd_lig)
out = execute_cmd(cmd)
#for i in out.split("\n"):
#    print i

#
#   Generar grid y añadir solvent
#

prefix_in = os.path.dirname(b_bd_lig)+"/"+target+"_"+os.path.basename(b_bd_lig)

charge = float(subprocess.check_output(['tail', '-1', prefix_in+"_complex.eng"]).split(" ")[1].strip())
cmd = '{} editconf -f {} -o {} -d {} -bt {} -c'.format(
    GMX,
    prefix_in + "_complex.gro",
    prefix_in + "_complex_box.gro",
    PADDING,
    TYPE_GRID
)

execute_cmd(cmd)

cmd = '{} solvate -cp {} -o {} -p {}'.format(
    GMX,
    prefix_in + "_complex_box.gro",
    prefix_in + "_complex_solv.gro",
    prefix_in + "_complex.top"
)
execute_cmd(cmd)


#
#
#   Añadimos iones
#
os.environ['padding_grid'] = str(PADDING)
os.environ['file_conf_tpr'] =  prefix_in+"_tpr.mdp"
cmd = 'source {}'.format(F_TPR_CONFIG)
execute_cmd(cmd)
cmd = '{0} grompp -f  {1} -c {2} -p {3} -o {4} -po {5} -maxwarn {6}'.format(
    GMX,
    prefix_in+"_tpr.mdp",
    prefix_in + "_complex_solv.gro",
    prefix_in + "_complex.top",
    prefix_in + "_complex_ions.tpr",
    prefix_in+"_tpr_po.mdp",
    MAX_WARNINGS
)
#print cmd
#exit()
execute_cmd(cmd)
charge = int(charge)
if charge > 0:
    ions = "-nname CL  -nn "+str(abs(int(charge))) 
elif charge < 0:
    ions = "-pname NA  -np "+str(abs(int(charge)))
else:
    ions = ""
cmd = 'echo {} | {} genion -s {} -o {} -p {} {}'.format(
    SOLVATION,
    GMX,
    prefix_in + "_complex_ions.tpr",
    prefix_in + "_complex_solv_ions.gro",
    prefix_in + "_complex.top",
    ions
)
execute_cmd(cmd)



#
#   Minimizacion apra relajar el sistema
#

os.environ['step_min'] = str(STEPS_MIN)
os.environ['padding_grid'] = str(PADDING)
os.environ['file_conf_min'] = prefix_in +"_min.mdp"
F_TPR_CONFIG = 'ShuttleMol/external_sw/gromacs/config_files/minimization.sh'

cmd = 'source {}'.format(F_TPR_CONFIG)
execute_cmd(cmd)
cmd = '{} grompp -f {} -c {} -p {} -o {} -maxwarn {}'.format(
    GMX,
    os.environ['file_conf_min'],
    prefix_in + "_complex_solv_ions.gro",
    prefix_in + "_complex.top",
    prefix_in + "_complex_min.tpr",
    MAX_WARNINGS
)

execute_cmd(cmd)
cmd = '{} mdrun -v -deffnm {} {} -g {}'.format(
    GMX,
    prefix_in +"_complex_min",
    CORES,
    prefix_in +"_min.log",
)
execute_cmd(cmd)

#
#   Simulation DM
#
#   Configuracion
os.environ['file_conf_md'] = prefix_in + "_complex_md.mdp"
os.environ['int_md'] = "0"
os.environ['step_md'] = str(STEPS_MD)
os.environ['write_data'] = "1"
os.environ['nstlist'] = "20"
os.environ['padding_grid'] = str(PADDING)
os.environ['tc_grps'] = "Protein Non-Protein"
os.environ['tau_t'] = "0.1 0.1"
os.environ['ref_t'] = "300 300"
os.environ['comm_grps'] = "Protein"
cmd = 'source {}'.format(F_MD_CONFIG)
execute_cmd(cmd)
#
#   Ejecucion
#
cmd = '{0} grompp -f {1} -c {2}  -p {3} -o {4} -maxwarn {5}'.format(
    GMX,
    prefix_in + "_complex_md.mdp",
    prefix_in + "_complex_min.gro", #"_complex_solv_ions.gro",
    prefix_in + "_complex.top",
    prefix_in + "_complex_md.tpr",
    MAX_WARNINGS
)
execute_cmd(cmd)
f = open (prefix_in + "_complex_md.mdp",'a')
f.write("energygrps = Protein L01")
f.close()
cmd = '{} mdrun -deffnm  {} {} -g {}'.format(
    GMX,
    prefix_in + "_complex_md",
    CORES,
    prefix_in + "_complex_md.log"
)
execute_cmd(cmd)
#
# Rerun
#
cmd = '{0} grompp -f {1} -c {2}  -p {3} -o {4} -maxwarn {5}'.format(
    GMX,
    prefix_in + "_complex_md.mdp",
    prefix_in + "_complex_md.gro", #"_complex_solv_ions.gro",
    prefix_in + "_complex.top",
    prefix_in + "_complex_md_re.tpr",
    MAX_WARNINGS
)
execute_cmd(cmd)
cmd = '{} mdrun -deffnm  {} {} -g {} -rerun {}'.format(
    GMX,
    prefix_in + "_complex_md_re",
    CORES,
    prefix_in + "_complex_md_re.log",
    prefix_in + "_complex_md.xtc"
)
execute_cmd(cmd)
#
#   Traemos energia
#
groups ="Coul-SR:Protein-L01 \n Coul-14:Protein-L01  \n LJ-SR:Protein-L01 \n LJ-14:Protein-L01 \n\n"
cmd = "echo \"{}\" > {}".format(
    groups,
    prefix_in + "_complex_md_re_groups.txt"

)
execute_cmd(cmd)
cmd = '{} energy -f {} -s {} -o {} -sum < {}'.format(
    GMX,
    prefix_in + "_complex_md_re.edr",
    prefix_in + "_complex_md_re.tpr",
    prefix_in + "_complex_md_re.xvg",
    prefix_in + "_complex_md_re_groups.txt"

)

out = execute_cmd(cmd)
#print out
for line in out.split("\n"):
    if str(line).startswith("Total"):

        line = re.sub(' +',' ',line)
        score = line.split(" ")[1]          #score global de la prueba
#
#   Borramos ficheros temporal
#
cmd = 'rm {} '.format(prefix_in + "*")
execute_cmd(cmd)


print("Score: "+score)
