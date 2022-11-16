#!/usr/bin/env python
# -*- coding: utf-8 -*-
#       En Dm limpiar el complejo con 
#       cat  VS_GR_naCHTvsNLRC4CHAINL_modificado_91_complex_md.pdb |grep -v 'SOL\|NA\|ACE\|NME'
#

import glob
import sys
import os
import subprocess
import shutil
import json
from interactions_colors import COLORS_INTERACTIONS
from add_paths import *
add_paths()

PYTHON_RUN = "python "
SCRIPT_CONVERT_TO = "ShuttleMol/extra_shuttlemol/used_by_shuttlemol/convert_to.py"
path_ext_sw = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),'external_sw/')
img_ubuntu = os.path.join(path_ext_sw,'img_ubuntu','ubuntu.simg')
F_PREPARE_COMPLEX = "singularity exec "+img_ubuntu+ ' python /opt/prepare_complex.py {} {} {}'
F_EDIT_SESSION = "singularity exec "+img_ubuntu+ ' python /opt/edit_sessions.py {}'
CONVERT_TO = '{} {} {} {} -b'
#PLIP_CMD = "ShuttleMol/external_sw/plip/plipcmd"
F_PLIP_CMD =  "singularity exec "+img_ubuntu+ ' plipcmd  '+" -f {} -o {} -y -t --maxthreads {} "
#PLIP_RUN = PLIP_CMD+" -f {} -o {} -y -t --maxthreads {}"
CORES = 1


def help():
    print("Error: Parameters")
    print("1º [ Receptor | Complex ] ")
    print("2º [ Ligand  | None ]")
    print("3º Output prefix")
    exit()

def execute_cmd(cmd):
    print (cmd)
    return subprocess.check_output(cmd, shell=True)

def run_plip(out_prefix):
    cmd = 'sed -i -e "s/\*\*\*/LIG/g" '+ out_prefix + '_interactions_complex.pdb'
    subprocess.check_output(cmd, shell=True)

    cmd = 'sed -i -e "s/>/ /" '+ out_prefix + '_interactions_complex.pdb'
    subprocess.check_output(cmd, shell=True)

    cmd = 'sed -i -e "s/</ /" '+ out_prefix + '_interactions_complex.pdb'
    subprocess.check_output(cmd, shell=True)    
    
    
    cmd = F_PLIP_CMD.format(
        out_prefix + "_interactions_complex.pdb",
        out_prefix + "_complex",
        CORES
    )
    # print (cmd)
    # exit()
    try:
        execute_cmd(cmd)
        # subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        # raise subprocess.CalledProcessError
        # plip lo manda todo a una carpeta, lo movemos a donde nos interesa y borramos esa carpeta
        pattern = glob.glob(out_prefix + '_complex/*.pse')
        fname, _ = os.path.splitext(pattern[0])

        #for ext in ('.pse', '.pml'):
        #    shutil.move(fname + ext, out_prefix + '_interactions' + ext)
        ext = ".pse"
        shutil.move(fname + ext, out_prefix + '_interactions' + ext)

        # El report va con otro nombre
        shutil.move(out_prefix + '_complex/report.txt', out_prefix + '_interactions.txt')

        # Eliminar el directorio de trabajo
        #shutil.rmtree(out_prefix + '_complex')
    except Exception as e:
        print("ERROR: Plip", out_prefix + '_complex')
        print (e)


def convert_molecule(file_1, file_pdb):
    cmd = CONVERT_TO.format(PYTHON_RUN, SCRIPT_CONVERT_TO, file_1, file_pdb)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)


def convert_name_plip(name_plip):
    return name_plip.replace("**", "").replace(" ", "").replace("-", "")


def read_file_plip_txt( prefix_interactions_plip ):
    """
        1º
        retorn:
        dicionario:
            clave = typo energia
            dict[key]][0]tipo de columna las validas son: 1 2 3 y penultima coords lig ultima coords prot
            dict[key]][1] scores
            dict[key]][2]
            ...
        2º Guarda las energias en un fichero json
    """
    file_txt = prefix_interactions_plip + '.txt'
    file_json = prefix_interactions_plip + '.json'
    table_scores = {}
    key = ""
    start_token = "+=======+"
    mid_token = "+-------+"
    title_token = "**"

    table_scores['interactions_shape'] = COLORS_INTERACTIONS
    table_scores['interactions_groups'] = {}
    if os.path.isfile(file_txt):
        f = open(file_txt)
        for line in f:
            if title_token in line:
                line = line.replace(title_token, "").strip()
                key = convert_name_plip(line)
            if line.startswith("|") and mid_token not in line and start_token not in line:
                line = line.replace(' ', '').strip()[1:-1]
                if key not in table_scores['interactions_groups']:
                    table_scores['interactions_groups'][key] = {}
                    table_scores['interactions_groups'][key]['interactions'] = []
                if len(table_scores['interactions_groups'][key]) == 1:
                    table_scores['interactions_groups'][key]['legend'] = line
                elif line != table_scores['interactions_groups'][key]['legend']:
                    table_scores['interactions_groups'][key]['interactions'].append(line)

        f.close()

    parsed = json.loads( json.dumps(table_scores))

    with open(file_json, 'w') as outfile:
        json.dump(parsed, outfile,indent=4, sort_keys=True)

    return table_scores


if __name__ == "__main__":
    if len(sys.argv) != 4:
        help()

    rec = sys.argv[1]
    lig = sys.argv[2]
    out_prefix = sys.argv[3]
    if lig == "None":
        lig = None

    if lig != None and not os.path.isfile(lig) or not os.path.isfile(rec):
        print ("ERROR: ligando o recptor no existe")
        exit()
    rec_name, rec_ext = os.path.splitext(rec)

    
    if rec_ext != ".pdb":
        if not os.path.isfile(rec_name+".pdb"):        
            convert_molecule(rec, rec_name + ".pdb")
        rec = rec_name + ".pdb"
    if lig == None:
        shutil.copy(rec, out_prefix + "_interactions_complex.pdb")
    else:
        lig_name, lig_ext = os.path.splitext(lig)
        if lig_ext != ".pdb":
            convert_molecule(lig, lig_name + ".pdb")
            lig = lig_name + ".pdb"
        cmd = F_PREPARE_COMPLEX.format(rec, lig, out_prefix+"_interactions_complex.pdb")
        execute_cmd(cmd)
        # prepare_complexes_pdb(rec, lig, out_prefix+"_interactions_complex.pdb")

    run_plip(out_prefix)

    read_file_plip_txt(out_prefix+"_interactions")
    cmd = F_EDIT_SESSION.format(out_prefix+"_interactions.pse")
    
    execute_cmd(cmd)
    #edit_sessions(out_prefix+"_interactions.pse")
