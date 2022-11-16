#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#   Description: Junta diferentes VSs
# ______________________________________________________________________________________________________________________
import argparse
import glob
import json
import os
from os.path import dirname, join
import os.path
from collections import OrderedDict
from subprocess import Popen, PIPE, STDOUT
import datetime
from shutil import copyfile
import tarfile

FORMAT_OUT='lst_{}{}.txt'
PYTHON_RUN = "python " 
JOIN_CL_SESSIONS = "ShuttleMol/extra_shuttlemol/results/join/join_cl_json_vs_session.py"
F_JOIN_SESSIONS = PYTHON_RUN + ' '+JOIN_CL_SESSIONS+ ' -f {} -d {} -r {} -o {} -v' 


parser = argparse.ArgumentParser()
parser.add_argument('d_energies', type=str, nargs='+', help='Directorio de la pruebas (minimo 2) ')
parser.add_argument('-c', '--cutoff', default=500, type=int)
parser.add_argument('-r', '--receptor', default='', type=str)

args = parser.parse_args()
if len(args.d_energies) == 1:
    parser.print_help()
    exit()



def read_energies(dir):
    dct = {}
    for file in glob.glob(dir+"*.json"):    	    
        with open(file) as json_file:
            data = json.load(json_file)
        key = os.path.splitext(os.path.basename(data['file_ori_query']))[0]
        try:
            dct[key] = float(data['global_score'])
        except:
            os.remove(file)
            print ("Erro: "+file)        
    lst = sorted(dct.items(), key=lambda x: x[1])
   
    dd = OrderedDict(sorted(dct.items(), key=lambda x: x[1])[:args.cutoff])
    return dd

def make_tarfile(source_dir, output_filename):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def execute_cmd(cmd):
		p = Popen(cmd, stdout = PIPE,  stderr = STDOUT, shell = True)
		for line in iter(p.stdout.readline, b''):
			print (line.strip())
		p.stdout.close()
		p.wait()
		

directories = '' # usada para llamar a join_cl
all = OrderedDict()
first_sw = ""
name_out= ''
for dir in args.d_energies:
	
    index = dir.find('VS_')
    directories += '{} '.format(dirname(dirname(dir)))
    sw = dir[index+3:index+5]
    name_out += '{}_'.format(sw) 
    all[sw] = read_energies(dir)
    if first_sw == "":
        first_sw = sw
name_out = name_out[:-1]


prefix_out = dirname(dirname(args.d_energies[0]))
prefix_out = prefix_out[prefix_out.find(first_sw)+len(first_sw):] 
file_out =  FORMAT_OUT.format(name_out, prefix_out)


header = " Rank ".join(all.keys())+" Rank Molecule"
header = ", ".join(all.keys())

f_out = open(file_out,'w')
print (header)
f_out.write('{}\n'.format(header))
rank = 1;
for molecule in all[first_sw]:
    if rank > args.cutoff:
        break
    score = all[first_sw][molecule]
    score_aux = ""
    for k, v in all.items():

        if k != first_sw:
            if v.has_key(molecule):
                score_aux += '{} {} '.format( v.keys().index(molecule), v[molecule])
            else:
                score_aux += "-- --"
    print (' {} {} {} {} '.format(rank, score,score_aux, molecule))
    f_out.write(' {} {} {} {} \n'.format(rank, score,score_aux, molecule))
    rank += 1
f_out.close()


if (args.receptor):
	if (os.path.isfile(args.receptor)):
		out_join = file_out[:file_out.rindex("_")]	
		cmd = F_JOIN_SESSIONS.format(file_out, directories, args.receptor, out_join)		
		f_cl = output = '{}_{}'.format(out_join, datetime.date.today())		
		print (cmd)
		execute_cmd(cmd)
		copyfile(file_out, join(f_cl, file_out))
		make_tarfile(f_cl, '{}.tar.gz'.format(f_cl))		
	else:
		print ("Error no existe el receptor")