#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   Author: Jorge de la Peña García
   Email:  jorge.dlpg@gmail.com
   Description: Recibe un fichero con las listas cruzadas de VS y las carpetas y generar una sesion pymol
   Madatory Modules:
            [ python3-lxml | python-lxml ]
"""
from collections import OrderedDict
import argparse
import os
from os.path import join, basename, splitext, isfile
import re
import sys
import fnmatch
from shutil import copyfile
import subprocess
import json
import datetime



PYTHON_RUN = 'python'
PLIP_SCRIPT = PYTHON_RUN + " " + join('ShuttleMol', 'extra_shuttlemol', 'used_by_shuttlemol', 'create_ligand_plip.py {} {} {}')
PYMOL_SCRIPT = PYTHON_RUN + " " + join('ShuttleMol', 'extra_shuttlemol', 'used_by_shuttlemol', 'create_ligand_pymol.py {} {} {} {} {}')
PML_HEAD_PYMOL = PYTHON_RUN + " ShuttleMol/extra_shuttlemol/used_by_shuttlemol/create_header_pml.py -t {}"

def read_all_files( header):
    all_files = {}
    for sw in header:

        path = get_string_pattern(args.dirs, '_{}_'.format(sw))
        if not path:    #si no se le ha psado la carpeta dara error
            print ("Error no folder exists for {}".format(sw))
            exit()              
        lst_files = []
        for root, dirs, files in os.walk(path):                
            for name in files:            
                lst_files.append(join(root, name))
                #if fnmatch.fnmatch(name, pattern):                
                #    result.append(join(root, name)) 
        all_files[sw] =  lst_files
    return all_files


def find(option_files, pattern):
    """
        find file
        https://stackoverflow.com/questions/1724693/find-a-file-in-python
        :param pattern:
        :param path:
        :return:
    """    
    result = []
    for i in option_files:
        if fnmatch.fnmatch(i, pattern):                
            result.append( i)
    return result


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def read_file(file_cross):
    """
        Receives a file with the cross_list_vvs format and returns a list with each line in an array
        :param str file_cross:
        :return: list with each line in an array
    """
    lst_cross = []
    with open(file_cross) as fp:    
        for line in fp:            
            if line.startswith(' '): # if you start with this mark it means that it is the header                
                line = re.sub(' +', ' ', line.strip())
                lst_cross.append(line.split(" "))    
    return lst_cross


def filter_date(lst_cross, num):
    """
         Filters the data with a certain method
        :param lst_cross: list with each line in an array
        :param int num: number of tests that can fail
        :return:
    """
    num = num * 2
    data = []
    for line in lst_cross:
        if (line.count('--')) == num:
            data.append(line)
    return data


def classification_strategy(lst_data, method):
    """
         Classification Strategy
        :param lst_data:
        :param header:
        :return: dictionary sorted by score
    """
    dict = {}
    for line in lst_data:
        

        sum = getattr(sys.modules[__name__], method)(line)        
        if sum not in dict:
            dict[sum] = line
        else: 
            while(sum not in dict):
                sum = str(sum)+1
            dict[sum] = line
            #print ("Error classification_strategy")
            #exit()
    dict = OrderedDict(sorter(dict)[:args.max_results])

    if args.verbose:
        for k, v in dict.items():
            print ('{:8.2f}'.format(k), v)
    return dict


def sorter(data):
    """
        Sort dictinary
        :param data: dict
        :return:
    """
    return sorted(data.items(), key=lambda kv: kv[0], reverse=False)


def by_rank(line):    
    """
        Adds the rankings of each program
        :param line:
        :return:
    """
    sum = 0
    cnt = 0
    for i in range(0, len(line[:-1]) ):
        if i % 2 == 0 and not '--' in line[i]:
            cnt += 1
            # print (line[i])
            sum += int(line[i])
    return sum / cnt


def by_score(line):    
    """
        adds the scores of each program
        :param line:
        :return: sum
    """
    sum = 0
    cnt = 0
    for i in range(0, len(line[:-1])):
        if i % 2 != 0 and not '--' in line[i]:
            cnt += 1
            sum += float(line[i])
    return sum / cnt


def make_folder(folder):
    """
        if it does not exist, create the folder
        :param str folder: name folder
        :return:
    """
    if not os.path.isdir(folder):
        os.makedirs(folder)


def params():
    """
        Collects the input parameters
        :return:
    """
    parser = argparse.ArgumentParser(
        description='Receive a file with VS cross lists and folders and generate a pymol session.',
        epilog='Copyright 2019 Autor bajo licencia GPL v3.0'
    )

    parser.add_argument('-f', '--file', type=argparse.FileType('r'), required=True, help='output cross_list_vs')
    parser.add_argument('-d', '--dirs', action='append', nargs='+',
                        help='Root folder generated by shuttlemol example: VS_AD_5hex_rec_DB-201016_PQ_31_-4_-65_2018-11-26',
                        required=True, type=is_dir)
    #parser.add_argument('-m', '--method', choices=['by_score', 'by_rank'], required=True, help='method of ordering')
    parser.add_argument('-r', '--receptor', required=True,  type=argparse.FileType('r'), help='receptor')
    parser.add_argument('-o', '--output', required=True, help='Output folder')
    parser.add_argument('-m', '--max_results', default=100, help='Max Results', type=int)
    parser.add_argument('-v', '--verbose', help='Verbose', action='store_true')
    
    
    return parser.parse_args()


def cp_files(files, folder):
    """
        Copies the files to the specified folder

        :param [] files:
        :param str folder:
        :return:
    """
    for file in files:
        copyfile(file, join(folder, basename(file) ))


def get_string_pattern(lst, pattern):
    """
        Returns the string of the list based on the pattern

        :param pattern:
        :return: 
    """
    for i in lst:
        if pattern in i:
            return i
    return False


def get_file_molecule(list):
    for ext in ['.mol2', '.pdbqt', '.pdb']: # know extensions for molecules
        query = get_string_pattern(files, ext)
        if query:
            return query
    print ("Error: get_file_molecule: Molecules not found ")
    exit()


def write_file(lst_pml, file):
    """
        "rites a list of lines to a plain text file

        :param lst_pml: list of lines
        :param file: output file
        :return:
    """
    f = open(file, "w")
    for i in lst_pml:
        f.write(i)
    f.close()


def add_receptor(receptor, folder_molecules):
    """
        Copy the protein to the folder of molecules and return the necessary lines pea insert it into the pml
        :param receptor:
        :param folder_molecules:
        :return:
    """

    cp_files([receptor], folder_molecules)
    receptor = "./../"+join(folder_molecules, basename(receptor))
    cmd = PML_HEAD_PYMOL.format(receptor)
    return subprocess.check_output(cmd, shell=True).decode('UTF-8')
    """
    name_protein = splitext(basename(receptor))[0]
    line = 'cmd.spectrum("b","rainbow_rev", \"tal\")\n' \
           + 'cmd.load(\'../' + receptor + '\',\'' + name_protein + '\')\n' \
           + 'cmd.set(\'transparency\', 0.6)\n' \
           + 'cmd.hide(\'everything\', \'' + name_protein + '\')\n' \
           + 'cmd.show_as(\'surface\', \'' + name_protein + '\')\n' \
           + 'cmd.color(\'green\', \'' + name_protein + '\')\n' \
           + 'cmd.hide("(all and hydro and (elem C extend 1))")\n'
    return line
    """


def read_json(file):
    """
        Reads a json file and returns its dictionary
        :param str file:
        :return:
    """
    if file:
        with open(file, 'r') as f:
            datastore = json.load(f)
    return datastore


def execute_command(cmd):
    """
        Execute commands

        :param cmd: command
        :return:
    """
    if args.verbose:
        print (cmd)
    return subprocess.check_output(cmd, shell=True).decode('UTF-8')


if __name__ == "__main__":
    print ("")
    args = params()
    args.dirs = args.dirs[0]
    args.output = '{}_{}'.format(args.output, datetime.date.today())
    with open(args.file.name) as f:
        header = f.readline().strip().replace(',','').split(" ")

    out_molecs = join(args.output,'Molecules', '')

    make_folder(args.output)
    make_folder(out_molecs)
    #
    #   Sort data
    #
    all_files = read_all_files(header)

    lst_cross = read_file(args.file.name)
    for cnt_discard in range(len(args.dirs)-1):   

        for method in ['by_score', 'by_rank']:
            #cnt_discard = 1
            
            lst_data = filter_date(lst_cross, cnt_discard)
            pml_file = join(args.output, '{}_{}.pml'.format(cnt_discard, method))
            
            dict_sort = classification_strategy(lst_data, method)
            
            #
            #   copy files
            #
            pml_lst = [add_receptor(args.receptor.name, out_molecs)]
            cnt_cluster = 1
          
            for k, v in dict_sort.items():
                print (k, v)
                group = []
                g_ranks = []
                for sw in range(len(header)):                                        
                    if v[sw*2] != '--':
                        g_ranks.append(header[sw] +"_"+ v[sw*2] )
                        if args.verbose:
                            print ('Search: *_{}_*'.format(v[-1]))
                        files = find(all_files[header[sw]] ,'*_{}_*'.format(v[-1]))
                        
                        cp_files(files, out_molecs)                    
                        query = get_file_molecule(files)
                        name_query = basename(query)
                        
                        query = join(out_molecs, name_query)
                        prefix = join(out_molecs, splitext(name_query)[0])
                        
                        out_json_plip = '{}_interactions.json'.format(splitext(query)[0])

                        if not get_string_pattern(files, 'interactions'):     
                            if not isfile(out_json_plip): #  if not exists interactions file                        
                                execute_command( PLIP_SCRIPT.format(args.receptor.name, query, prefix) )
                        
                        interactions_file = prefix + '_interactions.json'
                        energy_json = prefix + '.json'
                        cmd = PYMOL_SCRIPT.format('../'+query, energy_json, interactions_file, 0, False)
                        pml_lst.append( execute_command(cmd))
                        
                        score = round(float(read_json(energy_json)['global_score']), 2)

                        group.append ('{}_{}'.format(splitext(name_query)[0], score ))        
                        
                        
                
                pml_lst.append("cmd.group('CL_{} {} ( {} )', '{}')\n".format(cnt_cluster, round(k, 2), ' '.join(g_ranks), ' '.join(group)))
                cnt_cluster += 1
            
            write_file(pml_lst, pml_file)











