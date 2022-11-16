#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Genera Topologias para  gromacs
#
import getopt
from os import sys, path
sys.path.append(path.dirname(path.abspath(__file__)))
from options.Queries import Queries
from options.Target_query import Target_query
from options.Target_queries import Target_queries
from options.Target_one_query import Target_one_query
from options.Target import Target
from options.ConfigHolder import ConfigHolder
from options.components.tools import *

try:
    opts, args = getopt.getopt(sys.argv[1:], "h:t:q:p:g:")
    if len(opts) < 2:
        print("len menor 2")
        print_error_parameters('ERROR')
except getopt.GetoptError:
    print("except")
    print_error_parameters('ERROR')

profile = "TARGET_QUERY" #defecto
dir_queries = ""
target = ""
gmx = ""

for opt, arg in opts:
    if opt == '-h':
        print_error_parameters('')
    elif opt in ("-t"):
        target = arg
    elif opt in ("-q"):
        dir_queries = arg
    elif opt in ("-p"):
        profile = arg.upper()
    elif opt in ("-g"):
        gmx = arg.lower()


cfg = ConfigHolder(target, dir_queries, profile, gmx)
cmd = "{0} {1} {2}".format (cfg.python_run, cfg.check_protein, cfg.target)

if profile == "DNA_QUERY":
    profile = "TARGET_QUERY"
if profile ==  "BIPHSIC_SYSTEMS":
    profile = "QUERIES"
if profile not in cfg.lst_profiles:
    print(bcolors.FAIL + "\nERROR: profile not found"+ bcolors.ENDC +'\n')
    print(print_error_parameters(''))
elif profile != "QUERIES": #si es queries no tiene target y engañamos poniendolo a 1
    cfg.target_chains = execute_cmd(cmd)[0]
else:
    cfg.target_chains = [1]


if len(cfg.target_chains) > 2: # maximo 99 cadenas, creo que es imposible tener tantas
    print('\n')
    print(bcolors.FAIL + "ERROR: Target " + target + bcolors.ENDC)
    print(bcolors.FAIL + "ERROR: En la cadena de aminoacidos: "+bcolors.ENDC)
    print(bcolors.FAIL + "ERROR: " + cfg.target_chains + bcolors.ENDC +'\n')
    exit()
else:

    dynamic_class = profile[0].upper()+profile[1:].lower()
    klass = globals()[dynamic_class]
    instance = klass(cfg)
    instance.execute()




