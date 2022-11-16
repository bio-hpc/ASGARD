# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   Author: Jorge de la Peña García
   Email:  jorge.dlpg@gmail.com
   Description: Bucle para convertir carpetas de moleculas

"""
import argparse
import os
from glob import glob
from os.path import join, splitext
import subprocess
import sys

PYTHON_RUN = 'python'
F_SCRIPT_CONVERT = PYTHON_RUN + ' ShuttleMol/extra_shuttlemol/used_by_shuttlemol/convert_to.py {} {}'
BABEL_RUN = 'babel '
F_BABEL = BABEL_RUN + ' {} -opdbqt {}'

"""
    Ejecuta un comadno
"""
ON_POSIX = 'posix' in sys.builtin_module_names
from threading  import Thread
try:
    from queue import Queue, Empty
except ImportError:
    from Queue import Queue, Empty  # python 2.x
def enqueue_output(out, queue):
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()
def execute_cmd(cmd):
    print (cmd)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1,
                            close_fds=ON_POSIX)
    # stdout = proc.stdout.read()
    q = Queue()
    stdout = Thread(target=enqueue_output, args=(proc.stdout, q))
    stdout.daemon = True  # thread dies with the program
    stdout.start()
    stderr = proc.stderr.read()

    if stderr:
        print ("ERROR: {}".format(stderr.strip()))
        print ('ERROR CMD:  {}'.format(cmd))
    else:
        print ('{}'.format(cmd))

    return stdout


class readable_dir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir = values.strip()
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest, prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

parser = argparse.ArgumentParser()
parser.add_argument('folder_in',  action=readable_dir, help="molecule folder")
parser.add_argument('ext_in', choices=['mol2', 'pdbqt', 'sdf', 'oeb'], help="Input extension [ mol2 | pdbqt | sdf | oeb ]")
parser.add_argument('ext_out', choices=['mol2', 'pdbqt', 'sdf', 'oeb'], help="Output extension [ mol2 | pdbqt | sdf | oeb ]")
args = parser.parse_args()


folder_pattern = join(args.folder_in, '*'+str(args.ext_in))
for molecule in glob(folder_pattern):
    mol = '{}.{}'.format(splitext(molecule)[0], args.ext_out)
    execute_cmd (F_SCRIPT_CONVERT.format(molecule, mol))

if args.folder_in .endswith('/'):
    args.folder_in = args.folder_in[:len(args.folder_in)-1]

folder_pattern = join(args.folder_in, '*'+str(args.ext_out))

print ("\n"+F_BABEL.format(folder_pattern,  args.folder_in+'_all.pdbqt'))





