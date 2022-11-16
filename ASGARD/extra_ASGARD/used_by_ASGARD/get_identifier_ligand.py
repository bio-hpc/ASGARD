#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os, re
ALLOW_EXT = [ ".pdbqt", ".mol2" ]

def get_id_line(lst_lines, token):
    for c_line in  range(0, len(lst_lines)):
        if  token in lst_lines[c_line]:
            return c_line
    return -1


if len(sys.argv) != 2:
    print ("Error Parameters")
    print ("1º ligand")
    exit()
file = sys.argv[1]
file_name, file_ext = os.path.splitext(file)

if not file_ext in ALLOW_EXT:
    print ("ERROR extenxion not allowed")
    exit()
if not os.path.isfile(file):
    print("ERROR file not exits")
    exit()
lst_lines = []
with open(file, "r") as fileHandler:
    for line in fileHandler:
        lst_lines.append( line.strip() )

if file_ext == ".mol2":
    token = "@<TRIPOS>MOLECULE" # despues de esta linea esta el identificador de la molecula
    c_line = get_id_line(lst_lines, token)
    if c_line != -1 and lst_lines[c_line+1].strip() != "":
        print (lst_lines[c_line+1].strip())
        exit()

elif file_ext == ".pdbqt":
    token="COMPND"
    c_line = get_id_line(lst_lines, token)
    if c_line != -1 and lst_lines[c_line].replace(token,"").strip() != "":
        print (lst_lines[c_line].replace(token,"").strip())
        exit()

file_name = os.path.basename(file_name)
print( re.sub("[B|V][D|S]_[A-Z][A-Z]_","",file_name) ) # Eliminamos los prefijos de shuttlemol
