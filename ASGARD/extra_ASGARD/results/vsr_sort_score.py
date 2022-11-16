#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Se le pasan varias carpetas de VSR generados con el mismo conjunto de proteinas y genera uin excell
#   Author: Jorge de la Peña García
#
import sys
import glob
import os
import json
FOLDER_ENERGY = "/energies/"
EXT_ENERGIES = ".json"

if len(sys.argv) < 3:
    print("Error: parameters")
    print("1º file coords.txt of VSR")
    print("2º VSR folders f.e [ VSR_AD_1 VSR_AD_2 | VSR_AD_* ]")
    exit()
file_coords = sys.argv[1]
lst_folder_vsr = []
lst_aux = {}    # ses un hash apra que se guarde el scroe y este ordenado
lst_n_receptors = {} # mapea los cosdigos de las proteinas con su nombre
lst_lines = {}
for i in range(2,  len(sys.argv)):
    lst_folder_vsr.append(sys.argv[i])


#
#   Lista de proteinas
#
f = open(file_coords)
for i in f:
    aux = i.strip().split(":")
    name_rec = os.path.splitext(aux[0])[0]
    lst_n_receptors[name_rec] = aux[len(aux)-1]
    lst_aux [name_rec] = "-"

f.close()

#
#  cabecera xml
#
line = ""
line_aux = ""
for prot , score in lst_aux.items():
    line_aux +=" ;"+lst_n_receptors[prot]
    line +=" ; "+ prot
print line_aux
print line
#
#   Iteramos por todas las carpetas
#
for folder in lst_folder_vsr:
    file_rec = ""
    name_lig = ""
    name_lig = folder.split("_")[2]
    for prot, score in lst_aux.items():
        lst_aux[prot]="--"
    for f_energy in glob.glob(folder+FOLDER_ENERGY+'*'+EXT_ENERGIES):
        filename = os.path.basename(f_energy)
        filename = filename[ filename.find("_")+1: ]
        filename = filename[ filename.find("_")+1: ]
        file_rec = filename[:filename.find(name_lig)-1]
        with open(f_energy) as f:
            data = json.load(f)
        score = data['global_score']
        lst_aux[file_rec] = score

    line = ""
    for prot , score in lst_aux.items():
        line +=" ; "+ score
    print name_lig, line 
    
    #lst_lines[name_lig] = line

#for n_lig, scores in lst_lines.items():   
#    print n_lig +" ; " +scores
    

