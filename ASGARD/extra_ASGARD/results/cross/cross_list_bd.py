# !/usr/bin/env python
# -*- coding: utf-8 -*-
#   NO SE SUSA
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#   Description: Examina varios Bds proteina ligando, 
#           ojo siempre el mismo target y query y genera una tabla constrastando resultados.
#    
import sys
import json
import numpy
import math
import os
import re                               #usado para ordenar
from collections import OrderedDict

key_pat = re.compile(r"^(\D+)(\d+)$")   #usado para ordenar el diccionario
MAX_DISTANCE= 5
FAIL = '\033[91m'       #colores para error
WARNING = '\033[93m'
ENDC = '\033[0m'
FILE_CLUSTERS="{0}/{1}_clusters.json"
"""
    Muestra los errores
"""
def print_error( lst_text ):    
    print ("")
    if isinstance(lst_text, list) :
        for text in lst_text:
            print(("{}Error:{} {} {}".format(FAIL,WARNING,text,ENDC)))
    else:
        print(("{}Error:{} {} {}".format(FAIL,WARNING,lst_text,ENDC)))
    print ("")
    exit()


def key(item):
    m = key_pat.match(item[0])
    return m.group(1), int(m.group(2))

def read_json(file ):
    with open(file) as json_file:  
        data = json.load(json_file)
    
    data = sorted(data.items(),  key=key)
    
    for i in   data  : 
        print (i)

 

    return data
        
def distance_clusters(c_1, c_2):
    return math.sqrt((c_1[0] - c_2[0])**2 + (c_1[1] - c_2[1])**2+ (c_1[2] - c_2[2])**2)  
    

    

"""
    Main
"""
format_text_10 = '{:>10}'
format_text = '{:>20}'
format_text_coord='{:>25}'
if (len(sys.argv) < 2 ):
    print_error(['Debe introudcir las carpetas BDs','BD_[A-Z][A-Z]_1le0_GLA*']) 
array_bd_cl_folders = []
array_bd_cl = OrderedDict()
lst_sw = [] 
for i in range(1, len(sys.argv)):
    if not "tar.gz" in sys.argv[i] and os.path.exists(sys.argv[i]+"/energies" ):        
        array_bd_cl_folders.append(sys.argv[i]) #almacena todo los nombres de BD directorios que se le pasan
#
#   Leemos los ficheros de clusters
#

header = format_text_coord.format('Cluster')
for i in array_bd_cl_folders:
    sw = os.path.basename( i).split("_")[1]     
    header += format_text.format(sw)
    lst_sw.append(sw)
    print (sw)
    if os.path.isfile(FILE_CLUSTERS.format(i,os.path.basename(i) ) ):

        data = read_json (FILE_CLUSTERS.format(i,os.path.basename(i) ) )

        for j in data:   #SE AÑADE LA CARPETA Y EL SOFTWARE
            j[1]['sw']=sw        
            j[1]['folder'] = i
                    
        array_bd_cl[sw] = data
    else:
        print ("File josn no exists "+FILE_CLUSTERS.format(i,os.path.basename(i) ))

print  ("")
print (header)
ofset=""
lst_fill_sw = []

lst_cluster_json = {}

for sw, lst_cluter in list(array_bd_cl.items()):
    
    ofset  = [format_text.format('--') for num in lst_fill_sw] 
    ofset = ofset = "".join(ofset) if len(ofset)>0 else ""    
    lst_fill_sw.append(sw)     
    for nl in range (0 , len (lst_cluter)):
        lst_cluter_keep = []
        cl = lst_cluter[nl]
        lst_cluter_keep.append(cl[1])
        
        row = format_text_coord.format(str(cl[1]['coords'])[1:len(str(cl[1]['coords']))-1]) + str(ofset)           
        row += format_text.format(cl[0])
        cnt_cluster = 1
        row_scores = format_text_coord.format('')+ str(ofset)   
        score = (round(cl[1]['score'],2))
        row_scores += format_text.format(str(score))
        for sw_2, lst_cluter_2 in list(array_bd_cl.items()):  
            
            if sw_2 not in lst_fill_sw:
                find = False
                for cl_2 in lst_cluter_2:                                                          
                    distance = distance_clusters( cl[1]['coords'], cl_2[1]['coords'])
                    if  distance < MAX_DISTANCE:
                        score = (round(cl_2[1]['score'],2));
                        row += format_text.format(cl_2[0] )
                        row_scores += format_text_10.format(str( round(distance, 1)))
                        row_scores += format_text_10.format(str(score))
                        cnt_cluster += 1
                        lst_cluter_2.remove(cl_2)
                        lst_cluter_keep.append(cl_2[1])
                        find = True
                        break
                if find == False:
                    row += format_text.format('--')    
                    row_scores += format_text_10.format("")    
                    row_scores += format_text_10.format('--')    
        if cnt_cluster > 1:            
            lst_cluster_json[cl[0]] = lst_cluter_keep
            print (row)
            print (row_scores)
            print ("")

f_json = (os.path.basename( array_bd_cl_folders[0]).split("_")[2])+".json"
with open(f_json, 'w') as json_file:
    json.dump(lst_cluster_json, json_file,sort_keys=True)