#
#   Crea los txt antiguso de lanzadsor
# ________________________________________________________________________
import sys
import glob
import json
import os
import ntpath

if len(sys.argv) != 2:
    print("")
    print("Debe indicar una carpeta generada por Shuttlemol ejem:")
    print("print  ShuttleMol/extra_shuttlemol/get_txt.py BD_AD_1le0_GLA_2018-07-25/")
    print ("")
    exit()

data= []
directory = sys.argv[1]
path_txt = directory+"/txt/"
if not os.path.exists(path_txt):
    os.makedirs(path_txt)

for f_json in glob.glob(directory+"/energies/*.json"):
    with open(f_json) as f:
        d = json.load(f)
    if d['global_score'].strip() != "":
	    f_out = os.path.basename(f_json)
	    f_out,_  = os.path.splitext(f_out)
	    s_out = float( d['global_score'] ),  float(d['coords'][0]),float(d['coords'][1]),float(d['coords'][2]), str(d['file_ori_query']), int(d['num_execution'])
	    file = open(os.path.join(path_txt,f_out+".txt"), "w")
	    file.write(str(s_out).replace(",","").replace("'","")[1:-1]+"\n")
	    file.close()




