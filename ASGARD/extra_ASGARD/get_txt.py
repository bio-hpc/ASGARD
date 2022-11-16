#
#   Muestra por pantalla los txt antiguos de lanzador
# ________________________________________________________________________
import sys
import glob
import json
from collections import OrderedDict

if len(sys.argv) != 2:
    print("")
    print("Debe indicar una carpeta generada por Shuttlemol ejem:")
    print("print  ShuttleMol/extra_shuttlemol/get_txt.py BD_AD_1le0_GLA_2018-07-25/")
    print ("")
    exit()

data={}
cnt_json=0
directory = sys.argv[1]
for f_json in glob.glob(directory+"/energies/*.json"):
    with open(f_json) as f:
        d = json.load(f)
    data [cnt_json] = [ float( d['global_score'] ),  float(d['coords'][0]),float(d['coords'][1]),float(d['coords'][2]),str(d['file_orig']), int(d['num_execution']) ]
    cnt_json += 1
data=OrderedDict(sorted(data.items(), key=lambda t: t[1][0]))

for i, k in data.items():
    print (str(k)[1:-1])
