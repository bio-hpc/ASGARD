#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#    Comprueba la distancia entre [ ligando | receptor ] y un punto
#    Es util para gromacs
#    depende de:
#         "Shuttlemol/extra_shuttlemol/used_by_shuttlemol/standar_file_coords.py"
#
import sys
import subprocess
import math
STANDAR_FILE_COORDS = "ShuttleMol/extra_shuttlemol/used_by_shuttlemol/standar_file_coords.py"

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Error: Usage")
        print ("1º ligand")
        print("2º point 2,3,2 ")
        print("python ShuttleMol/extra_shuttlemol/distance_ligand_point.py ligand.mol2 5,5,5")
        exit()
    query = sys.argv[1]
    coords = sys.argv[2]
    x = float(coords.split(",")[0])
    y = float(coords.split(",")[1])
    z = float(coords.split(",")[2]) 
    cmd = 'python '+STANDAR_FILE_COORDS+ ' '+ query 
    out = subprocess.check_output(cmd, shell=True)
    x_sum = 0
    y_sum = 0
    z_sum = 0
    cnt_row = 0
    for i in out.split("\n"):
        if not i.startswith("#") and i.strip() != "":
            aux = i.split(":")
            x_sum += float(aux[0])
            y_sum += float(aux[1])
            z_sum += float(aux[2])
            cnt_row += 1
    x_center = x_sum / cnt_row
    y_center = y_sum / cnt_row
    z_center = z_sum / cnt_row
    distance = math.sqrt ( pow( ( x - x_center ), 2 ) + pow ( ( y - y_center ), 2) + pow( ( z - z_center ), 2) )

    print(distance)
