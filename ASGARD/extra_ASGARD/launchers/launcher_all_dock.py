#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Jorge de la Peña García
#	Email: 	jorge.dlpg@gmail.com
#	Description: Utiliza shuttlemol con diferentes progrmas de docking, antes de lanzar las pruebas comprueba qu esistan los ficheros necesarios
#
import sys
import os
import subprocess

JOBS_PER_OPT = 10 # Numero de jobs que se utilizaran para cada prueba 
SM_CALL = './sm.sh -t {} -q {} -o {} -s {} -j {} -hi y {}' #el ultimo cambpo sirve para las coordenads en VS
FAIL = '\033[91m'       #colores para error
WARNING = '\033[93m'
ENDC = '\033[0m'
""" 
    Difernetes opciones con los programas que tienen y sus extensiones
"""
OPTS = {
        'BD':{
            'AD':{ 'target_ext' :'pdbqt', 'query_ext': 'pdbqt'},
            'AK':{ 'target_ext' :'pdbqt', 'query_ext': 'pdbqt'},
            'LF':{ 'target_ext' :'mol2', 'query_ext': 'mol2'},
            'FB':{ 'target_ext' :'mol2', 'query_ext': 'mol2'},
            'FR':{ 'target_ext' :'pdb', 'query_ext': 'oeb'},
        },
        'VS':{
            'AD':{ 'target_ext' :'pdbqt', 'query_ext': 'pdbqt'},
            'AK':{ 'target_ext' :'pdbqt', 'query_ext': 'pdbqt'},
            'LF':{ 'target_ext' :'mol2', 'query_ext': 'mol2'},
            'FB':{ 'target_ext' :'mol2', 'query_ext': 'mol2'},
            'FR':{ 'target_ext' :'oeb', 'query_ext': 'oeb'},                      
        }
}

"""
 Comprueba si existen los ficheros para la prueba de dock
"""
def get_files(sw, t, q):
    new_t = '{}.{}'.format(os.path.splitext(t)[0],OPTS[opt][sw]['target_ext'])    
    new_q = '{}.{}'.format(os.path.splitext(q)[0],OPTS[opt][sw]['query_ext'])    
    if not os.path.isfile(new_t):
        print_error('El fihcero {} para {} no existe'.format(new_t, sw))                
    if not os.path.isfile(new_q):        
        print_error('El fihcero {} para {} no existe'.format(new_q, sw))        
    return new_t, new_q
"""
    Pinta los errores
"""
def print_error( lst_text ):    
    print ("")
    if isinstance(lst_text, list) :
        for text in lst_text:
            print ("{}Error:{} {} {}".format(FAIL,WARNING,text,ENDC))
    else:
        print ("{}Error:{} {} {}".format(FAIL,WARNING,lst_text,ENDC))
    print ("")
    exit()
#https://www.programcreek.com/python/example/2244/subprocess.CalledProcessError
def run_command(command, wait=False):

    try:
        if (wait):

            p = subprocess.Popen(
                [command], 
                stdout = subprocess.PIPE,
                shell = True)
            p.wait()
        else:
            p = subprocess.Popen(
                [command], 
                shell = True, 
                stdin = None, stdout = None, stderr = None, close_fds = True)

        (result, error) = p.communicate()
        
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            "common::run_command() : [ERROR]: output = %s, error code = %s\n" 
            % (e.output, e.returncode))

    return result 
"""
    Main
"""
if len(sys.argv) < 4:
    print_error(
        ['Debe introudcir un recptor, una proteina y opcion y en caso de VS tambien las coordenas', '','por ejemplo:',
         'python ShuttleMol/extra_shuttlemol/launcher_all_dock.py targets/test/1le0.pdbqt  queries/test/GLA.mol2 BD ', 
         'python ShuttleMol/extra_shuttlemol/launcher_all_dock.py targets/test/1le0.pdbqt  queries/test/GLA.mol2 VS -2 6 2  ' 
         ]
        )
    
    
target = sys.argv[1]
query = sys.argv[2]
opt = sys.argv[3].upper()
coords_VS = ""
if opt in OPTS:
    if opt == "VS":
        if  len(sys.argv) != 7:
            print_error ([ "Para VS debe indicar las coordendas x, y, z", '','por ejemplo','python ShuttleMol/extra_shuttlemol/launcher_all_dock.py targets/test/1le0.pdbqt  queries/test/GLA.mol2 VS -2 6 2  ' ])
        else:
            coords_VS = "-x {} -y {} -z {}".format(sys.argv[4], sys.argv[5], sys.argv[6])


    #
    #   Compruebo si existen todos los ficheros para utilizar los diferentes programas de docking 
    #
    for sw in OPTS[opt].keys():
        get_files(sw, target, query)            
    #
    #   ejecuto sm
    #
    for sw in OPTS[opt].keys():
        t, q  = get_files(sw, target, query)   
        if opt == "VS":           
            q = os.path.dirname(q)+"/"    
        cmd = SM_CALL.format(t, q, opt, sw, JOBS_PER_OPT, coords_VS)  
        run_command (cmd)
        #out = subprocess.check_output(cmd ,encoding='UTF-8', shell=True)
        #print (out)
            


else:
    print_error('Opcion desconocida')
    

