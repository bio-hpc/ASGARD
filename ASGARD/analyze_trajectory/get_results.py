#!/usr/bin/env python
# -*- coding: utf-8 -*-
from GetResults import *
import sys
import os
os.environ['GMX_MAXBACKUP'] = "-1"
os.environ['GMXLIB'] = os.getcwd()+"/ASGARD/external_sw/gromacs/force_field/" #necesario para g_mmpbsa

if len(sys.argv) != 4:
    print("\nDebe introducir:")
    print("Ruta absoluta de la carpeta de la prueba mas el sufijo")
    print("Profile: [ TARGET_QUERY | TARGET | TARGET_QUERIES | QUERIES | DNA_QUERY ]  ")
    print ("gromacs_run [ gmx | gmx_mpi ]" )
    print("Ejemplo:")
    print("python /home/alejandro/ASGARD/VS_GR_mpro_ritonavir_results_2023_01_02/molecules/VS_GR_mpro_ritonavir TARGET_QUERY gmx_mpi\n")
    exit()
try:
    out_molec = sys.argv[1]
    profile = sys.argv[2].upper()
    gromacs = sys.argv[3]
except ValueError:
    print ("Error con la entra de datos")
    exit()
#
# Creamos obj de configuracion
#
cfg = ConfigHolder(out_molec, profile, gromacs)
Resume(cfg)

#
# JObs en paralalelo
#
GraphStepFluctuation(cfg)
GraphsInteractionsTargetQueries(cfg)
GraphDssp(cfg)
GraphSasa(cfg)
#
# Fin jobs paralelo
#
GraphStabilization(cfg)
GraphRmsd(cfg)
GraphDistance(cfg)
#GraphHelix(cfg)
GraphGyrate(cfg)
TableMultiMolecule(cfg)
#
# Generar latex
#
GenerateLatex( cfg)

