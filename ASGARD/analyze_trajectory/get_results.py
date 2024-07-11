#!/usr/bin/env python
# -*- coding: utf-8 -*-
from GetResults import *
import sys
import os
os.environ['GMX_MAXBACKUP'] = "-1"
os.environ['GMXLIB'] = os.getcwd()+"/ASGARD/external_sw/gromacs/force_field/" #required for g_mmpbsa

if len(sys.argv) < 4:
    print("\nDebe introducir:")
    print("Absolute path of the folder with the calculations + the suffix")
    print("Profile: [ TARGET_QUERY | TARGET | TARGET_QUERIES | QUERIES | DNA_QUERY ]  ")
    print ("gromacs_run [ gmx | gmx_mpi ]" )
    print("Example:")
    print("python /home/alejandro/ASGARD/VS_GR_mpro_ritonavir_results_2023_01_02/molecules/VS_GR_mpro_ritonavir TARGET_QUERY gmx_mpi\n")
    exit()
try:
    out_molec = sys.argv[1]
    profile =   sys.argv[2].upper()
    gromacs =   sys.argv[3]
    ligand = None
    reference = None

    for arg in sys.argv[4:]:
        if arg.startswith("ligand="):
            ligand = arg.split("=")[1]
        elif arg.startswith("reference="):
            reference = arg.split("=")[1]
except ValueError:
    print ("Error con la entra de datos")
    exit()


#
# Create configuration object
#
cfg = ConfigHolder(out_molec, profile, gromacs, ligand, reference)
Resume(cfg)

#
# parallel JObs
#
GraphStepFluctuation(cfg)
GraphsInteractionsTargetQueries(cfg)
GraphDssp(cfg)
GraphSasa(cfg)
#
# End of parallel jobs
#
GraphStabilization(cfg)
GraphRmsd(cfg)
GraphDistance(cfg)
#GraphHelix(cfg)
GraphGyrate(cfg)
TableMultiMolecule(cfg)
#
# Generate latex
#
GenerateLatex( cfg)

