#!/usr/bin/env python
# -*- coding: utf-8 -*-
import commands
import sys
import shutil #copiar ficheros
from os import listdir
import os
import re
import numpy as np
def GenerarBoxIons():
	execute "echo 0| $gromacs5 editconf${mpi} -f ${outMolec}_complex.gro -o ${outMolec}_box.gro -c -princ -d ${paddingGrid} -bt ${typeGrid}"    # Generacion grid 
    execute "$gromacs5 genbox${mpi} -cp ${outMolec}_box.gro -cs spc216.gro  -o ${outMolec}_solv.gro -p ${outMolec}_topol.top "   				#rellena la grid  de disolvente
    execute "mv ${outMolec}_complex.gro ${folderOutUcm}"
    execute "mv ${outMolec}_box.gro ${folderOutUcm}"
    execute "mv ${folderMolec}\#${salida}_topol.top.1\#* ${folderOutUcm}${salida}_topol.top.Cube"

#
#   Genera Topologias
# 
if len(sys.argv)!=3:
	print "El script necesita:"
	print "1º fichero proteina"
	print "2º directorio ligando"
	exit()

#
#	Main

#
#	Recogemos variables
#
protein=sys.argv[1]
folderLigs=sys.argv[2]
filenameProt, fileExtensionProt = os.path.splitext(protein) 
outComplex,_=os.path.splitext(filenameProt)
outComplex=os.path.basename(outComplex)
#
#Ficheros topologials
#
gro=dirAtucal+filenameProt+".gro"
top=dirAtucal+filenameProt+".top"
itp=dirAtucal+filenameProt+"_porse.itp"

if result == "False":
	print ""
	print "ERROR: con la proteina "+protein
	print "LA cadena de aminoacidos esta rota"
	print ""
	exit()
else:
	charge=crearTopologiaProt(protein)
	for fichero in listdir(folderLigs):
		filenameLig, fileExtensionLig = os.path.splitext(fichero)
		if (fileExtensionLig==".mol2"):
			checlNameResLig(folderLigs+fichero)
			atoms, resu, listaGro=crearLigando(filenameLig, fileExtensionLig, folderLigs+filenameLig+"/", folderLigs)
			addLigandTopol(charge,atoms,resu,listaGro, gro, top,folderLigs+outComplex+"_"+filenameLig,filenameLig,folderLigs)
			getChargeItpLigand(folderLigs+filenameLig,folderLigs+outComplex+"_"+filenameLig, charge )
print charge
#
#	Se borran ficheros inutiles que se crearon
#
os.remove(gro)		
os.remove(itp)		
os.remove(top)
