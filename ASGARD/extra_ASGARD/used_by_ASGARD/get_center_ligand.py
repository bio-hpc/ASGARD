#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author: Jorge de la Peña García
# Email:  jorge.dlpg@gmail.com
# Description: Se le pasa por parametro un ligando devuelve el centro
#_________________________________________________________________________________________________________________________________________________________________________
import os,sys,re

ll=() #vector con las coordenadas x,y,z
valoresCoordenadas=[0,0,0]
nAtom=0;
#
# 	Dependiendo si es fichero pdb/pdbqt o mol las columanas son diferentes
#
def get_com():
	global nAtom;
	if extension == ".mol2":  #los ficheros mol 2 tienen  las cordenadas en 2 3 4
		if contadorDeModelo== 1 and linea.find("@<TRIPOS>ATOM")==-1:
			nAtom+=1;
			ll=linea.split( );
			valoresCoordenadas[0]+=float (ll[2]);
			valoresCoordenadas[1]+=float(ll[3]);
			valoresCoordenadas[2]+=float(ll[4]);
	else:
		if contadorDeModelo==1 and ('ATOM ' in linea or 'HETATM ' in linea ) : ##5, 6, 7 son la poscion de xyz del atomo	
			nAtom+=1;
			x=linea[30:38]
			y=linea[38:46]
			z=linea[46:54]
			valoresCoordenadas[0]+=float(x);
			valoresCoordenadas[1]+=float(y);
			valoresCoordenadas[2]+=float(z);
#
#	Main
#
if len(sys.argv) != 2:
	print ("escript get_ccenter_ligand.py")
	print ("Se especificar el fihero de query para calcular su centro")
	exit();
ligand = sys.argv[1]
filename, extension = os.path.splitext(ligand)

f=open(ligand);
contadorDeModelo=0;
for linea in f:
	#cadena exacta len(re.findall('\\bROOT\\b', linea))>0
	if  "ENDMDL" in linea or "MODEL" in linea or len(re.findall('\\bROOT\\b', linea))>0 or "TORSDOF" in linea or "@<TRIPOS>ATOM" in linea or "@<TRIPOS>BOND" in linea:
		contadorDeModelo+=1;
	linea=linea[:-1]
	if linea!= "":
		get_com()
f.close()

valoresCoordenadas[0] = round (valoresCoordenadas[0] / nAtom ,3);
valoresCoordenadas[1] = round (valoresCoordenadas[1] / nAtom ,3);
valoresCoordenadas[2] = round (valoresCoordenadas[2] / nAtom ,3);
print (str(valoresCoordenadas[0])+":"+str(valoresCoordenadas[1])+":"+str(valoresCoordenadas[2]))


