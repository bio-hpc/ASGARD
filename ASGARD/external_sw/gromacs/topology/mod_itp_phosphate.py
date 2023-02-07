#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#	Se encarga de leer una topologia del ligando y ponerlo los angulos diaedros inpropios en las formaciones 0 P 0 0
#  FOSFATO= PO4 -<fosforo + 4 oxigenos


import sys
import shutil #copiar ficheros
import re

#
#  Formatos
# 
formatInp='{0:>6}{1:>7}{2:>7}{3:>7}{4:>7}{5:>9}{6:>10}'
def readFile(topLig):
	lstItp = []
	f = open(topLig)
	for i in f:
		lstItp.append(i[:-1])
	f.close()
	return lstItp
def findAtoms(lstItp):
	lstAtom = []
	contador = 0
	tokken1 = "[ atoms ]"
	tokken2 = "[ bonds ]"
	fosfaos = []
	for i in lstItp:
		if i.find(tokken2) != -1:
			contador += 1
		if contador == 1 and not i.startswith(";") and not i.strip() == "":
			aux=re.sub(' +',' ',i).strip().split(" ")
			lstAtom.append(aux)
			if aux[4].startswith("p") or aux[4].startswith("P"):
				fosfaos.append(aux)
		if i.find(tokken1) != -1:
			contador += 1
	return lstAtom, fosfaos
def finParis(lstItp):
	lstBonds = []
	contador = 0
	tokken1 = "[ bonds ]"
	tokken2 = "[ pairs ]"
	for i in lstItp:
		if i.find(tokken2)!=-1:
			contador+=1
		if contador==1 and not i.startswith(";") and not i.strip()=="":
			aux=re.sub(' +',' ',i).strip().split(" ")
			lstBonds.append(aux)
		if i.find(tokken1)!=-1:
			contador+=1
	return lstBonds
def findAtomsNeighbor(lstBonds,fosfato):
	lstNeighbor=[]
	#lstNeighbor.append(fosfato[0])
	for i in lstBonds:
		if i[0]==fosfato[0]:
			lstNeighbor.append(i[1])
		elif i[1]==fosfato[0]:
			lstNeighbor.append(i[0])
	return lstNeighbor

def carbonNeighbo(lstNeighbor, lstBonds, lstAtom):
	for i in lstNeighbor:
		for j in lstBonds:
			if i==j[0]:
				for z in lstAtom:
					if z[0]==j[1] and (z[4].startswith("C") or z[4].startswith("C")):
						return i
			elif i==j[1]:
				for z in lstAtom:
					if z[0]==j[0] and (z[4].startswith("C") or z[4].startswith("C")):
						return i
def generarItp(lstItp,lstAddDiedrals,topoliLg):
	shutil.copyfile(topoliLg, topoliLg+".bk" ) #copia del original por si acaso
	contador = 0
	f = open(topoliLg, 'w')
	tokken1 = "[ dihedrals ]"
	for i in lstItp:
		if i.startswith(tokken1):
			contador = 1
		if contador >= 1:
			contador += 1
		if contador == 5:
			for j in lstAddDiedrals:
				f.write(j+"\n")
		f.write(i+"\n")
	f.close()

		


#
#Escreibe un fichero
#
def writeFile(lista, file):
	f= open(file, 'w')	
	for i in lista:
		f.write(str(i)+"\n")
	f.close()

if len(sys.argv)!=2:
	print ("debe inteorucir una topologia del ligando para agregar los angulos diaedros en los fosfamos necesarios")
	exit()

#
#	Recogemos variables
#
angulo = "-35.00"
muelle = "90.93200"
lstAddDiedrals = []
topoliLg = sys.argv[1]

lstItp = readFile(topoliLg)
lstAtom, fosfatos = findAtoms(lstItp)
lstBonds = finParis(lstItp)

for i in fosfatos:
	lstNeighbor=findAtomsNeighbor(lstBonds,i)
	atomJoinC=carbonNeighbo(lstNeighbor,lstBonds,lstAtom)
	lstNeighbor.remove(atomJoinC)
	lstNeighbor = sorted(lstNeighbor)
	if len(lstNeighbor) == 3: #si tiene 3 oxigenos
		lstAddDiedrals.append(formatInp.format(i[0],atomJoinC,lstNeighbor[1],lstNeighbor[2],2,angulo,muelle ))
		lstAddDiedrals.append(formatInp.format(i[0],atomJoinC,lstNeighbor[2],lstNeighbor[0],2,angulo,muelle ))
		lstAddDiedrals.append(formatInp.format(i[0],atomJoinC,lstNeighbor[0],lstNeighbor[1],2,angulo,muelle ))
if (len(lstAddDiedrals)>0):
	generarItp(lstItp,lstAddDiedrals,topoliLg)

	#print lstNeighbor

