#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#	Date:	13/03/2020
#   Description: Descarge un pdb
# ______________________________________________________________________________________________________________________
import argparse
import sys
import argparse
import urllib
try:
    from six.moves import urllib
except:
    import urllib.request

D_URL_PDB = "https://files.rcsb.org/download/"
F_URL_PDB = D_URL_PDB + "{}" 

parser = argparse.ArgumentParser(description='Download pdb.')
parser.add_argument('pdb_code', help='Pdb Code, example 3s5z')
parser.add_argument('-o', '--output', help='Pdb Code, example 3s5z')
args = parser.parse_args()
pdb = args.pdb_code+".pdb"
if args.output:
	out = args.output
else:
	out = pdb
try:
	url = F_URL_PDB.format(pdb)
	print  (url)
	response = urllib.request.urlopen(url)	
	data = response.read()

	f = open(out, 'w')
	f.write (data.decode('utf-8')) 
	f.close()
except Exception as e:
	print ("Error pdb not found")
	
