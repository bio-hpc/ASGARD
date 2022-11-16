# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#   Description: busca dentro de sdf generados por Ls los nombres y scores de las molecuals
#          
try:
	from urllib.request import urlopen
	import urllib as urllib
except ImportError:
	from urllib2 import urlopen
	import urllib2 as urllib
import argparse
from operator import itemgetter
import json
from bs4 import BeautifulSoup
import subprocess

SCRIPT_SCOPUS = 'ShuttleMol/extra_shuttlemol/scopus/scopus_search.sh'
CUTOFF = 0.7
TOKEN_SCORE = '> <Relative Pharmacophore-Fit Score>'
TOKEN_NAME = '> <DRUGBANK_ID>'
TOKENS = [TOKEN_SCORE, TOKEN_NAME]
# http get https://www.ebi.ac.uk/chembl/api/data/compound_record/07443 --json
API_DRUGBANK = 'https://www.drugbank.ca/drugs/'
FOLDER_OUT='test/'
KEY_WORD='Acetylcholinesterase' 
parser = argparse.ArgumentParser()
parser.add_argument('file_sdf', type=file)
args = parser.parse_args()

cnt = 0
aux = []
lst_all = []
with open(args.file_sdf.name) as file:
	read = False
	for line in file:
		if read:
			aux.append(line.strip())			
			if len(aux) > 1:
				if float(aux[1]) > CUTOFF:
					if aux not in lst_all:
						lst_all.append(aux)
					aux = []
		if line.strip() in TOKENS:
			read = True
		else:
			read = False
cnt_error= 0

for i in sorted(lst_all, key=itemgetter(1), reverse = True):	

		#if i[0] == 'DB07694':
	#	cnt_error+=1
	#if cnt_error>0:		
	request = urllib.Request('{}{}'.format(API_DRUGBANK, i[0]))
	response = urlopen(request)
	data = response.read()
	soup = BeautifulSoup(str(data), "html.parser")
	aux = soup.findAll("dd", {"class": "col-md-10 col-sm-8"})
	name = aux[1].text.split(" ")[0].strip()
	print ('{} {}'.format( i[0], name.encode('ascii', 'ignore')) )
	if name == i[0]:
		print ('{}  {}'.format(name, aux[0].text.encode('ascii', 'ignore')))
		out = FOLDER_OUT+i[0]
		cmd = '{} "{}+{}" "{}.txt"'.format(SCRIPT_SCOPUS, aux[0].text.encode('ascii', 'ignore'), KEY_WORD, out)
		output = subprocess.check_output(
			cmd,
			shell=True,
			stderr=subprocess.STDOUT,
			
		)
		#print (output)
			
		
