#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Antonio Jesus Banegas Luna
#   Email:  ajbanegas@gmail.com
#	Date:	13/03/2020
#   Description: Busca en el api de ChEMBL las palabras indicadas y descarga los compuestos en formato mol.
#	API url: https://www.ebi.ac.uk/chembl/api/
# ______________________________________________________________________________________________________________________
import argparse
import json
import math
import os
import urllib
import urllib2

try:
    from six.moves import urllib
except:
    import urllib.request

# Escribe el contenido de un fichero de texto
def save_file(path, content):
	f = open(path, 'w')
	f.write(content)
	f.close()

# Flujo principal
FILE_PATH='{}/{}.mol'
URL_API_CHEMBL='https://www.ebi.ac.uk'
URL_API_CHEMBL_SEARCH =URL_API_CHEMBL+'/chembl/api/data/activity/search.json?q={}'
URL_API_CHEMBL_COMP =URL_API_CHEMBL+'/chembl/api/data/molecule/{}?format=json'

parser = argparse.ArgumentParser()
parser.add_argument('search_words', help="Keyterms to search for")
parser.add_argument('output_dir', help="Output directory")
args = parser.parse_args()
search = args.search_words
output_dir = args.output_dir

print "Query: "+search

# NOTA: No se elimina el directorio de salida por si contiene otros ficheros de utilidad

# Comprobar si el directorio de salida existe. Si no, crearlo
if not os.path.isdir(output_dir):
	os.mkdir(output_dir)
	print "Created Output Directory: "+output_dir

# Formatear la cadena de busqueda
str_search = ""
if len(search.split('+')) > 1:
    for i in search.split('+'):
        str_search += '{}+'.format(urllib.parse.quote(i))
    str_search = str_search[:-1]
else:
    str_search = urllib.parse.quote(search)
url = URL_API_CHEMBL_SEARCH.format(str_search, 0)
 
print "Formatted Query: "+str_search

try:
    response = urllib.request.urlopen(url)
    data = response.read()      # a `bytes` object
    text = data.decode('utf-8') # a `str`; this
    json_data = json.loads(text)

    totalResults = int(json_data['page_meta']['total_count'])
    itemsPerPage = int(json_data['page_meta']['limit'])
    nextPage = json_data['page_meta']['next']

    print "Total Results: "+str(totalResults)
    print "Items Per Page: "+str(itemsPerPage)

    # Iterar mientras haya una pagina siguiente
    counter = 0
    compounds = {}

    while nextPage != None:
	url = URL_API_CHEMBL+nextPage
	response_i = urllib.request.urlopen(url)
	data = response_i.read()
	text = data.decode('utf-8')
	json_data = json.loads(text)
	nextPage = json_data['page_meta']['next']	

	# Para cada entrada buscar el identificador del compuesto
	for entry in json_data['activities']:
		chembl_id = entry['molecule_chembl_id']

		# Si ya lo hemos descargado, lo saltamos
		if chembl_id in compounds:
			continue

		# Descargar el mol del compuesto
		url2 = URL_API_CHEMBL_COMP.format(chembl_id)
		response_c = urllib.request.urlopen(url2)
		data2 = response_c.read()
		json_data2 = json.loads(data2.decode('utf-8'))
		mol = json_data2['molecule_structures']['molfile']
		path = FILE_PATH.format(output_dir, chembl_id)

		save_file(path, mol)
		counter += 1
		compounds[chembl_id] = 1

	print str(len(json_data['activities']))+" Results Saved (Total: "+str(counter)+")"

    print "Search Complete!"
    print str(counter)+" Results Written Into "+output_dir

except Exception as e:
    print e
