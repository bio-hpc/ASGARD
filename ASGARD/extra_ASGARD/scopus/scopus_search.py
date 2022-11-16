#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Antonio Banegas
#   Email:
#   Date:   7/3/2020
#   Description: Busca en el api de scopus las palabras indicadas y devueleve el número de resultados
#   API url: https://dev.elsevier.com/api_docs.html
# ______________________________________________________________________________________________________________________
import argparse
import json
import math
import os
import sys
import pandas as pd

# Flujo principal
try:

    from six.moves import urllib
except:
    import urllib.request

API_KEY = "220684153fe6faa1b6994202cfec8fe3"  # Hprez
FORMAT_SEARCH = "TITLE-ABS-KEY({})"
URL_API_SCOPUS = 'https://api.elsevier.com/content/search/scopus?query={}&apiKey={}'
# Filtrar por los campos deseados acelera la busqueda
URL_API_SCOPUS_ITER = 'https://api.elsevier.com/content/search/scopus?query={}&apiKey={}&start={}&sort=citedby-count&field=title,publicationName,doi,citedby-count'


# Comprueba si la propiedad existe. Si no existe, devuelve vacio. Si existe, devuelve su valor
def get_property(entry, key):
    if key not in entry or entry[key] is None:
        return ""
    return entry[key].encode("utf-8")


parser = argparse.ArgumentParser()
parser.add_argument('search_words', help="Keyterms to search for")
parser.add_argument('output_file', help="Output file")
args = parser.parse_args()
search = args.search_words
output_file = args.output_file

print("Query: " + search)

# Eliminar el fichero de salida si ya existe
try:
    os.remove(output_file)
except Exception as e:
    print("Error while deleting file ", output_file)
    print(e)

str_search = ""
if len(search.split('+')) > 1:
    for i in search.split('+'):
        str_search += '{}+'.format(urllib.parse.quote(i))
    str_search = str_search[:-1]
else:
    str_search = urllib.parse.quote(search)
url = URL_API_SCOPUS.format(FORMAT_SEARCH.format(str_search), API_KEY)

print("Formatted Query: " + str_search)

try:
    response = urllib.request.urlopen(url)
    data = response.read()  # a `bytes` object
    text = data.decode('utf-8')  # a `str`; this
    json_data = json.loads(text)
    totalResults = int(json_data['search-results']['opensearch:totalResults'])
    itemsPerPage = int(json_data['search-results']['opensearch:itemsPerPage'])

    print("Total Results: " + str(totalResults))
    print("Items Per Page: " + str(itemsPerPage))

    # Iterar con cada indice
    start = 0
    while start < totalResults:

        print("From index " + str(start) + " out of " + str(totalResults) +
              " items")

        url = URL_API_SCOPUS_ITER.format(
            FORMAT_SEARCH.format(str_search), API_KEY, str(start))
        response_i = urllib.request.urlopen(url)
        data = response_i.read()  # a `bytes` object
        text = data.decode('utf-8')  # a `str`; this
        json_data = json.loads(text)
        output = []

        # Procesar la lista de entradas
        for entry in json_data['search-results']['entry']:
            # Parsear y guardar la salida
            info = {
                'Citations': get_property(entry, 'citedby-count'),
                'Title': get_property(entry, 'dc:title'),
                #'Year': '',
                #'Author(s)': '',
                'Journal': get_property(entry, 'prism:publicationName'),
                #'Volume': get_property(entry,'prism:volume'),
                #'Issue': get_property(entry,'prism:issueIdentifier'),
                #'Pages': get_property(entry,'prism:pageRange'),
                'DOI': get_property(entry, 'prism:doi')
            }

            output.append(info)

        # Se guarda tras cada iteracion para evitar problemas de memoria si hay muchos resultados.
        df = pd.DataFrame(output)
        df.to_csv(
            output_file,
            sep='\t',
            na_rep="none",
            index=False,
            mode="a",
            header=False)

        start += len(output)

        print("Added " + str(len(output)) + " entries\n")
        

    print("Search Complete!")
    print(str(totalResults) + " Results Written Into " + output_file )

except Exception as e:
    print(e)
