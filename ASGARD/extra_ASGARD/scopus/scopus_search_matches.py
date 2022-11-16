#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#	Date:	7/2/2109
#   Description: Busca en el api de scopus las palabras indicadas y devueleve el número de resultados
#	API url: https://dev.elsevier.com/api_docs.html
# ______________________________________________________________________________________________________________________
import argparse
import json
import urllib

try:
    from six.moves import urllib
except:
    import urllib.request

API_KEY = "220684153fe6faa1b6994202cfec8fe3" #	Hprez
FORMAT_SEARCH = "TITLE-ABS-KEY({})"
URL_API_SCOPUS ='https://api.elsevier.com/content/search/scopus?query={}&apiKey={}'
parser = argparse.ArgumentParser()   
parser.add_argument('search_words', help="Palabras a buscar")
args = parser.parse_args()
search = args.search_words
str_search = ""
if len(search.split('+')) > 1:
    for i in search.split('+'):
        str_search += '{}+'.format(urllib.parse.quote(i))
    str_search = str_search[:-1]
else:
    str_search = urllib.parse.quote(search)
url = URL_API_SCOPUS.format(FORMAT_SEARCH.format(str_search), API_KEY)
try:
    response = urllib.request.urlopen(url)
    data = response.read()      # a `bytes` object
    text = data.decode('utf-8') # a `str`; this
    json_data = json.loads(text)
    print (json_data['search-results']['opensearch:totalResults'])
except Exception as e:
    print ('0')
