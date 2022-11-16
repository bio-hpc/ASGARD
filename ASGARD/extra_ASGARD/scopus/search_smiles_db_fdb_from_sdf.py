#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#   Date:   7/2/2109
#   Description: Lee ficheros sdf y busca en deugbank o FDB
# ______________________________________________________________________________________________________________________
import argparse
import json
import urllib
from glob import glob
from bs4 import BeautifulSoup
import re
try:
    from six.moves import urllib
except:
    import urllib.request
FOLDER_MOLECULES = 'molecules'
TAG_END_MOL = '$$$$'
TAG_SCORE = '> <Score>'
TAGS_NAME = ['DATABASE_ID', 'VillaPharma_Code'
             ]  #special case for the name of the molecules VP and drugbank
F_DB_URL = "https://www.drugbank.ca/drugs/{}"
F_FDB_URL = "https://foodb.ca/compounds//{}"


def clean_firs_line_molecule(molecule):
    """
        The sdf files cannot start a blank line, this function removes that line if it exists.
        :return:
    """
    if molecule[0] == '':
        molecule.pop(0)
    return molecule


def get_score_molecule(molecule):
    """
        Return score fo molecule
        :param molecule:
    """
    for cnt_line in range(len(molecule)):
        if TAG_SCORE in molecule[cnt_line]:
            return float(molecule[cnt_line + 1])
    return 0


def get_name_molecule(molecule):
    """
        return the name of the molecule
        :param str molecule:
    """
    for cnt_line in range(0, len(molecule)):
        for tag_name in TAGS_NAME:
            if tag_name in molecule[cnt_line]:
                return (molecule[cnt_line + 1])
    if molecule[0] != '':
        return molecule[0]
    else:
        print("Error No name for molecule")
        exit(0)


def read_molecules(dict_molecules, f_sdf):
    """
        Reads an sdf file, separates the moleulas by the delimiter and if score is better than cut-off stores it in a list
        :param f_sdf:
        :return:
    """

    f = open(f_sdf)
    molecule = []
    for line in f:
        molecule.append(line[:-1])
        if TAG_END_MOL in line:
            molecule = clean_firs_line_molecule(molecule)
            name_molecule = get_name_molecule(molecule)
            #score = get_score_molecule(molecule)
            dict_molecules[name_molecule] = molecule

            molecule = []
    f.close()
    return dict_molecules, f_sdf


def get_number_article(str_search):
    API_KEY = "220684153fe6faa1b6994202cfec8fe3"  # Hprez
    FORMAT_SEARCH = "TITLE-ABS-KEY({})"
    URL_API_SCOPUS = 'https://api.elsevier.com/content/search/scopus?query={}&apiKey={}'
    """
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
    """

    url = URL_API_SCOPUS.format(FORMAT_SEARCH.format(str_search), API_KEY)

    try:
        response = urllib.request.urlopen(url)
        data = response.read()  # a `bytes` object
        text = data.decode('utf-8')  # a `str`; this
        json_data = json.loads(text)
        return (json_data['search-results']['opensearch:totalResults'])
    except Exception as e:
        return ('0')


def drugbank(name, molecule):
    try:
        lst_fields = ['Name', 'Groups', 'CAS number', 'SMILES']
        url = F_DB_URL.format(name)
        response = urllib.request.urlopen(url)
        data = response.read()  # a `bytes` object
        text = data.decode('utf-8')  # a `str`; this
        soup = BeautifulSoup(text, "html.parser")
        regex = re.compile("col-sm-")
        content_lis = soup.find_all(attrs={'class': regex})
        dict = {}
        for i in range(len(content_lis)):
            if (content_lis[i].text in lst_fields):
                #if (content_lis[i].text == "Name"):
                dict[content_lis[i].text] = content_lis[i + 1].text.replace(
                    ",", "_")
                #else:
                #    dict[content_lis[i].text]=content_lis[i+1].text
        return url, dict["Name"], get_number_article(
            dict["Name"] + "+coronavirus"), dict['Groups'], dict[
                'CAS number'], dict['SMILES']
    except:
        return "--", "--", "--", "--", "--", "--"


def fdb(name, value):
    lst_fields = ['FooDB Name', 'Groups', 'CAS Number', 'Isomeric SMILES']
    url = F_FDB_URL.format(name.strip())
    try:
        response = urllib.request.urlopen(url)
        data = response.read()  # a `bytes` object
        text = data.decode('utf-8')  # a `str`; this
        soup = BeautifulSoup(text, "html.parser")
        #regex = re.compile("-")
        content_lis = soup.find_all("tr")
        dict = {}
        for i in range(len(content_lis)):
            try:
                if (content_lis[i].find("th").text in lst_fields):

                    #if (content_lis[i].find("th").text == "FooDB Name"):
                    dict[content_lis[i].find("th").text] = content_lis[i].find(
                        "td").text.replace(",", "_")
                    #else:
                    #    dict[content_lis[i].find("th").text] =  content_lis[i].find("td").text

            except:
                pass
        return url, dict["FooDB Name"], get_number_article(
            dict["FooDB Name"] + "+coronavirus"), "--", dict[
                'CAS Number'], dict['Isomeric SMILES']
    except:
        return "--", "--", "--", "--", "--", "--"


dict_molecules = {}
path = "ligands/*.sdf"
#path="*.sdf"
for lig in glob(path):
    read_molecules(dict_molecules, lig)
f_out = open("compounds_ls_coronavirus.csv", 'w')
f_out.write('{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(
    'Score', 'Code', 'Url', 'Name', 'Match', 'Group', 'CAS Number', 'Smiles'))

for key, value in dict_molecules.items():
    if key.startswith("DB"):
        line = str(get_score_molecule(value)) + "," + key + "," + ','.join(
            drugbank(key, value))
    elif key.startswith("FDB"):
        score = str(get_score_molecule(value))
        line = '{},{},{}'.format(score, key.strip(), ','.join(fdb(key, value)))
    else:
        print("Error " + key)
        #exit();
    print(line)
    f_out.write('{}\n'.format(line))

f_out.close()
