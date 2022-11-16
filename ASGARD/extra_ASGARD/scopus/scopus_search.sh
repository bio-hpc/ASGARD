#!/bin/bash
#
#   Author: Jorge de la Peña García
#   Email:  jorge.dlpg@gmail.com
#	Date:	7/2/2109
#   Description: Busca en el api de scopus las palabras indicadas 
#	API url: https://dev.elsevier.com/api_docs.html
# ______________________________________________________________________________________________________________________

#http pos 'https://api.elsevier.com/content/search/scopus?query=title(cancer+inhibitor)+AND+openaccess(1)&apiKey=cd77a791f069b01e0d81c5a17c1c1890&1'
# cat tmp_search.json |grep "dc:title\|citedby-count" |paste - - |awk -F\" '{print $4"\t"$8}'  |sort -k1 -n -r |uniq

#API_KEY="cd77a791f069b01e0d81c5a17c1c1890" #OJO la api esta registrada como jpena@ucam.edu
API_KEY="220684153fe6faa1b6994202cfec8fe3" #	Hprez
#search=$1
#search="title(cancer+inhibitor)+AND+openaccess(1)"
#search="title(Humulone)"
#search="all(Humulone)"
#search="TITLE-ABS-KEY(receptor+diabetes)"
#out_put_file='tmp_search.json' 
search="TITLE-ABS-KEY(${1})"
out_put_file=$2
echo $search
echo $out_put_file



http post 'https://api.elsevier.com/content/search/scopus?query='${search}'&apiKey='${API_KEY} | python -m json.tool > ${out_put_file}
echo "https://api.elsevier.com/content/search/scopus?query='${search}'&apiKey='${API_KEY}"

items_per_page=`cat ${out_put_file}  |grep -i itemsPerPage |awk -F\" '{print $4}'`
total_results=`cat ${out_put_file}  |grep -i totalResults |awk -F\" '{print $4}'`
iterations=`echo "scale=0 ; $total_results / $items_per_page" | bc`
iterations=`expr $iterations + 1`
echo "Query: $search"
echo "Total Results: $total_results"
echo "Items Per Page: $items_per_page"
echo "Iterations: $iterations"

for i in `seq $iterations`;do
	http get "https://api.elsevier.com/content/search/scopus?query=${search}&apiKey=${API_KEY}&start=${i}" | python -m json.tool >> ${out_put_file}
	echo "http post \"https://api.elsevier.com/content/search/scopus?query=${search}&apiKey=${API_KEY}&start=${i}\" | python -m json.tool >> ${out_put_file}"
done

echo "Query: $search"
echo "Total Results: $total_results"
echo "Items Per Page: $items_per_page"
echo "Iterations: $iterations"
filename="${out_put_file%.*}"
 cat ${out_put_file} |grep "dc:title\|citedby-count" |paste - - |awk -F\" '{print $4"\t"$8}'  |sort -k1 -n -r |uniq > ${filename}_filter.json
