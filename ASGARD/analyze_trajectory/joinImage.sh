#!/bin/sh
###_________________________________________________________________________________________________________---
###
###	Genera imagenes en una tabla de X*2 y lo guarda en montage.png
###	ejemplo: ./joinImage.sh imageA imageB
###__________________________________________________________________________________________________________
letraImagen=65 #letra 65 numero 1
iamgenes=""
conatorFilas=0
filas=0
out=""
columnas=2

while [ $1 ];do
	if [[ $1 != *"Estabilizacion"* ]] && [[ $1 != *"RMSDS"* ]];then
		#letra=$letraImagen # NUMEROS
		letra=`echo $letraImagen | awk '{printf("%c",$1)}'` #letra que se le asignara, en orden de entrada
		convert -size 2048x1536 ${1} -thumbnail 800x600 ${1}${letra}A.png									#redimensionando la imagen para que todas sean iuales
		convert ${1}${letra}A.png  -background transparent -gravity center -extent 840x640 ${1}${letra}B.png	 #pone el tamaño del lienzo en 800x600
	#	convert ${letra}B.png -gravity northwest -background YellowGreen  -pointsize 40 -annotate +420+600  ${letra} ${letra}.png #se coloca la letra
		convert ${1}${letra}B.png -gravity northwest -background YellowGreen  -pointsize 30 -annotate +50+50  ${letra} ${1}${letra}.png #se coloca la letra
	
		rm ${1}${letra}A.png		#se borra la imagen Aux
		rm ${1}${letra}B.png		#se borra la imagen Aux
		imagenes=$imagenes" "${1}${letra}.png
		letraImagen=`expr $letraImagen + 1`
		conatorFilas=`expr $conatorFilas + 1`
	else
		out=$1
	fi

	shift #paso aprametro 
done
if [ `echo "$conatorFilas % 2" | bc` -ne 0 ]; then 
	conatorFilas=`expr $conatorFilas + 1`
fi
echo $out
filas=$((conatorFilas / columnas))

if [ `echo "$conatorFilas % $columnas" |bc` -ne 0 ];then #si la dividison tiene resto quiere decir que necesita una fila mas
	filas=`expr $filas + 1`
fi
montage  -geometry +3+2  $imagenes -tile ${columnas}x${filas} $out".jpeg"
convert -resize 620x $out.jpeg $out-620.jpeg
convert -resize 1024x $out.jpeg $out-1024.jpeg
rm $imagenes
