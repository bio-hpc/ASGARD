import sys, commands, math
#
#	Mide la distnacia maxima del ligando busca x1, y1, z1 maximo y x1, y1, z1 minimo y calcula la diagonal
#	para generar cubos en pymol
#	depende de:
#		standarFileCoords.py
#	Se le pasa x, y ,z y devueleve los residudos cerca de 3 A de la proteina
#
StandarCoords="lanzador/scriptsLanzador/standarFileCoords.py"
def distance(a,b):
	x=math.pow((a[0]-b[0]),2) 
	y=math.pow((a[1]-b[1]),2)
	z=math.pow(( a[2]-b[2]),2)
	suma=x+y+z;
	return round(math.sqrt(suma),3 )

if len(sys.argv)!=5:
	print "debe introducir proteina x, y ,z"
	exit()
fichero=sys.argv[1]
centerX=sys.argv[2]
centerY=sys.argv[3]
centerZ=sys.argv[4]
comando="python "+StandarCoords+" "+fichero
lineas=commands.getoutput(comando)
x=[]
y=[]
z=[]
lineas=lineas.split(" ")
for f in lineas:
	print f
"""	
	aux=lineas.split(":")
	x.append(float(aux[0]))
	y.append(float(aux[1]))
	z.append(float(aux[2]))
distanciaTotal=0;
for i in range(len(x)):
	a=[x[i],y[i],z[i]]
	for j in range(len(x)):
		b=[x[j],y[j],z[j]]
		d=distance(a,b)
		if d>distanciaTotal:
			distanciaTotal=d
comando="rm "+fichero+"TMP"
commands.getoutput(comando)
print distanciaTotal
"""
