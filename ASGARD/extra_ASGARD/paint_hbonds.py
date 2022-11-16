###________________________________________________________________________________________________________________________________________
###
###     Genera fichero pymol de hbonds entre proteina y ligando, se tienen en cuenta la distancia y los angulos.
###     python paint_hbonds.py proteina ligando y genera un .pml hbonds_proteina-ligando .pml
##      Muestra por pantalla los enlaces.
###_________________________________________________________________________________________________________________________________________

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys,  os
import pymol
from pymol import cmd
import numpy as np
import numpy.linalg as la
import math
from chempy import cpv
pymol.finish_launching()
paint_hbond=[]
count_hbond=0
seleccion=""
tablaVector1=[] #coordenadas ded x,y,z x,y,z y distancia del enlace
def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def get_round(i):
    #return round(i,3)
    return i

	
def get_atom(x,y,z):
	mymapLig = cmd.get_session("lig",1,1)['names'][0]
	mymapProt= cmd.get_session("prot",1,1)['names'][0]
	a=[]
	count=1
	countAtom=1
	atomo=""
	paint=""

	for a  in cmd.get_model('prot').atom:
		if a.coord[0]==x and a.coord[1]==y and a.coord[2]==z:
			atomo=a.resn+a.resi+"("+a.name+")-"
			paint=a.resi+" " +a.name
			
	for a  in cmd.get_model('lig').atom:
		if a.coord[0]==x and a.coord[1]==y and a.coord[2]==z:
			aa=""
			paint=a.id
			for numero in a.name:
				if not numero.isdigit():
					aa+=numero
			atomo="("+aa+str(a.id)+")" #+str(x) +" "+str(y) +" "+str(z)
	return atomo, paint

def get_results(r_object,state,names): #angleP o angleL
	a=[]
	astr=[]
	countador=1
	c=0
	global count_hbond
	global seleccion
	resFinal=""
	for obj in r_object: 
            try:
            	points = obj[5][2][state-1][1]
                for i in points:       
                    a.append(get_round(i)) #lista que se va llenando con las coordenadas x,y,z y x,y,z       
               
                    if countador % 6 ==0:
                    	for z in tablaVector1: 
                            #print a
                    		
                    		if z[3] == a[0] and z[4] ==a[1] and z[5]==a[2] and z[6] >= distanciaMinimaHidrogeno:	
                    			
                        		v1=[]
                        		v2=[]
                        		v1.append(get_round((z[0]-z[3])))
                        		v1.append(get_round((z[1]-z[4])))
                        		v1.append(get_round((z[2]-z[5])))
                        		v2.append(get_round((a[3]-a[0])))
                        		v2.append(get_round((a[4]-a[1])))
                        		v2.append(get_round((a[5]-a[2])))
                        		theta=py_ang(v1, v2)
                        		#print "hoasa"
                        		if theta>=math.radians(120):# and theta<=3.14: #2.44
                        		#print math.degrees(theta)
	                        	#	result
						if names=="HBA":
							a, pa =get_atom( z[3],z[4],z[5])
							result="["+a
							a, pb =get_atom( z[0],z[1],z[2])
							result+=a+","+str(round(z[6]))+","+str(int(math.degrees(theta)) )+"]"
							pa=pa.split()
							seleccion+=str(pa[0])+"+"


							paint_hbond.append("cmd.dist(\"hbond"+str(count_hbond)+"\",\"(resi "+str(pa[0])+" and name "+str(pa[1])+")\",\"(id "+str(pb)+" and lig)\")")
							count_hbond=count_hbond+1
		                                       # result, p="["+get_atom( z[3],z[4],z[5] ) # hidrogenos de proteina
                        		               # result, p=get_atom( z[0],z[1],z[2] )+","+str(z[6])+","+str(int(math.degrees(theta)) )+"]" # atomos pesados ligando
						else:
			                        	a, pa =get_atom( z[0],z[1],z[2])
			                        	
                                                        result="["+a
                                                        a, pb =get_atom( z[3],z[4],z[5])
                                                        result+=a+","+str(round(z[6]))+","+str(int(math.degrees(theta)) )+"]"
                                                        pa=pa.split()
                                                        seleccion+=str(pa[0])+"+"
							#print pa[0], pa[1]
                                                        paint_hbond.append("cmd.dist(\"hbond"+str(count_hbond)+"\",\"(resi "+str(pa[0])+" and name "+str(pa[1])+")\",\"(id "+str(pb)+" and lig)\")")
                                                        #print paint_hbond[count_hbond]
                                                        count_hbond=count_hbond+1
							#result="["+get_atom( z[0],z[1],z[2] ) #tarigo el hidrogen
			                        	#result+=get_atom( z[3],z[4],z[5] ) +","+str(z[6])+","+str(int(math.degrees(theta)))+"]"
			                        resFinal+=result+";"
	                        	

                    	a=[]
                    countador= countador+1
                if points is None:
                   raise ValueError
            except (KeyError, ValueError):
                continue
	return resFinal

def get_raw_distances(names='',  mol='', state=1, selection='all', quiet=1,resFinal=''):
    '''
DESCRIPTION
    Get the list of pair items from distance objects. Each list item is a
    tuple of (index1, index2, distance).
    Based on a script from Takanori Nakane, posted on pymol-users mailing list.
    http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10143.html
ARGUMENTS
    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}
    state = integer: object state {default: 1}
    selection = string: atom selection {default: all}
SEE ALSO
    select_distances, cmd.find_pairs, cmd.get_raw_alignment
    '''
	
    state, quiet = int(state), int(quiet)

    if state < 1:
        state = cmd.get_state()
    valid_names = cmd.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                #print ' Error: no such distance object:', name
                raise CmdException

    mol = cmd.get_names_of_type('object:molecule')             
    raw_objects = cmd.get_session(names, 1, 1, 0, 0)['names'] ##aqui casca y funciona a veces.    
    xyz2idx = {}
    cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,index)',  space=locals())
    c2=[]# atomo del aminoacido dos (debe ser siempre un H)
    a=[] # coord de los atomos   
    contadorPuentes=0;
    r = []

    countador=1
    counterAnglep=1
    for obj in raw_objects:

        try:
            points = obj[5][2][state - 1][1]
            for i in points:
                a.append(get_round(i))
                if countador % 6 ==0:
                    tablaVector1.append([])
                    tablaVector1[contadorPuentes].append(a[0])
                    tablaVector1[contadorPuentes].append(a[1])
                    tablaVector1[contadorPuentes].append(a[2])
                    tablaVector1[contadorPuentes].append(a[3])
                    tablaVector1[contadorPuentes].append(a[4])
                    tablaVector1[contadorPuentes].append(a[5])
                    contadorPuentes+=1;
                    a=[]
                countador= countador+1

            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        contadorPuentes=0;
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                tablaVector1[contadorPuentes].append(round ( cpv.distance (xyz1, xyz2 ) , 2 ) )
                contadorPuentes+=1
                #r.append((xyz2idx[xyz1], xyz2idx[xyz2], cpv.distance(xyz1, xyz2)))
                if not quiet:
                    print ' get_raw_distances:', r[-1]
                    a.append(r[-1])
            except KeyError:
                if quiet < 0:
                    print ' Debug: no index for', xyz1, xyz2
		
        if anglepF==0:
        	resFinal=get_results(cmd.get_session("anglep", 1, 1, 0, 0)['names'],state,names)           
        if anglelF==0:
        	resFinal+=get_results(cmd.get_session("anglel", 1, 1, 0, 0)['names'],state,names)  
      	#print resFinal
      	#exit()

        return resFinal
###_________________________________________________________________________________________________________________
###
###						MAIN
###______________________________________________________________________________________________________________________
distanciadeEnlaceMaxima=1.1
distanciaPuenteHidrogeno=2.5
distanciaMinimaHidrogeno=1.5
resFinal=""
prot = sys.argv[1]					#fichero de la proteina
lig = sys.argv[2]					#fichero del ligando
hbond_file="hbond_"+os.path.splitext(prot)[0]+"-"+os.path.splitext(lig)[0]+".pml"
#print hbond_file
anglepF=0
anglelF=0
HBAF=0
HBDF=0
#cargo proteina y ligando
pymol.cmd.load(prot,"prot")
pymol.cmd.load(lig, "lig")
##aniado hydros a la proteina
#pymol.cmd.h_add("prot & (don.|acc.)")
#pymol.cmd.h_add("lig & (don.|acc.)")

pymol.cmd.select("don","(hydro and (neighbor elem n,f,s,o))")
pymol.cmd.select("acc","(elem o or (elem n and not (neighbor hydro)))")
pymol.cmd.select("donor","(elem n,o,f,s and (neighbor hydro))")
##selecciono los hidrogenos de los donores de la prot y el acepotor del ligando
pymol.cmd.dist("HBA", "(\"lig\" and \"acc\")","(\"prot\" and \"don\")", distanciaPuenteHidrogeno)
##selecciono los hidrogenos de los donores de lig y el acepotor de la prot
pymol.cmd.dist("HBD", "(\"prot\" and \"acc\")","(\"lig\" and \"don\")" ,distanciaPuenteHidrogeno)

#selecciono los donores del ligando H y su corresponiente atomo para hacer el puente
pymol.cmd.dist("anglel","(\"lig\" and \"don\")","(\"lig\" and \"donor\")", distanciadeEnlaceMaxima)
#selecciono los donores de la proteinaH y su corresponiente atomo para hacer el puente
pymol.cmd.dist("anglep","(\"prot\" and \"don\")","(\"prot\" and \"donor\")", distanciadeEnlaceMaxima)
state=0

##
##apanio se podria mejorar es apra saber si las selecciones contienen
##	el 0 es correcto el 1 fallo
try: 
	get_results(cmd.get_session("anglep", 1, 1, 0, 0)['names'],state,"HBA")
	anglepF=0
except:
	anglepF=1
try: 
	get_results(cmd.get_session("anglel", 1, 1, 0, 0)['names'],state,"HBA")
	anglelF=0
except:
	anglelF=1
try: 
	resFinal=get_raw_distances("HBA",resFinal)
	#pymol.cmd.count_atoms("HBA")
	#HBAF=0
except:
	print "no se encuentran HBA"

tablaVector1=[]
	#HBAF=1
try: 
	resFinal+=get_raw_distances("HBD",resFinal)
	#pymol.cmd.count_atoms("HBD")
	#HBDF=0
except:
	print "no se encuentra HBD"
	#HBDF=1

#print "Selecciones : ",HBAF, HBDF, anglepF, anglelF
###
###________________________________________________________________________________________________
###
#if HBAF==0:
#	resFinal=get_raw_distances("HBA",resFinal)
#tablaVector1=[]
#print "dadas"+ str( HBDF)
#if HBDF== 0:
#	resFinal+=get_raw_distances("HBD",resFinal)

#print "ErTAS SEGURA?????"

#resFinal+=get_raw_distances("HBD",resFinal)

print resFinal #resFinal[:-2] #quitop el ultimo ;
f = open(hbond_file,"w")
f.write("cmd.load(\""+prot+"\", \"prot\")\n")
f.write("cmd.load(\""+lig+"\", \"lig\")\n")
for i in paint_hbond:
	f.write(i+"\n")
f.write("select aa, resi "+seleccion[:-1]+" or lig\n")
f.write("create amino, aa\n")
f.write("cmd.show(\"sticks\"    ,\"amino\")\n")
f.write("cmd.label(\'\'\'(name CA+C1*+C1\' and (byres(amino)))\'\'\',\'\'\'\"%s-%s\"%(resn,resi)\'\'\')\n")
f.write("hide lines\n")
f.write("cmd.set(\"valence\", 1)\n")
f.write("cmd.set(\"stick_radius\", 0.2)\n")
f.write("orient amino\n")
f.write("bg_color white\n")
f.write("set depth_cue=0\n")
f.write("set ray_trace_fog=0\n")
f.write("cmd.ray()\n")
#f.write("time.sleep(1)\n")
f.write("png "+hbond_file+".png, width=1200, height=1200, dpi=300, ray=1\n")
#f.write("cmd.quit()\n")
f.close()

pymol.cmd.save("prueba.pse")
pymol.cmd.quit()



"""
                        	if z[3] == a[0]: # and z[4] ==a[1] and z[5]==a[2]:	
                        		v1=[]
                        		v2=[]
                        		v1.append(get_round((z[0]-z[3])))
                        		v1.append(get_round((z[1]-z[4])))
                        		v1.append(get_round((z[2]-z[5])))
                        		v2.append(get_round((a[3]-a[0])))
                        		v2.append(get_round((a[4]-a[1])))
                        		v2.append(get_round((a[5]-a[2])))
                        		theta=py_ang(v1, v2)
                        		if theta>=2.44 and theta<=3.14:
                        			print theta 
                        			print "HBA: ", z[0], z[1], z[2], z[3], z[4], z[5], z[6], "angleP ", a[3], a[4], a[5]
                                   #print "cat "+mol[2]+".tmp"+"|grep "+get_str(a[0])+"|grep "+ get_str(a[1]) +"|grep "+ get_str(a[2]) +"|awk -F: '{print $6 $7}'"

"""


#dir = os.getcwd()					#directorio de trabajo
#pext = os.path.splitext(prot)[0]	#proteina sin extensuion
#lext = os.path.splitext(lig)[0]		#ligando sin extension
#dirp = dir+"/"+pext+"_h.pdb"		#aniade la ruta 
#prot_h = pext+"_h.pdb"
#pymol.cmd.save(dirp,"prot")
#pymol.cmd.delete("prot")
#pymol.cmd.load(prot_h,"prot")
#pymol.cmd.load(prot_h)
#pymol.cmd.load(lig)
#coordsAnglep=[]
#c=[] # atomo de los aminoacidos que intervine en el hbond
#l=[] # atomo del ligando interviene en el hbondo
#p=[] # cod de tres letras del aa interviene en hbond
##ordena la prot y el lig
#pymol.cmd.sort("prot extend 1")
#pymol.cmd.sort("lig extend 1")
#pymol.cmd.save("prot.pdb","prot")
#pymol.cmd.save("lig.pdb","lig")
