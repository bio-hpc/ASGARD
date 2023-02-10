#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Generate protein-ligand rmsd graphs if they exist
#
#
import os ,re
import numpy as np
import shutil
"""
class GraphsGeneral(object):
    def __init__(self, executeComand,cfg):
        self.cfg=cfg
        self.execute=executeComand
    #
	#	Se generan los ficheros de coordenadas xvg para las graficas de outData (RAPIDA)
	#
    def generarXvg(self,num,gr, out):
        comando="echo \""+num+" "+num+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+gr+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+out
        self.execute.run(comando)
    def generateXvgAnalyze(self,inp, out, step):
		comando=self.cfg.gromacs+" "+self.cfg.graph+"analyze"+self.cfg.mpi+" -f "+inp+" -dist "+out +" -bw "+str(step)
		self.execute.run(comando)
	#
	#	Genera el xvg del rmsf por residuo, se le indica el numero de grupo y la salida
	#	
	def gereateRmsdF(self,num, out)	:
		if self.cfg.perfil!="TWOLIGS":
			comando="echo \""+num+" "+num+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+"rmsf"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+out+" -res"
			self.execute.run(comando)
		else:
			comando="echo \""+num+" "+num+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+"rmsf"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+out
			self.execute.run(comando)
	#
	#	Genera Graficas de la proteina y el ligando, rmsd, distribucion, fluctuacion
	#
	def generateGraphProtinLigand(self, inp1, inp2,out,titulo):
		comando=self.cfg.python+" "+self.cfg.grapCALig+" "+inp1 +" "+inp2+" "+out+" "+titulo
		self.execute.run(comando)
	def graficaStander(self,inpu, titulo):
		comando=self.cfg.python+" "+self.cfg.graficaStandar+" "+ inpu +" "+titulo
		self.execute.run(comando)
	def graficaDistribucion(self,fichero,subTitleX,subtitlY,titulo):
		comando=self.cfg.python +" "+self.cfg.GrapAnalizeDistribution+" "+fichero + " \""+subTitleX+"\" "+"\""+subtitlY+"\" \""+titulo+"\""
		self.execute.run(comando)
	def grapAnalizeDistribution2Daata(self,fichero1,fichero2,out,subTitleX,subtitlY,titulo):
		comando=self.cfg.python +" "+self.cfg.grapAnalizeDistribution2Daata+" "+fichero1 +" "+fichero2 +" "+out+ " \""+subTitleX+"\" "+"\""+subtitlY+"\" \""+titulo+"\""
		self.execute.run(comando)

	#
	#	Genera grafica de RMSD y RMSDF de la proteina, si 
	#	
	def generateGraphRmsd(self):
		if self.cfg.graphRmsd:
			#
			#	Grafica prot Rmsd y RMSDF
			#
			outData=[self.cfg.outTxt+"_rmsdProt",self.cfg.outTxt+"_rmsdFProt"]
			self.generarXvg(self.cfg.protein[self.cfg.nomProtein],"rms", outData[0]+"R.xvg")									 #genera el fichero xvg de rms para de la proteina			
			self.generateXvgAnalyze(outData[0]+"R.xvg", outData[0]+"A.xvg", 0.001)									 #genera el fichero xvg de distribucion de la proteina llmado ...A.xvg
			self.gereateRmsdF(self.cfg.protein[self.cfg.nomProtein],outData[1]+"R.xvg")										 #genera el fichero xvg de rmdfF por residuo de la proteina
			self.graficaStander(outData[1]+"R.xvg", "\"RMSD Fluctuacion: "+self.cfg.nomProtein+"\"")				 #genera grafica fluctuacion proteina
			
			if len(self.cfg.ligandos)>0:																			 #si existen ligandos
				for nomLig,gnumLig  in self.cfg.ligandos.items():													 #para cada ligando
					
					out=self.cfg.outTxt+"_"+nomLig+"_rmsdLig"														 #salida rmsd
					outF=self.cfg.outTxt+"_"+nomLig+"_rmsdFLig"														 #sqilda Rmsdf
					self.generarXvg( gnumLig ,"rms",out+"R.xvg") 									 				 #genera el fichero xvg de rms para ligando
					self.generateXvgAnalyze(out+"R.xvg", out+"A.xvg", 0.001)						 				 #genera el fichero xvg de poblacion para ligando
					self.generateGraphProtinLigand( outData[0]+"R.xvg", out+"R.xvg",out+"_Protein_ligand.png","\"RMSD "+self.cfg.nomProtein+" vs "+nomLig+"\"") 				#graficca proteina ligando rmsd
					self.grapAnalizeDistribution2Daata( outData[0]+"A.xvg", out+"A.xvg",out+"_Protein_ligandMSD.png","Frecuecy","RMSD","Distribution "+self.cfg.nomProtein+" VS "+nomLig) 				#grafioca proteina ligando rmsd

					self.cfg.tools.joinImage(out+"_Protein_ligand.png",out+"_Protein_ligandMSD.png",out+"_Protein_ligand.png")	#unimos imagenes de la prot y el ligando
					comando="echo \""+gnumLig+" "+gnumLig+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+"rmsf"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+outF+"A.xvg"#generamos grafica fluctuaicon lig
					self.execute.run(comando)																																				#
					self.graficaStander(outF+"A.xvg", "\"RMSD Fluctuacion: "+nomLig+"\"")										  #grafica rmsf del ligando
					self.cfg.tools.joinImage(outData[1]+"R.png",outF+"A.png",outF+"Prot.png")										  #unimos proteina con remdf del ligando
			else:
				os.rename(outData[1]+"R.png",outData[1]+".png")	#si no existe ligando renombramos la graifca para que se vea
				self.graficaStander(outData[0]+"R.xvg","\"RMSD "+self.cfg.nomProtein+"\"" )
				self.graficaDistribucion(outData[0]+"A.xvg","RMSD","Distancia" ,"Distribution RMSD "+self.cfg.nomProtein )
				self.cfg.tools.joinImage(outData[0]+"R.png",outData[0]+"A.png",self.cfg.outTxt+"_rmsdProt.png")
							
	#
	#	Genera grafica de helices (RAPIDA)
	#
	def generateGraphHelix(self):
		if self.cfg.getHelix:
			helixfolder=self.cfg.outTxt+"Helix/"

			self.checkDirectorio(helixfolder)
			comando="echo "+self.cfg.protein[self.cfg.nomProtein]+" | "+self.cfg.gromacs+" "+self.cfg.graph+"helix"+self.cfg.mpi+ " -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -n "+self.cfg.grids+"_index_prot_lig_sol.ndx"
			self.execute.run(comando)
			os.rename("helicity.xvg",self.cfg.outTxt+"_Helix.xvg")
			os.chdir(self.cfg.dirActual)
			#shutil.rmtree(helixfolder)  #decomentar cuando no sea TMI
			comando=self.cfg.python+" "+self.cfg.graoGR+" "+self.cfg.outTxt+"_Helix.xvg"
			self.execute.run(comando)
			#################################################################################
			#	Grafica auxiliar solo para TMI Borrar Para otrsa dinamicas
			#
			filename, file_extension = os.path.splitext(self.cfg.graoGR)
			comando=self.cfg.python+" "+filename+"Aux.py "+self.cfg.outTxt+"_Helix.xvg 110 140" 
			self.execute.run(comando)
			comando=self.cfg.python+" "+filename+"Aux.py "+self.cfg.outTxt+"_Helix.xvg 385 415" 
			self.execute.run(comando)
			self.cfg.tools.joinImage(self.cfg.outTxt+"_Helix"+"_110_140_TMI.png" ,self.cfg.outTxt+"_Helix"+"_385_415_TMI.png" ,self.cfg.outTxt+"_Helix_TMI.png" )
			shutil.rmtree(helixfolder)  
			###################################################################################
	#
	#	(RAPIDA)
	#
	def graficaGyrate(self):	
		if (self.cfg.getGyrate):
			comando="echo "+self.cfg.protein[self.cfg.nomProtein] +"| "+self.cfg.gromacs+" "+self.cfg.graph+"gyrate"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+self.cfg.outTxt+"_Gyrate.xvg"
			self.execute.run(comando)
			comando=self.cfg.python+" "+self.cfg.graoGR+" "+self.cfg.outTxt+"_Gyrate.xvg"
			self.execute.run(comando)
	#
	#	Grafica de energia proteina ligando (RAPIDA)
	#
	def graficaEnergiaLigProt(self): #no se sua
		if self.cfg.energiLigProt:#la dejamos para el final
			pass
			#print "Vacio"
			#if self.cfg.nomLigando!="":
			#	comando="echo \""+self.cfg.gNumProt+" "+self.cfg.nLig+"\" |"+self.cfg.gromacs+" "+self.cfg.graph+"energy"+self.cfg.mpi+ " -f "+self.cfg.edrSimulacion+" -s "+self.cfg.tprMin+" -o "+self.cfg.outTxt+"_EnergyGlobalProtLig.xvg  < "+self.cfg.outTxt+"EnergyGlobalProtLig.txt" 
			#	self.execute.run(comando)
			#	os.remove(self.cfg.outTxt+"EnergyGlobalProtLig.txt")
			#	comando=self.cfg.python+" "+self.cfg.graficaStandar+" "+ self.cfg.outTxt+"_EnergyGlobalProtLig.xvg"
			#	self.execute.run(comando)
	#
	#	Grafica centro de masas distancia prot lig , tambien genera la distribucion
	#
	def graficaDist(self): 

		if self.cfg.distance: 
			for nomLig,gnumLig  in self.cfg.ligandos.items():
				if self.cfg.versionGromcas=="4":
					comando="echo \""+self.cfg.protein[self.cfg.nomProtein]+" "+gnumLig+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+"dist"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg"
				else:
					comando=self.cfg.gromacs+" "+self.cfg.graph+"distance"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -oall "+self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg -select 'com of group "+str(self.cfg.protein[self.cfg.nomProtein])+" plus com of group "+str(gnumLig)+" ' "			
				self.execute.run(comando)
				comando=self.cfg.python+" "+self.cfg.grapDistance+" "+self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg \": "+self.cfg.nomProtein+" "+ nomLig+"\" "
				self.execute.run(comando)
				self.generateXvgAnalyze(self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg", self.cfg.outTxt+"_"+nomLig+"_DistanceB.xvg",0.001)
				#comando=self.cfg.python +" "+self.cfg.grapAnalizeDistance+" "+self.cfg.outTxt+"_"+nomLig+"_DistanceB.xvg "+self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg \"Distribution: "+self.cfg.nomProtein+" vs "+nomLig+"\"" 
				self.graficaDistribucion(self.cfg.outTxt+"_"+nomLig+"_DistanceB.xvg","Frequency","Distance",self.cfg.nomProtein+" vs "+nomLig)
				self.cfg.tools.joinImage(self.cfg.outTxt+"_"+nomLig+"_DistanceA.png",self.cfg.outTxt+"_"+nomLig+"_DistanceB.png", self.cfg.outTxt+"_"+nomLig+"_Distance.png" )
	
	#
	#	TableMultiLigand
	#	Genera una tabla  con los datos de media, varianza del ligando con la proteina, se queda sin guardar por que no estoy seguro que lo usemos mucho
	#
	def generateTableMultiLigand(self,ftmp):
		f=open(self.cfg.outTxt+"_TableMultiligand.tex","w")
		f.write('\\begin{center} '+"\n")
		f.write('\\begin{tabular}{ l r r r r r r r r r r }'+"\n")
		f.write('\multicolumn{11}{c} {Distance Protein} \\\\'+"\n")
		f.write('\hline \n' )
		#f.write('Receptor  &  Ligand  &    Mean & Variance  &  Inicio   &  Final    &    7A   &     6A    &    5A   &     4A    &    3A\\\\'+"\n")
		firstLine=True
		for i in ftmp:

			f.write(i+"\\\\ \n")
			if firstLine:
				f.write('\hline \n' )
				firstLine=False
		f.write('\hline \n')
		f.write('\\end{tabular}'+"\n")
		f.write('\\end{center}'+"\n")

		f.close()


	#
	#	Genera una tabla tex para ver la aproximacion de los ligandos con la proteina
	#	
	def tableMultiliLiigand(self):
		if self.cfg.multipleLigand:
			datosPrint=[]
			formatTable='{0:>10}&{1:>10}&{2:>10}&{3:>10}&{4:>10}&{5:>10}&{6:>10}&{7:>10}&{8:>10}&{9:>10}&{10:>10}'
			datosPrint.append(formatTable.format("Receptor","Ligand","Mean","Variance","Inicio","Final", "7A","6A","5A","4A","3A" ))
			#
			#	Para cada ligando se calcula la distancia varizanza media inicio final
			#
			for nomLig,gnumLig  in self.cfg.ligandos.items():
				#ESTO NO HACE FALTA SI PREVIAMENTE SE HAN GENERADO LAS GRAFICAS DE DISTANCIA
				if not os.path.exists(self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg"):
					if self.cfg.gromacs=="":
						comando="echo \""+self.cfg.protein[self.cfg.nomProtein]+" "+gnumLig+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+"dist"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg"
					else:
						comando=self.cfg.gromacs+" "+self.cfg.graph+"distance"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -oall "+self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg -select 'com of group "+str(self.cfg.protein[self.cfg.nomProtein])+" plus com of group "+str(gnumLig)+" ' "
					self.execute.run(comando)
				distanciaA=[]
				for cont in range(5):
					d="\times"
					for i in self.cfg.ligandosCoords[nomLig]:
						for j in self.cfg.protinCoords:
							auxLig=re.sub(' +',' ',i).strip().split(" ")
							auxProt=re.sub(' +',' ',j).strip().split(" ")
							if len(auxProt)>3 and auxProt[0]!="TITLE" and auxProt[0]!="REMARK" and auxProt[0]!="CRYST1":
								dist=self.cfg.tools.distancia([float(auxLig[5]),float(auxLig[6]),float(auxLig[7])],[float(auxProt[5]),float(auxProt[6]),float(auxProt[7])]  )
								#print self.cfg.distanceLigadProt-cont
								if dist < self.cfg.distanceLigadProt-cont:
									d="\checkmark"
									break
					distanciaA.append(d)
				ftmp=open(self.cfg.outTxt+"_"+nomLig+"_DistanceA.xvg")
				auxTmp=[]
				for i in ftmp:
					if not i.startswith("@") and not i.startswith("#"):
						aux=re.sub(' +',' ',i).strip().split(" ")
						auxTmp.append(float(aux[1]))
				ftmp.close()
				
				datosPrint.append(formatTable.format(self.cfg.nomProtein,nomLig,round(np.mean(auxTmp),3),round(np.var(auxTmp),3),round(auxTmp[0],3),round(auxTmp[len(auxTmp)-1],3),distanciaA[0],distanciaA[1],distanciaA[2],distanciaA[3],distanciaA[3]  ) )
			self.generateTableMultiLigand(datosPrint)	
			#print "datos print"
			#for i in datosPrint:
			#	print i
			
	#
	#Borra el directorio si existe y lo vuelve a crear y se introduce en el
	#	
	def checkDirectorio(self, path):
		if os.path.exists(path):
			shutil.rmtree(path)
			os.mkdir(path)
		else:
			os.mkdir(path)
		os.chdir(path)
		

	


	#def standarGraph(self, inpu):
	#	comando=self.cfg.python+" "+self.cfg.graficaStandar+" "+ inpu
	#	self.execute.run(comando)
######comando="echo \""+self.cfg.protein["Protein"]+" "+gnumLig+"\" | "+self.cfg.gromacs+" "+self.cfg.graph+"dist"+self.cfg.mpi+" -f "+self.cfg.xtcSimu+" -s "+self.cfg.tprMin+" -o "+self.cfg.outTxt+"_Distance_Protein_"+nomLig+".xvg"

#python_packages_path = os.path.realpath("/mnt/scratch/users/ac_001_um/lanzador/lanzador//externalSw/python2-packages/") #esto debe ser más dinamico
#sys.path.append(python_packages_path)
#try:
#    import seaborn
#except ImportError:
#    pass
#class GenerateGraphs(object): 
"""generateGraphProtinLigand