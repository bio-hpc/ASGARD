#! /usr/bin/python
import ReadAmberFilesLight as raf
import numpy as np
import phenix_amber_interface as ext
import argparse

#=========================================================================#
#                                                                         #
# A regression test for phenix_amber_interface. vAla3.pdb coordinates are #
# jiggled by a 0.1Ang 3D-Gaussian. Steepest descent is applied with a     #
# default step size of 0.0001 Ang for 1000 iteration or until either the  #
# change in total energy is less than 0.001 kcal/mol or the change in     #
# gradients rms is less than 0.001 kcal/mol/A. Finally steps taken to     #
# reach termination condition and coordinate RMSD between starting (pre-  #
# jiggle) and minimized structures is reported.                           #
#                                                                         #
#========================================================================#



#currently std.dev of 0.15Ang per dimension shakes with with std. dev. 0.1Ang in 3D
def jiggle(coords2):
  for atom in range(len(coords2)):
    coords2[atom]=coords2[atom]+(.15*np.random.randn(3))
  return  coords2

#get gradients and energies from mdgx
def get_gradients_and_target(sites_cart,gradients,target, U):
  sites_cart_c=ext.VectorOfDouble()
  sites_cart_c.extend(i for i in sites_cart)
  gradients_c=ext.VectorOfDouble()
  gradients_c.extend(i for i in gradients)
  target_c=ext.VectorOfDouble()
  target_c.extend(i for i in target)
  ext.callMdgx(sites_cart_c, gradients_c, target_c, U)
  gradients=np.array(gradients_c)
  target=np.array(target_c)

  return gradients, target

#apply gradients
def steepest_descent(sites_cart,gradients,step=0.0001):
  sites_cart+=(gradients*step)
  return sites_cart

def run(pdbfile,prmtopfile,crdfile):
  #READ PDB FILE
  apdb=raf.pdb(pdbfile)
  box=apdb.Get_Box()
  coords=apdb.Get_Coords()
  natoms=natoms= len(apdb.ATOMlines)

  #JIGGLE COORDINATES
  coords=jiggle(coords)

  #INSTATIATE MDGX STRUCTS AND ARRAYS
  U=ext.uform(prmtopfile, crdfile)
  sites_cart=coords.flatten()
  gradients=np.zeros(len(sites_cart))
  target=np.zeros(10)

  #GET INITIAL GRADIENTS AND ENERGIES
  old_gradients,old_target=get_gradients_and_target(sites_cart, gradients, target, U) 

  #STEEPEST DESCENT ITERATIONS
  i=0
  while i<100:
    sites_cart=steepest_descent(sites_cart, old_gradients)

    new_gradients, new_target= get_gradients_and_target(sites_cart, gradients, target, U)     

    target_diff=old_target[0]-new_target[0]
    delta_gradient_rms=old_target[9]-new_target[9]
    print("%7.3f   %7.3f   %7.3f" %(new_target[0], target_diff, delta_gradient_rms))
    if delta_gradient_rms<0.001 or target_diff<0.001:
      break
    
    old_gradients, old_target=new_gradients, new_target
    
    i+=1        

  print("Minimization termination condition reached after %d steps." %i)

  #CALCULATE RMSD
  newcoords=np.reshape(sites_cart, (-1,3))
  newcoords=raf.KabschAlign(coords, newcoords)
  rmsd=raf.RMSD(coords,newcoords)
  print("RMSD between original and minimized structure: %7.4f" %rmsd)

  # optional uncomment to print out coords in Amber Restart File
  #~ raf.printrst7(newcoords, box, natoms, 'new_coords.rst7')

if __name__ == "__main__" :
        parser = argparse.ArgumentParser()
        parser.add_argument("-pdb", help="name of pdb file", default="vAla3.pdb")
        parser.add_argument("-prmtop", help="name of topology file", default="vAla3.prmtop")
        parser.add_argument("-crd", help="name of coordinate file", default="vAla3.rst7")
        args = parser.parse_args()
        run(args.pdb,args.prmtop, args.crd)



########################################################################
#~ def getstd():
  #~ R=[]
  #~ for iter in range(100):
    #~ i=.15*np.random.randn(3)
    #~ r=np.sqrt(i[0]**2+i[1]**2+i[2]**2)
    #~ R.append(r)
  #~ return np.array(R).std()
  #~ 
#~ x=0
#~ for jiter in range(100):
  #~ x+=getstd()
#~ print x/100.0  
  
#~ print np.sqrt(3)
#~ print np.sqrt(1.0/3.0)

#~ mean=[0,0,0]
#~ cov=[[2.25,0,0],[0,2.25,0],[0,0,2.25]]
#~ x=0
#~ for i in range (1000):
  #~ sample=np.random.multivariate_normal(mean,cov,100)
  #~ x+=np.std(np.array([np.linalg.norm(r) for r in sample]))
#~ print x/1000


