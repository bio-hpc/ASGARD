#! /usr/bin/python
import sys
import os
from numpy import *
import Scientific.IO.NetCDF
from Numeric import * 
from Scientific.IO import NetCDF as Net
from ReadAmberFiles import *

##################################
#ARGUMENTS

#these are the dimensions of the supercell (what PropPDB was given)
ix=3
iy=2	
iz=2

# number of asymmetric units in unit cell
asymunits=1

#this is the transpose of the invU matrix (xyz components of an "x", "y", "z" propagation of the unit cell)
# obsolete with the CompXfrm python routine introduce below
#~ xv=array([94.1000, 0, 0])
#~ yv=array([-47.0500, 81.4930, 0])
#~ zv=array([0., 0., 132.4])

#directory with the split trajectories
trajpath='../splittrajectories/'

# box dimensions of the crystal supercell. If pressure scaling is isotropic (as this script assumes) you really
# only need the first element (the a or x dimension of the supercell)
A=[81.72, 63.74, 68.46, 88.52 ,108.53 ,111.89] 
pdbwsymtry='../4lzt.pdb'




############################################
# calculate the xyz components of an "x", "y", "z" propagation of the unit cell
pdbfile=pdb(pdbwsymtry)
UCbox=pdbfile.Get_Box()
u,invu=CompXfrm(UCbox)
xv=transpose(invu)[0,:]
yv=transpose(invu)[1,:]
zv=transpose(invu)[2,:]


# Get symmetry rotation matrices (smtry) and translation vectors (tr)
smtry,tr =pdbfile.Get_SMTRY()

print 'Number of asymunits provide: %d' %asymunits
print 'Number of symmetry operations in pdb: %d' %(len(smtry))

i=1 #counter for the names of the unitcell trajectories
#for each unit cell in the supercell
for x in range(ix):
	for y in range(iy):
		for z in range(iz):
			#calculate the move vector that was applied
			xmove=(x*(xv[0])+y*(yv[0])+z*(zv[0]))
			ymove=(x*(xv[1])+y*(yv[1])+z*(zv[1]))
			zmove=(x*(xv[2])+y*(yv[2])+z*(zv[2]))
			MoveVector=array([xmove,ymove,zmove],dtype=float32)

			for j in range(asymunits):			
				#get symmetry operations for this asym unit
				s=float32(transpose(smtry[j]))
				t=float32(tr[j])
				
				print '\n\nStarting file %02d_%02d.nc' %(i,j+1)
				print 'Move vector is '+ str(MoveVector)
				print 'Symmetry rotation matrix is:'
				print s
				print 'Symmetry tranlation vector is %s' %(str(t))

				#copy the unit cell trajectory to a new taget file
				#read the netcdf file of each unitcell trajectory
				#coords stores the entire coordinates variable which is [frames,atoms,xyzcoords]
				#frames stores the number of frames in the trajectory
				orig_filename=trajpath+'%02d_%02d.nc' %(i,j+1)
				filename='RevSymm_%02d_%02d.nc' %(i,j+1)
				os.system('cp %s %s' %(orig_filename, filename))
				ofile = Net.NetCDFFile(filename, 'a')
				coords=ofile.variables['coordinates']
				frames=coords.shape[0]
				celllen=ofile.variables['cell_lengths']
				
				#iterate over each frame: rescale the atomic coordinates according to boxsize pressure scaling
				for frame in range(frames):
					a=celllen[frame,0] #get a-vector length in current frame
					scale=A[0]/a        
					coords[frame,:,:]=coords[frame,:,:]*scale

				#############################
				# TEST setup to check anisotropic pressure scaling impact
				#############################

				#~ for frame in range(frames):
					#~ a=celllen[frame,0] #get a-vector length in current frame
					#~ scale=A[0]/a        
					#~ coords[frame,:,0]=coords[frame,:,0]*scale
#~ 
					#~ a=celllen[frame,1] #get a-vector length in current frame
					#~ scale=A[1]/a        
					#~ coords[frame,:,1]=coords[frame,:,1]*scale					
					#~ 
					#~ a=celllen[frame,2] #get a-vector length in current frame
					#~ scale=A[2]/a        
					#~ coords[frame,:,2]=coords[frame,:,2]*scale					
					
				#reverse translate to original unit cell by subtracting the translation vector
				coords[:,:,:]=coords[:,:,:]-MoveVector
				#reverse symmetry operate to original asym unit by applying the symmetry translation and rotation
				coords[:,:,:]=dot( (coords[:,:,:]-t),linalg.inv(s) )
				#close file (write to netcdf file)
				ofile.close()
			i+=1





