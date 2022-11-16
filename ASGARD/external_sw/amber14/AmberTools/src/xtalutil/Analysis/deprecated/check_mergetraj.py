#!/usr/bin/env python
import sys
import os
from numpy import *
from ReadAmberFiles import *


########################################################################
# After merging a trajectory, checks if trajectory is centered on pdb.
# Arguments:
#     argv[1] - parmtop topology
#     argv[2] - trajectory
#	  argv[3] - unit cell rst7 file (amber parmtop atom order but pdb crystal
#				coordinates
#
# Return:
#	  prints error or exits with status 0
#########################################################################



##### Check centering
topo=prmtop(sys.argv[1])
masses=topo.Get_Masses()

B=nc(sys.argv[2])
coords=B.Get_Traj(100)
com1=COM(coords[0,:,:],masses)
com2=COM(coords[1,:,:],masses)
com3=COM(coords[2,:,:],masses)
com4=COM(coords[3,:,:],masses)
for i in range(1,4):
	
	if round(max(abs(com1-COM(coords[i,:,:],masses))),2) >=0.001:
		print "ERROR: centers of mass do not agree" 
		print "Frame 0 center of mass: "+str(com1) 
		print ("Frame %d center of mass: " %(i*100)) + str(COM(coords[i,:,:],masses))
		sys.exit(1)
#### Check pdb coordinates

pdbfile=rst7(sys.argv[3])
coords_pdb=pdbfile.Get_Coords()
coords_traj=coords[0,:,:]
for i in range(1,40,10):
	
	if max(abs(coords_pdb[i]-coords_traj[i])) >=0.001:
		print "ERROR: frame 0 coordinates do not agree with pdb"
		print ("PDB Atom %d coordinates: " %(i+1)) + str(coords_pdb[i])
		print ("Frame 0 Atom %d coordinates: " %(i+1)) + str(coords_traj[i])
		sys.exit(1)
		
