#! /usr/bin/python
import sys
import os
from numpy import *
import Scientific.IO.NetCDF
from Scientific.IO import NetCDF as Net
import time
import RunTime

#======================================================================#
#                                                                      #
# Prepare pdb frames from a trajectory for md2map.sh.                  #
# The trick is to make sure that column 13 of the pdb file has correct #
# element names. Run it once on just one frame (set frames=1) to get   #
# output of unrecognized element names. (For example, tip4p EPW's will #
# need to be removed.                                                  #
#                                                                      #
#======================================================================#


#======================================================================#
#
# SET VARIABLES
SC_TOPO='fav8_Ex2.prmtop'
SC_TRAJ='fit_wat.nc'
PDB_FILE='asu1.pdb'              #must contain Asym unit CRYST1 record
frames=5000                      #total no of PDB files to create (one for each frame)
startframe=1                     #first frame of trajectory to start with
offsetframe=2                    #no. of frames to skip between 
#
#======================================================================#

# get absolute paths and check existence of files
print "Using files:"
SC_TOPO=os.path.abspath(SC_TOPO)
SC_TRAJ=os.path.abspath(SC_TRAJ)
PDB_FILE=os.path.abspath(PDB_FILE)
for i in [SC_TOPO,SC_TRAJ,PDB_FILE]:
  if not os.path.isfile(i):
    print "ERROR: File %s does not exist." %i
    sys.exit()
  else:
    print "%s" %i 
print ""

#make directory
os.system('mkdir -p PDBData')
os.chdir('PDBData')

#get cell lenghts
ofile = Net.NetCDFFile('%s' %SC_TRAJ, 'r')
celllen=ofile.variables['cell_lengths']

# prepare the simulation pdbs
#GET ANGLES
f=open(PDB_FILE,'r')
while 1:
  line=f.readline()
  if line[0:6] == 'CRYST1':
    break
angles=" ".join(line.split()[4:])
    
start_time=time.time()
for j in range(frames):
	frame=startframe+offsetframe*j
	#REPORT TIME
	frac_complete=float(frame)/frames
	print "frame: %d" %frame
	print RunTime.main(start_time=start_time, frac_complete=frac_complete)
	#CPPTRAJ CALL TO GET PDB OF FRAME
	f=open('ctraj_frame', 'w')
	f.write('parm %s\n' %SC_TOPO)
	f.write('trajin %s %d %d 1\n' %(SC_TRAJ,frame,frame))
	f.write('trajout %04d.pdb pdb\n' %frame)
	f.close()
	os.system('cpptraj -i ctraj_frame >tmp')
	#ADD CELL PARAMETERS
	cell=celllen[frame-1]
	os.system('sed -i \'1i CRYST1   %6.3f   %6.3f   %6.3f %s  \' %04d.pdb' 
            %(cell[0],cell[1],cell[2],angles,frame))
	#REMOVE TIP4PEW VIRTUAL ATOM
	os.system('grep -v EPW %04d.pdb >tmp' %frame)
	os.system('mv tmp %04d.pdb' %frame)

#CLEAN UP
os.system('rm -rf ctraj_frame new.pdb tmp')	


