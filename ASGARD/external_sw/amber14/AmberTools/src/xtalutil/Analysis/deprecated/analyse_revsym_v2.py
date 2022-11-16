#! /usr/bin/python
import sys
import os
from numpy import *
from ReadAmberFiles import *

# SET ARGUMENTS
unitcells=12
asymunits=1
topo='../UC.prmtop'
supercell_trajectory='../mergetraj_centonpdb_nowat.nc'
pdb='../UC_centonpdb.rst7'
# NEED TO CHANGE ATOM SELECTIONS FOR RMSD AND B-FACTOR


########################################################################

# Get Frames
traj=nc(supercell_trajectory)
frames=traj.Get_Frames()

###############
#   Calculate B-factors
###############

f=open('ptraj_bfactor','w')
for i in range(unitcells):
	for j in range(asymunits):
		if i <=8:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		else:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))

f.write('atomicfluct out bfac_revsymm_calpha.dat :1-129@CA&!(@H=) byatom bfactor\n')
f.write('atomicfluct out bfac_revsymm_res.dat :1-129&!(@H=) byres bfactor\n')
f.close()
os.system('ptraj '+topo+' <ptraj_bfactor')
 
 
###############
# Calculate B-factors after rmsd fitting unit cells
# (shown to produce same results as supertrajectory method
###############

f=open('ptraj_bfactor','w')
for i in range(unitcells):
	for j in range(asymunits):
		if i <=8:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		else:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
f.write('rms first mass \':2-123@C,CA,N,O\'\n')
f.write('atomicfluct out bfac_rmsd_calpha.dat :1-129@CA&!(@H=) byatom bfactor\n')
f.write('atomicfluct out bfac_rmsd_res.dat :1-129&!(@H=)  byres bfactor\n')
f.close()
os.system('ptraj '+topo+' <ptraj_bfactor')

#~ ###############
#~ # Calculate B-factors after com fitting unit cells
#~ ###############
#~ 
#~ f=open('ptraj_bfactor','w')
#~ for i in range(unitcells):
	#~ for j in range(asymunits):
		#~ if i <=8:
			#~ if j <=8:
				#~ f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			#~ else:
				#~ f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		#~ else:
			#~ if j <=8:
				#~ f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			#~ else:
				#~ f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
#~ f.write('center mass origin\n')
#~ f.write('atomicfluct out bfac_com_calpha.dat :1-124@CA|:125-126@C1\' byatom bfactor\n')
#~ f.write('atomicfluct out bfac_com_res.dat \':1-126 & !@H=\' byres bfactor\n')
#~ f.close()
#~ os.system('ptraj '+topo+' <ptraj_bfactor')


###############
#   Calculate Lattice rmsd
###############

f=open('ptraj_CalcLatRmsd','w')
#~ f.write('reference ../UC.pdb\n')
for i in range(unitcells):
	for j in range(asymunits):
		if i <=8:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		else:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
f.write('rms first mass out rmsd_lat_sdcn.dat :3-127&!(@H=,CA,C,O,N) nofit time 0.01 \n')
f.write('rms first mass out rmsd_lat_bbone.dat :3-127@CA,C,N nofit time 0.01\n')
f.close()
os.system('ptraj '+topo+' <ptraj_CalcLatRmsd')


###############
#   Calculate Monomer rmsd
###############

f=open('ptraj_CalcUCRmsd','w')
#f.write('reference ../UC.pdb\n')
for i in range(unitcells):
	for j in range(asymunits):
		if i <=8:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		else:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
f.write('rms first mass out rmsd_monomer_sdcn.dat :3-127&!(@H=,CA,C,O,N) time 0.01 \n')
f.write('rms first mass out rmsd_monomer_bbone.dat :3-127@CA,C,N time 0.01\n')
f.close()
os.system('ptraj '+topo+' <ptraj_CalcUCRmsd')

def ProcessRmsd(ifile,ofile,frames,total_monomers):
	data=genfromtxt(ifile)
	assert frames%total_monomers == 0, "number of frames not a multiple of total_monomers" % id
	x=array([])
	for i in range(frames):
		for j in range(total_monomers):
			x=append(x,data[i+j*frames,1])
		f.write('%7.2f  %7.5f  %7.5  %7.5\n' %(data[i,0], x,mean(x),std(x)))
		x=array([])
	f.close()

UCrmsd('rmsd_lat_bbone.dat','rmsd_lat_bbone_table.dat',frames,unitcells*asymunits)
UCrmsd('rmsd_lat_sdcn.dat','rmsd_lat_sdcn_table.dat',frames,unitcells*asymunits)
UCrmsd('rmsd_monomer_bbone.dat','rmsd_monomer_bbone_table.dat',frames,unitcells*asymunits)
UCrmsd('rmsd_monomer_sdcn.dat','rmsd_monomer_sdcn_table.dat',frames,unitcells*asymunits)
	
def UCrmsd(ifile, ofile, frames, total_monomers):
	f=open(ofile,'w')
	data=genfromtxt(ifile)
	x=array([])
	for i in range(frames):
		for j in range(total_monomers):
			x=append(x,data[i+j*frames,1])
		rmsd=sqrt(mean(x**2))
		f.write('%7.2f  %7.5f \n' %(data[i,0], rmsd))
		x=array([])
	f.close()

	
UCrmsd('rmsd_lat_bbone.dat','rmsd_lat_bbone_UC.dat',frames,unitcells*asymunits)
UCrmsd('rmsd_lat_sdcn.dat','rmsd_lat_sdcn_UC.dat',frames,unitcells*asymunits)
UCrmsd('rmsd_monomer_bbone.dat','rmsd_monomer_bbone_UC.dat',frames,unitcells*asymunits)
UCrmsd('rmsd_monomer_sdcn.dat','rmsd_monomer_sdcn_UC.dat',frames,unitcells*asymunits)


 
#################
# Calculate average structure using rmsd
#################

f=open('ptraj_CalcAvgesAvg','w')
f.write('reference %s\n' %pdb)
for i in range(unitcells):
	for j in range(asymunits):
		if i <=8:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		else:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
f.write('rms reference mass :3-127&!(@H=)\n')
f.write('average average_rmsd.mdcrd \n')
f.close()
os.system('ptraj '+topo+' <ptraj_CalcAvgesAvg')

#~ #################
#~ # Calculate average structure using com
#~ #################
#~ 
#~ f=open('ptraj_CalcAvgesAvg','w')
#~ f.write('reference %s\n' %pdb)
#~ for i in range(unitcells):
	#~ for j in range(asymunits):
		#~ if i <=8:
			#~ if j <=8:
				#~ f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			#~ else:
				#~ f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		#~ else:
			#~ if j <=8:
				#~ f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			#~ else:
				#~ f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
#~ f.write('center mass origin\n')
#~ f.write('average average_com.mdcrd \n')
#~ f.close()
#~ os.system('ptraj '+topo+' <ptraj_CalcAvgesAvg')

#################
# Calculate average structure using reverse symmetry
#################

f=open('ptraj_CalcAvgesAvg','w')
for i in range(unitcells):
	for j in range(asymunits):
		if i <=8:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
		else:
			if j <=8:
				f.write('trajin RevSymm_0%d_0%d.nc\n' %(i+1,j+1))
			else:
				f.write('trajin RevSymm_0%d_%d.nc\n' %(i+1,j+1))
f.write('average average_revsym.mdcrd \n')
f.close()
os.system('ptraj '+topo+' <ptraj_CalcAvgesAvg')


###############
# Calculate average structure rmsd's
###############

f=open('ptraj_rmsd','w')
f.write('reference %s\n' %pdb)
f.write('trajin average_rmsd.mdcrd\n')
f.write('rms reference mass out avg_rmsd_sdcn.dat :3-127&!(@H=,CA,C,O,N) time 0.01 \n')
f.write('rms reference mass out avg_rmsd_bbone.dat :3-127@CA,C,N time 0.01 \n')
f.close()
os.system('ptraj '+topo+' <ptraj_rmsd')

#~ f=open('ptraj_rmsd','w')
#~ f.write('reference %s\n' %pdb)
#~ f.write('trajin average_com.mdcrd\n')
#~ f.write('rms reference mass out avg_com_sdcn.dat \':1-124 & !(@H=,CA,C,O,N)\' time 0.01 \n')
#~ f.write('rms reference mass out avg_com_bbone.dat \':1-124@CA & !(@H=)\' time 0.01 \n')
#~ f.close()
#~ os.system('ptraj '+topo+' <ptraj_rmsd')

f=open('ptraj_rmsd','w')
f.write('reference %s\n' %pdb)
f.write('trajin average_revsym.mdcrd\n')
f.write('rms reference mass out avg_revsym_sdcn.dat :3-127&!(@H=,CA,C,O,N) time 0.01 \n')
f.write('rms reference mass out avg_revsym_bbone.dat :3-127@CA,C,N time 0.01 \n')
f.close()
os.system('ptraj '+topo+' <ptraj_rmsd')
