#!/bin/sh

charmm < 01_minimize.inp >01_minimize.out

#Truncate crd file to 3d.p since Charmm does
#not seem to be able to read its own pdb file 
#so we have to use a crd.
./zero_3dp

charmm < 02_energy.inp >02_energy.out
charmm < 02_energy_shake.inp >02_energy_shake.out

cp 01_minimize.crd ../dhfr_min_charmm.crd
cp dhfr_cmap_pbc.psf ../

rm 01_minimize.out 01_minimize.crd *.psf 01_minimize.pdb initial_structure_from_charmm.pdb

