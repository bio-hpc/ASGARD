charmm < 01_minimize.inp >01_minimize.out
charmm < 02_energy_pdb.inp >02_energy_pdb.out
charmm < 02_energy_crd.inp >02_energy_crd.out
cp 01_minimize.pdb ../poly_pro_min_charmm.pdb
cp 01_minimize.crd ../poly_pro_min_charmm.crd
cp poly_pro_gas_all22.psf ../
charmm < 03_equil_crd.inp >03_equil_crd.out
cp 03_equil_crd.rst ../poly_pro_equil_charmm.rst
charmm < 04_nve_crd.inp >04_nve_crd.out

#clean up
rm -f *.pdb *.ene *.dcd *.rst *.dcd *.psf *.out *.crd
