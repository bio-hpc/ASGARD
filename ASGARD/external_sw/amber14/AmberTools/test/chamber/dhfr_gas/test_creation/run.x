charmm < 01_minimize.inp >01_minimize.out
charmm < 02_energy.inp >02_energy.out
cp 01_minimize.crd ../dhfr_min_charmm.crd
cp dhfr_gas_all22_prot.psf ../


rm *.out *.crd *.psf 01_minimize.pdb

