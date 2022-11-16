Comparison output are saved in save/ subdir.


make
  or
make serial
     tests:
	i) The log output from chamber
	ii) The prmtop and incprd generated from chamber
	iii) SANDER's energy and forces on this system via the charmm_gold output

make charmm
    tests: energy printout of charmm and sander match.

make gluc_energy.out
   requires charmm executable in your path or reset CHARMM in ../common.mk
   produces the energy output for checking, the psf and xplor psf
