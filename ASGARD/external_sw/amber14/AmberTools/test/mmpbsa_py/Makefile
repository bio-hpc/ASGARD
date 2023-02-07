# Makefile for the MMPBSA.py tests

test: COMPREHENSIVE
all: ANTE-MMPBSA GB PB ALA DECOMP NAB COMPREHENSIVE STABILITY QMMM 3DRISM

GB:
	cd 01_Generalized_Born && ./Run.GB

PB:
	cd 02_Poisson_Boltzmann && ./Run.PB

ALA:
	cd 03_Alanine_Scanning && ./Run.ALA

DECOMP:
	cd 04_Per_Residue_Decomp && ./Run.PerRes
	cd 05_Pairwise_Decomp && ./Run.Pairwise

NAB:
	cd 06_NAB_Nmode && ./Run.nmode

COMPREHENSIVE:
	cd 07_Comprehensive && ./Run.comprehensive

STABILITY:
	cd 08_Stability && ./Run.stability

APBS:
	cd 09_APBS_Poisson_Boltzmann && ./Run.PB

QMMM:
	cd 10_QM_MMGBSA && ./Run.QMMMGBSA

3DRISM:
	cd 11_3D-RISM && ./Run.rism3d

ANTE-MMPBSA:
	cd EstRAL_Files && ./Run.makeparms

clean:
	-cd EstRAL_Files && /bin/rm -f ?.top _AnteMMPBSA_* *.dif
	-cd 01_Generalized_Born && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 02_Poisson_Boltzmann && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 03_Alanine_Scanning && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 04_Per_Residue_Decomp && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 05_Pairwise_Decomp && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 06_NAB_Nmode && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 07_Comprehensive && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA*.dat *.dif *.csv
	-cd 08_Stability && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA*.dat *.dif *.csv
	-cd 09_APBS_Poisson_Boltzmann && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 10_QM_MMGBSA && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
	-cd 11_3D-RISM && /bin/rm -f _MMPBSA_* mmpbsa.* FINAL_*_MMPBSA.dat *.dif *.csv
