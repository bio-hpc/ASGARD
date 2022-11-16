#Defines specific to the AMBER install
include ../common.mk

#Defines local to the test
TESTNAME=trans-but-2-ene
PSF=$(TESTNAME).psf
COR=$(TESTNAME).coor
CHARMM_PAR=../../../../dat/chamber/par_all22_prot.inp
CHARMM_TOP=../../../../dat/chamber/top_all22_prot.inp
MAX_NUM_OF_PROC=2

SYSTEM= -psf $(PSF) -crd $(COR)
CMDLINE= -param $(CHARMM_PAR) -top $(CHARMM_TOP) $(SYSTEM) -cmap

serial: prmtop mdinfo
	@$(DACDIF) ./save/mdinfo ./mdinfo
	@$(DACDIF) -r 10E-15 ./save/prmtop ./prmtop 
	@$(DACDIF) ./save/inpcrd ./inpcrd
	@$(MAKE) clean

prmtop: $(PSF) $(COR)
	@$(CHAMBER) $(CMDLINE) > chamber.out
	@$(DACDIF) ./save/chamber.out ./chamber.out


#Serial
mdinfo: prmtop mdin
	@$(MD_ENGINE) -O
	@#$(FFCOMPARE) ./save/charmm_gold ./charmm_gold


#parallel
parallel: prmtop mdinfo.mpi
	@$(DACDIF) ./save/mdinfo ./mdinfo.mpi
	@$(MAKE) clean

mdinfo.mpi: prmtop mdin
	MPIRUN_CMD= $(MPIRUN) -n $(MAX_NUM_OF_PROC) $(MD_ENGINE_MPI) -O -inf mdinfo.mpi
	$(MPIRUN_CMD)



#awk script compare
compare: charmm.out mdinfo
	@awk -f ../comp_ene.awk charmm.out mdinfo
	@$(MAKE) clean

compare.parallel: charmm.out mdinfo.mpi
	@awk -f ../comp_ene.awk charmm.out mdinfo.mpi
	@$(MAKE) clean



clean:
	@rm -f chamber.out
	@rm -f prmtop inpcrd
	@rm -f restrt mdout mdinfo
	@rm -f charmm_gold
	@rm -f *.out



vmd:
	vmd $(COR)

