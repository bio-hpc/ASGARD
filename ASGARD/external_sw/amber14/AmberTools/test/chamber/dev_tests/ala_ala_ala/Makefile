#Defines specific to the AMBER install
include ../common.mk

#Defines local to the test
TESTNAME=ala_ala_ala
PSF=$(TESTNAME).psf
COR=$(TESTNAME).pdb
CHARMM_PAR=../../../../dat/chamber/par_all22_prot.inp
CHARMM_TOP=../../../../dat/chamber/top_all22_prot.inp
MAX_NUM_OF_PROC=2

SYSTEM= -psf $(PSF) -crd $(COR)
CMDLINE= -param $(CHARMM_PAR) -top $(CHARMM_TOP) $(SYSTEM) -cmap

#COMPARABLES=chamber.out prmtop inpcrd !TODO

serial: prmtop mdinfo
	@$(DACDIF) ./save/mdinfo ./mdinfo
	@$(DACDIF) -r 10E-15 ./save/prmtop ./prmtop
	@$(DACDIF) ./save/inpcrd ./inpcrd
	$(FFCOMPARE) ./save/charmm_gold ./charmm_gold
	@$(MAKE) clean

prmtop: $(PSF) $(COR)
	@$(CHAMBER) $(CMDLINE) > chamber.out
	@$(DACDIF) ./save/chamber.out ./chamber.out


#Serial
mdinfo: prmtop mdin
	@$(MD_ENGINE) -O
	@#$(DACDIF) -r 10E-5 ./save/charmm_gold ./charmm_gold

#parallel
parallel: prmtop mdinfo.mpi
	@$(DACDIF) ./save/mdinfo ./mdinfo.mpi
	@$(MAKE) clean

mdinfo.mpi: prmtop mdin
	MPIRUN_CMD= $(MPIRUN) -n $(MAX_NUM_OF_PROC) $(MD_ENGINE_MPI) -O -inf mdinfo.mpi
	$(MPIRUN_CMD)


clean:
	@rm -f chamber.out
	@rm -f prmtop inpcrd
	@rm -f restrt mdout mdinfo*
	@rm -f charmm_gold



#Local debugging

namd:   $(PSF) $(COR)
	$(NAMD) mdin.namd > namd.out

debug:
	#echo "break _write_prmtop.f:457" > gdb_commands
	echo "break main" > gdb_commands
	echo "run $(CMDLINE)" >> gdb_commands
	$(GDB) $(CHAMBER) -x ./gdb_commands


debug_gradients: prmtop mdin_debug_gradients
	$(MD_ENGINE) -i mdin_debug_gradients -O
	cat mdout

vmd:
	vmd $(COR)

update_saved_files:
	mv prmtop ./save/
	mv inpcrd ./save/
	mv chamber.out ./save/
	mv mdinfo ./save/
	mv ala_ala_ala_forces.dat ./save/charmm_gold
