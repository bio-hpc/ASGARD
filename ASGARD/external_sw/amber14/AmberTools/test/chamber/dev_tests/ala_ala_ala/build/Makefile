#Defines specific to the AMBER install
include ../../common.mk

#Defines local to the test
TESTNAME=ala_ala_ala
PSF=$(TESTNAME).psf
PSFX=$(TESTNAME).psf.xplor
COR=$(TESTNAME).pdb


all: charmm.out $(PSF) $(COR)

$(PSF): charmm.out
$(COR): charmm.out

charmm.out: $(TESTNAME).inp
	@$(CHARMM) < $(TESTNAME).inp > charmm.out
	@cp $(PSF) $(PSFX) ..
	@cp $(COR) ..
	@cp charmm_gold ../save/

#	@$(MAKE) clean

clean:
	@rm -f *.dat  *.pdb  *.psf  *.xplor  *.out charmm_gold

