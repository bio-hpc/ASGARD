#Defines specific to the AMBER install
include ../../common.mk

#Defines local to the test
TESTNAME=ala_ala_ala
PSF=$(TESTNAME).psf
COR=$(TESTNAME).pdb


$(PSF): charmm.out
$(COR): charmm.out


all: charmm.out clean

charmm.out: $(TESTNAME).inp
	@$(CHARMM) < $(TESTNAME).inp > charmm.out
	@cp $(PSF) ..
	@cp $(COR) ..
	@$(MAKE) clean

clean:
	@rm -f *.dat  *.pdb  *.psf  *.xplor  *.out

