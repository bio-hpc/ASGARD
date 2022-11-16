#Defines specific to the AMBER install
include ../../common.mk

#Defines local to the test
TESTNAME=trans-but-2-ene
PSF=$(TESTNAME).psf
COR=$(TESTNAME).coor


$(PSF): charmm.out
$(COR): charmm.out


all: charmm.out clean

charmm.out: $(TESTNAME).inp
	@$(CHARMM) < $(TESTNAME).inp > charmm.out
	@cp $(PSF) ..
	@cp $(COR) ..
	@cp charmm.out ..
	@$(MAKE) clean

clean:
	@rm -f *.dat  *.pdb  *.psf  *.xplor  *.out *.coor

