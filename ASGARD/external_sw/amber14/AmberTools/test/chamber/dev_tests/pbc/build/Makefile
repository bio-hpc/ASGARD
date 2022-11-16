#Defines specific to the AMBER install
include ../../common.mk




all: heat.test1.res heat.test2.res heat.test3.res
	cp waterbox.psf ..
	cp *.res ..
	cp *.out ..
	
waterbox.psf: generate_waterbox_psf.inp watbox222.crd
	$(CHARMM) < generate_waterbox_psf.inp > generate_waterbox_psf.out


heat.test1.res: waterbox.psf
	$(CHARMM) < test1.inp > test1.out

heat.test2.res: waterbox.psf
	$(CHARMM) < test2.inp > test2.out

heat.test3.res: waterbox.psf
	$(CHARMM) < test3.inp > test3.out


clean:
	rm -f test*.out *.res waterbox.cor  waterbox.psf
