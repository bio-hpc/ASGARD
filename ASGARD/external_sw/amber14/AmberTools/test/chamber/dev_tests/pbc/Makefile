#Defines specific to the AMBER install
include ../common.mk

CMDLINE_TEST1=     -cmap\
             -verbose\
             -top ../../../../dat/chamber/top_all22_prot.inp \
             -param ../../../../dat/chamber/par_all22_prot.inp \
             -psf waterbox.psf \
             -crd test1.res \
             -p prmtop.test1 \
             -inpcrd inpcrd.test1 \

CMDLINE_TEST2=     -cmap\
             -verbose\
             -top ../../../../dat/chamber/top_all22_prot.inp \
             -param ../../../../dat/chamber/par_all22_prot.inp \
             -psf waterbox.psf \
             -crd test2.res \
             -p prmtop.test2 \
             -inpcrd inpcrd.test2 \

CMDLINE_TEST3=     -cmap\
             -verbose\
             -top ../../../../dat/chamber/top_all22_prot.inp \
             -param ../../../../dat/chamber/par_all22_prot.inp \
             -psf waterbox.psf \
             -crd test3.res \
             -p prmtop.test3 \
             -inpcrd inpcrd.test3 \


all: test1 test2 test3
	@$(MAKE) clean

test1: mdout.test1 comp_ene.awk test1.res 
	@echo "=============="
	@echo "Running test 1"
	@echo "=============="
	awk -f comp_ene.awk test1.out mdout.test1 mdinfo.test1
	awk -f comp_ektot.awk test1.out mdinfo.test1
	
test2: mdout.test2 comp_ene.awk test2.res
	@echo "=============="
	@echo "Running test 2"
	@echo "=============="
	awk -f comp_ene.awk test2.out mdout.test2 mdinfo.test2
	awk -f comp_ektot.awk test2.out  mdinfo.test2

test3: mdout.test3 comp_ene.awk test3.res
	@echo "=============="
	@echo "Running test 3"
	@echo "=============="
	awk -f comp_ene.awk test3.out mdout.test3 mdinfo.test3
	awk -f comp_ektot.awk test3.out mdinfo.test3



mdout.test1: prmtop.test1 mdin.test1
	@$(MD_ENGINE)  -O -i mdin.test1 -o mdout.test1 -c inpcrd.test1 -p prmtop.test1 -inf mdinfo.test1
 
mdout.test2: prmtop.test2 mdin.test2
	@$(MD_ENGINE)  -O -i mdin.test2 -o mdout.test2 -c inpcrd.test2 -p prmtop.test2 -inf mdinfo.test2

mdout.test3: prmtop.test3 mdin.test3
	@$(MD_ENGINE)  -O -i mdin.test3 -o mdout.test3 -c inpcrd.test3 -p prmtop.test3 -inf mdinfo.test3



prmtop.test1: test1.res
	@$(CHAMBER) $(CMDLINE_TEST1) > chamber.test1.out

prmtop.test2: test2.res
	@$(CHAMBER) $(CMDLINE_TEST2) > chamber.test2.out

prmtop.test3: test3.res
	@$(CHAMBER) $(CMDLINE_TEST3) > chamber.test3.out




vmd: prmtop mdin
	~/code/VMD/bin/vmd -parm7 prmtop -crdbox mdcrd

clean:
	@rm -f prmtop*  chamber* 
	@rm -f inpcrd* restrt mdout mdinfo* mdcrd
	@rm -f *.dat
	@rm -f mdout*
