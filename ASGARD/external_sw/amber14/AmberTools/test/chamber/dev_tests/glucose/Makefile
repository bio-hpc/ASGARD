#Defines specific to the AMBER install
include ../common.mk
#PAR=par_allxx_sugar.inp
#TOP=top_allxx_sugar.rtf
#EXT=no
PAR=par_cbh1_cellulose.prm
TOP=top_cbh1_cellulose.rtf
EXT=yes

serial: prmtop mdinfo 
	@$(DACDIF) ./save/mdinfo ./mdinfo
	@$(DACDIF) -r 10E-15 ./save/prmtop ./prmtop
	@$(DACDIF) ./save/inpcrd ./inpcrd
	@$(DACDIF) ./save/chamber.out ./chamber.out
	@#$(FFCOMPARE) ./save/charmm_gold ./charmm_gold  #To add in when I get clearance on the FFCOMPARE script

#	@$(MAKE) clean

clean:
	@rm -f restrt prmtop mdout mdinfo inpcrd gluc_energy.out chamber.out *.dat charmm_gold glucmin_new.pdb

mdinfo: prmtop 
	$(MD_ENGINE) -O -c inpcrd -p prmtop


prmtop:     #gluc.psf glucmin.pdb 
	@$(CHAMBER) \
         -nocmap \
         -psf gluc.psf \
         -crd glucmin.pdb \
         -param $(PAR) \
         -top   $(TOP) \
         -p prmtop \
         -inpcrd inpcrd > chamber.out

charmm: gluc_energy.out mdinfo
	awk -f ../comp_ene.awk gluc_energy.out mdinfo

gluc.psf: save/gluc.psf
	/bin/cp save/gluc.psf .

gluc_energy.out: 
	$(CHARMM) par=$(PAR) top=$(TOP) ext=$(EXT) < gluc_energy.inp > gluc_energy.out

vmd: prmtop
	vmd -parm7 prmtop -rst7 inpcrd
