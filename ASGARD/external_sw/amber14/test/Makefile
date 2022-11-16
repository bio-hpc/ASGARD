include ../src/config.h

SHELL=/bin/sh

#Cuda Options
#------------
#Precision model to test (options are SPFP, SPDP, DPDP) - default is SPFP
PREC_MODEL=SPFP
#------------

######    Top-level targets, for general use:    ############################

test: clean is_amberhome_defined
	./test_amber_serial.sh

test.parallel: clean is_amberhome_defined
	./test_amber_parallel.sh

test.mic_offload: clean is_amberhome_defined
	./test_amber_mic_offload.sh

test.cuda: clean is_amberhome_defined
	./test_amber_cuda.sh

test.cuda.parallel: clean is_amberhome_defined
	./test_amber_cuda_parallel.sh

######    Intermediate-level targets for AmberTools:    #############

test.serial:  test.serial.MM test.serial.QMMM test.serial.sander.SEBOMD \
        test.serial.emil test.serial.sanderapi

test.serial.sander.MM: test.sander.BASIC.nopar test.sander.BASIC \
	test.sander.GB test.sander.LES test.sander.PIMD.partial test.amoeba \
	test.ncsu.serial test.pbsa test.sander.constph test.sander.RISM \
	test.sander.crgreloc test.serial.sander.AMD test.serial.sander.sgld \
	test.serial.sander.emap test.sander.CHARMM

test.serial.MM: test.serial.sander.MM test.nmode

test.serial.QMMM: \
	test.sander.QMMM.nopar test.sander.QMMM \
	test.sander.DFTB.nopar test.sander.DFTB \
	test.sander.EXTERN test.sander.ABFQMMM

test.serial.sander: test.serial.sander.MM test.serial.QMMM

######    Intermediate-level targets for Amber:    ##################

test.serial.pmemd: test.serial.pmemd.basic test.serial.pmemd.TI \
	test.serial.pmemd.sgld test.pmemd.amoeba test.serial.pmemd.emap \
	test.serial.emil.pmemd finished.serial

test.parallel.pmemd: test.parallel.pmemd.basic test.parallel.pmemd.TI \
	test.pmemd.amoeba.MPI test.parallel.pmemd.sgld  test.parallel.pmemd.emap \
	test.parallel.emil.pmemd test.pmemd.REM finished.parallel

test.mic_offload.pmemd: test.mic_offload.pmemd.basic test.mic_offload.pmemd.TI \
	test.mic_offload.pmemd.sgld test.mic_offload.pmemd.emap \
	test.mic_offload.emil.pmemd test.mic_offload.pmemd.REM finished.mic_offload 

test.serial.cuda: test.pmemd.cuda

test.parallel.cuda: test.pmemd.cuda.MPI

test.parallel.at:  test.parallel.MM test.parallel.QMMM test.parallel.EVB \
	test.parallel.sander.SEBOMD test.parallel.emil

######   More specialized rules; not for general use ################

test.parallel.4proc:
	-cd rem_gb_4rep && ./Run.rem
	-cd PIMD/full_cmd_water/equilib && ./Run.full_cmd
	-cd PIMD/full_cmd_water/start && ./Run.full_cmd
	-cd PIMD/full_cmd_water/restart && ./Run.full_cmd
	-cd PIMD/full_nmpimd_water && ./Run.full_pimd
	-cd PIMD/full_pimd_amoeba && ./Run.full_pimd_amoeba
	-cd PIMD/full_pimd_ln_water && ./Run.full_pimd_ln
	-cd PIMD/full_pimd_nhc_water && ./Run.full_pimd_nhc
	-cd PIMD/full_pimd_ntp_water && ./Run.full_pimd_ntp
	-cd PIMD/full_rpmd_water && ./Run.full_rpmd
	-cd ncsu && ./run-parallel.sh

test.parallel.MM: test.sander.BASIC.MPI test.sander.GB.MPI \
	test.sander.REM test.sander.TI test.sander.LES.MPI \
	test.neb test.ncsu.parallel test.sander.CHARMM.MPI \
	test.sander.PIMD.MPI.partial test.sander.PIMD.MPI.full \
	test.sander.constph.MPI  test.sander.RISM.MPI \
	test.sander.crgreloc.MPI test.parallel.sander.AMD \
	test.parallel.sander.sgld test.parallel.sander.emap  

test.parallel.QMMM: \
	test.sander.QMMM.MPI test.sander.DFTB.MPI test.sander.EXTERN.MPI \
	test.sander.ADQMMM.MPI test.sander.ABFQMMM.MPI

test.parallel.EVB:  test.sander.EVB

######    Individual targets, for specialized use:    ######################

test.sander.BASIC.nopar:
	-cd dmp && ./Run.dmp
	-cd adenine && ./Run.adenine
	-cd variable_14 && ./Run.variable_14_igb1

test.sander.BASIC:
	-cd cytosine && ./Run.cytosine
	-cd nonper && ./Run.nonper
	-cd nonper && ./Run.nonper.belly
	-cd nonper && ./Run.nonper.belly.mask
	-cd nonper && ./Run.nonper.min
	-cd nonper && ./Run.cap
	-cd nonper && ./Run.nonper.nocut
	-cd tip4p && ./Run.tip4p
	-cd tip4p && ./Run.tip4p_mcbar
	-cd tip4p && ./Run.tip4p_nve
	-cd tip5p && ./Run.tip5p
	-cd tip5p && ./Run.tip5p_nve
	-cd 4096wat && ./Run.pure_wat
	-cd 4096wat && ./Run.pure_wat_nmr_temp_reg
	-cd 4096wat && ./Run.vrand
	-cd 4096wat && ./Run.frcdmp
	-cd 4096wat_oct && ./Run.pure_wat_oct
	-cd dhfr && ./Run.dhfr && ./Run.dhfr prmtop_an
	-cd dhfr && ./Run.dhfr.noboxinfo
	-cd dhfr && ./Run.dhfr.noshake
	-cd dhfr && ./Run.dhfr.min
	-cd dhfr && ./Run.dhfr.min_prcg
	-cd dhfr && ./Run.dhfr.lmodxmin
	-cd ff14ipq && ./Run.ff14ipq
	-cd variable_14 && ./Run.variable_14_ntb1
	-cd gact_ips && ./Run.ips
	-cd gact_ips && ./Run.ipsnve
	-cd gact_ips && ./Run.ips_dfft
	-cd gact_ips && ./Run.ips_vacuum
	-cd polarizable_water && ./Run.pol_wat
	-cd ubiquitin && ./Run.ubiquitin
	-cd dna_pol && ./Run.dna_pol
	-cd trx && ./Run.trx
	-cd trx && ./Run.trx.cpln
	-cd rdc && ./Run.dip
	-cd rdc && ./Run.csa
	-cd rdc && ./Run.csa_min
	-cd pcsa && ./Run.dip
	-cd pcsa && ./Run.csa
	-cd pcsa && ./Run.csa2
	-cd pcsa && ./Run.csa_min
	-cd tgtmd/change_target && ./Run.tgtmd
	-cd tgtmd/change_target.rms && ./Run.tgtmd
	-cd tgtmd/change_target.ntr && ./Run.tgtmd
	-cd tgtmd/conserve_ene && ./Run.tgtmd
	-cd tgtmd/minimize && ./Run.tgtmin
	-cd tgtmd/PME && ./Run.tgtPME
	-cd mtmdtest && ./Run.mtmdtest
	-cd trajene && ./Run.trajene
	-cd trajene_box && ./Run.trajene
	-cd trajene_netcdf && ./Run.trajene $(NETCDF)
	-cd alp && ./Run.alp
	-cd umbrella && ./Run.umbrella
	-cd noesy && ./Run.noesy
	-cd jar && ./Run.jarz
	-cd plane_plane_restraint && ./Run.dinuc_plpt
	-cd plane_plane_restraint && ./Run.dinuc_pln
	-cd iwrap2 && ./Run.iwrap2
	-cd idecomp1 && ./run_idecomp1.sh
	-cd idecomp4 && ./run_idecomp4.sh
	-cd bintraj && ./Run.bintraj $(NETCDF)
	-cd netcdf && $(MAKE) -k test OPT=$(NETCDF)
#    Constant surface tension (csurften > 0)
	-cd csurften && ./Run.csurften_z-dir
#    Semi-isotropic pressure coupling (csurften > 0 AND ntp = 3)
	-cd csurften && ./Run.csurften_z-dir_npt_3
#    Solvation Energy via LIE theory
	-cd linear_interaction_energy && ./Run.lie
#    NMR Opt Test Cases
	-cd nmropt && $(MAKE) all
#    LJ 12-6-4 potential test
	-cd lj_12-6-4 && ./Run.12-6-4

test.sander.BASIC.MPI:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; $(MAKE) -k test.sander.BASIC
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd 4096wat && ./Run.column_fft
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd jar_multi && ./Run.jarz
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd ti_eth2meth_gas && ./Run.test1
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd ti_eth2meth_gas && ./Run.test2
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd ti_ggcc && ./Run.test1
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd ti_ggcc && ./Run.test2

test.sander.GB:
	-cd gb_rna && ./Run.gbrna
	-cd gb_rna && ./Run.gbrna.min
	-cd gb_rna && ./Run.gbrna.ln
	-cd gb_rna && ./Run.gbrna.ips
	-cd gbsa_xfin && ./Run.gbsa
	-cd circ_dna && ./Run.circdna
	-cd gb2_trx && ./Run.trxox
	-cd gb7_trx && ./Run.trxox && ./Run.trxox prmtop_an
	-cd gb7_trx && ./Run.trxox_md
	-cd gb8_trx && ./Run.trxox
	-cd gb8_trx && ./Run.trxox_md && ./Run.trxox_md prmtop_an
	-cd gbneck2nu_test && ./Run.1hji_md
	-cd alpb_trx && ./Run.trxox
	-cd gb1_cox2 && ./Run.cox2

test.serial.sander.SEBOMD:
	-cd sebomd/AM1/dimethylether       && ./min.csh && ./md.csh
	-cd sebomd/AM1/methionine          && ./min.csh && ./md.csh
	-cd sebomd/AM1/n-methyl-acetamide  && ./min.csh && ./md.csh
	-cd sebomd/AM1/nucleoside          && ./min.csh && ./md.csh
	-cd sebomd/AM1/water1              && ./min.csh && ./md.csh
	-cd sebomd/AM1D/dimethylether      && ./min.csh && ./md.csh
	-cd sebomd/AM1D/methionine         && ./min.csh && ./md.csh
	-cd sebomd/AM1D/n-methyl-acetamide && ./min.csh && ./md.csh
	-cd sebomd/AM1D/nucleoside         && ./min.csh && ./md.csh
	-cd sebomd/AM1D/water1             && ./min.csh && ./md.csh
	-cd sebomd/DC/water32              && ./min.csh && ./md.csh
	-cd sebomd/Ewald-CM1/Water2        && ./md.csh
	-cd sebomd/Ewald-ConstantCharge/Water2   && ./md.csh
	-cd sebomd/Ewald-Mulliken/Water2   && ./md.csh
	-cd sebomd/MinimumImage/Water2     && ./md.csh
	-cd sebomd/method3/alaninedipeptide  && ./min.csh && ./md.csh
	-cd sebomd/method3/ethanol_and_water && ./min.csh && ./md.csh
	-cd sebomd/MNDO/dimethylether      && ./min.csh && ./md.csh
	-cd sebomd/MNDO/methionine         && ./min.csh && ./md.csh
	-cd sebomd/MNDO/n-methyl-acetamide && ./min.csh && ./md.csh
	-cd sebomd/MNDO/nucleoside         && ./min.csh && ./md.csh
	-cd sebomd/MNDO/water1             && ./min.csh && ./md.csh
	-cd sebomd/RM1/dimethylether       && ./min.csh && ./md.csh
	-cd sebomd/RM1/methionine          && ./min.csh && ./md.csh
	-cd sebomd/RM1/n-methyl-acetamide  && ./min.csh && ./md.csh
	-cd sebomd/RM1/nucleoside          && ./min.csh && ./md.csh
	-cd sebomd/RM1/water1              && ./min.csh && ./md.csh
	-cd sebomd/PIF2/dimethylether      && ./min.csh && ./md.csh
	-cd sebomd/PIF2/water1             && ./min.csh && ./md.csh
	-cd sebomd/PIF2/doublewater        && ./min.csh && ./md.csh
	-cd sebomd/PIF2/ethanol_and_water  && ./min.csh && ./md.csh
	-cd sebomd/PIF3/water1             && ./min.csh && ./md.csh
	-cd sebomd/PIF3/doublewater        && ./min.csh && ./md.csh
	-cd sebomd/PIF3/ethanol_and_water  && ./min.csh && ./md.csh
	-cd sebomd/PM3/alaninedipeptide    && ./min.csh && ./md.csh
	-cd sebomd/PM3/dimethylether       && ./min.csh && ./md.csh
	-cd sebomd/PM3/methionine          && ./min.csh && ./md.csh
	-cd sebomd/PM3/n-methyl-acetamide  && ./min.csh && ./md.csh
	-cd sebomd/PM3/nucleoside          && ./min.csh && ./md.csh
	-cd sebomd/PM3/water1              && ./min.csh && ./md.csh

test.parallel.sander.SEBOMD:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd sebomd/DC/water32      && ./min.csh && ./md.csh

test.sander.constph:
	cd cnstph/implicit && ./Run.cnstph
#	cd cnstph/explicit && ./Run.cnstph

test.sander.constph.MPI:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph/implicit && ./Run.cnstph
#	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph/explicit && ./Run.cnstph
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph_remd/TempRem && ./Run.cnstph_remd
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph_remd/pHREM && ./Run.pHremd
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph_remd/Explicit_pHREM && ./Run.pHremd

test.sander.GB.MPI:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; $(MAKE) -k test.sander.GB

test.sander.CHARMM:
	-cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.min
	-cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.md
	-cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.min
	-cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.md
	-cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.min
	-cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.md
	-cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.min
	-cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.md

test.sander.CHARMM.MPI:
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; $(MAKE) -k test.sander.CHARMM

#  following are parallel only:
test.sander.TI:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd sodium && ./Run.sodium
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd pheMTI && ./Run.0
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd pheMTI && ./Run.1
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd pheMTI && ./Run.lambda0
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd pheMTI && ./Run.lambda1
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd campTI && ./Run.campTI
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd softcore && ./Run_sc
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd ti_decomp; ./Run.ti_decomp
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd ti_decomp; ./Run.ti_decomp_sc

#  following are parallel only:
test.sander.EVB:
	-cd evb/malon_bond_umb_GAFF && ./Run.evb
	-cd evb/malon_bond_umb_GAFF_MORSIFY && ./Run.evb
	-cd evb/malon_bond_umb_GAFF_MORSIFY_XVDW && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_nmpimd_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_pimd_ld_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_pimd_nhc_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_EXP-XCH && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_EXP-XCH_pimd_ld_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_EXP-XCH_pimd_nhc_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_EXP-XCH_nmpimd_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_GAUSS-XCH && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_GAUSS-XCH_pimd_ld_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_GAUSS-XCH_pimd_nhc_full && ./Run.evb
	-cd evb/malon_dbonds_umb_GAFF_MORSIFY_XVDW_GAUSS-XCH_nmpimd_full && ./Run.evb
	-cd evb/malon_emap_GAFF_MORSIFY && ./Run.evb
	-cd evb/malon_emap_GAFF_MORSIFY_XVDW && ./Run.evb
	-cd evb/malon_bond_umb_dg_UFF_3DG && ./Run.evb
	-cd evb/malon_bond_umb_dg_UFF_3DG_nmpimd_full && ./Run.evb
	-cd evb/malon_bond_umb_dg_UFF_3DG_pimd_ld_full && ./Run.evb
	-cd evb/malon_bond_umb_dg_UFF_3DG_pimd_nhc_full && ./Run.evb
	-cd evb/malon_bond_umb_dg_UFF_3DG_pimd_nhc_full_TI-MASS && ./Run.evb
	-cd evb/malon_dbond_umb_dg_UFF_3DG_pimd_nhc_full_TI-MASS && ./Run.evb
	-cd evb/malon_dbonds_umb_dg_UFF_3DG_qi_full_2D-PMF && ./Run.evb
	-cd evb/malon_dbonds_umb_dg_UFF_3DG_qi_full_corrF && ./Run.evb
	-cd evb/poh_dbonds_umb_dg_UFF_3DG && ./Run.evb
	-cd evb/poh_dbonds_umb_dg_UFF_9DG && ./Run.evb
	-cd evb/poh_dbonds_umb_dg_UFF_9DG_pimd_ld_full && ./Run.evb
	-cd evb/poh_dbonds_umb_dg_UFF_9DG_pimd_nhc_full && ./Run.evb
	-cd evb/poh_dbonds_umb_dg_UFF_9DG_nmpimd_full && ./Run.evb
	-cd evb/poh_dbonds_umb_dg_UFF_9DG_nmpimd_full_TST-freqf && ./Run.evb

test.sander.EVBDEBUG:
	-cd evb_debug/hcn_dg_cart && ./Run.evb
	-cd evb_debug/hcn_dg_cart_irc && ./Run.evb
	-cd evb_debug/hcn_dg_cart_irc_dgelsd && ./Run.evb
	-cd evb_debug/hcn_dg_cart_irc_dgelsy && ./Run.evb
	-cd evb_debug/hcn_dg_int_irc_grad && ./Run.evb
	-cd evb_debug/poh_EVB_CART_DIIS && ./Run.evb
	-cd evb_debug/poh_EVB_INT_DIIS && ./Run.evb
	-cd evb_debug/poh_gfchk_GMRES && ./Run.evb
	-cd evb_debug/poh_gfchk_GMRES_K3_R32_UFF_dyn && ./Run.evb
	-cd evb_debug/poh_gfchk_GMRES_K5_R32_UFF_dyn && ./Run.evb
	-cd evb_debug/poh_gfchk_GMRES_K7_R32_UFF_dyn && ./Run.evb
	-cd evb_debug/poh_gfchk_GMRES_K9_R32_UFF_dyn && ./Run.evb

test.sander.LES:
	-cd LES_noPME && ./Run.LESmd
	-cd LES_noPME && ./Run.LESmd.rdiel
	-cd LES && ./Run.PME_LES
	-cd LES_CUT && ./Run.LES
	-cd LES_TEMP && ./Run.2temp
	-cd LES_GB && ./Run.LES

test.sander.LES.MPI:
	export TESTsanderLES=$(AMBERHOME)/bin/sander.LES.MPI; $(MAKE) -k test.sander.LES

test.sander.QMMM.nopar:
	-cd qmmm2/input_validation && ./Run.invalid_atom_ids
	-cd qmmm2/input_validation && ./Run.non_unique_qm_atoms
	-cd qmmm2/pure_QM_MD && ./Run.nma
	-cd qmmm2/pure_QM_MD && ./Run.nma_md
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_0
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_2
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_3
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_4
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_5
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_6
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_diag_routine_7
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_qmgb0
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_qmgb1
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_qmgb2
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_qmgb2_saltcon
	-cd qmmm2/pure_QM_MD && ./Run.nma_md_qmgb3
	-cd qmmm2/pure_QM_MD && ./Run.C6H4FCL_md
	-cd qmmm2/pure_QM_MD_shake && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_tight_p && ./Run.nma
	-cd qmmm2/pure_QM_MD_tight_p && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_PDDGPM3 && ./Run.nma
	-cd qmmm2/pure_QM_MD_PDDGPM3 && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_PDDGPM3 && ./Run.C6H4FCL_md
	-cd qmmm2/pure_QM_MD_PDDGPM3_08 && ./Run.nma
	-cd qmmm2/pure_QM_MD_PDDGPM3_08 && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_PDDGMNDO && ./Run.nma
	-cd qmmm2/pure_QM_MD_PDDGMNDO && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_PDDGMNDO && ./Run.C6H4FCL_md
	-cd qmmm2/pure_QM_MD_AM1 && ./Run.nma
	-cd qmmm2/pure_QM_MD_AM1 && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_AM1 && ./Run.C6H4FCL_md
	-cd qmmm2/pure_QM_MD_AM1_ADDMM && ./Run.nma
	-cd qmmm2/pure_QM_MD_AM1_ADDMM && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_RM1 && ./Run.nma
	-cd qmmm2/pure_QM_MD_RM1 && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_RM1 && ./Run.C6H4FCL_md
	-cd qmmm2/pure_QM_MD_MNDO && ./Run.nma
	-cd qmmm2/pure_QM_MD_MNDO && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_MNDO && ./Run.C6H4FCL_md
	-cd qmmm2/pure_QM_MD_PM6 && ./Run.nma
	-cd qmmm2/pure_QM_MD_PM6 && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_NUMER_DERIV && ./Run.nma
	-cd qmmm2/pure_QM_MD_NUMER_DERIV && ./Run.nma_md
	-cd qmmm2/pure_QM_MD_NUMER_DERIV && ./Run.nma_md_qmgb0
	-cd qmmm2/pure_QM_MD_NUMER_DERIV && ./Run.nma_md_qmgb2
	-cd qmmm2/pure_QM_MD_NUMER_DERIV && ./Run.nma_md_qmgb2_saltcon
	-cd qmmm2/1-2-ethanediol_PM3CARB1 && ./Run.PM3CARB1_md
	-cd qmmm2/lysine_PM3 && ./Run.lysine
	-cd qmmm2/lysine_PM3 && ./Run.lysine_md
	-cd qmmm2/lysine_PM3 && ./Run.lysine_md_qmmm_int0
	-cd qmmm2/lysine_PM3 && ./Run.lysine_md_qmmm_int2
	-cd qmmm2/lysine_PM3 && ./Run.lysine_md_lnk_dis_neg
	-cd qmmm2/lysine_PM3 && ./Run.lysine_md_lnk_method2
	-cd qmmm2/lysine_PM3_qmgb0 && ./Run.lysine
	-cd qmmm2/lysine_PM3_qmgb0 && ./Run.lysine_md
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md_igb2
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md_igb5
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md_igb7
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md_adjust_q1
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md_adjust_q2
	-cd qmmm2/lysine_PM3_qmgb2 && ./Run.lysine_md_lnk_dis_neg
	-cd qmmm2/lysine_PM3_qmgb2_alpb && ./Run.lysine_md
	-cd qmmm2/lysine_PM3_calc_mulliken && ./Run.lysine
	-cd qmmm2/lysine_PM3_calc_mulliken && ./Run.lysine_md
	-cd qmmm2/lysine_PM3_calc_dipole && ./Run.lysine
	-cd qmmm2/lysine_PM3_calc_dipole && ./Run.lysine_md
	-cd qmmm2/lysine_PM3_calc_dipole && ./Run.lysine_qmmm
	-cd qmmm2/lysine_AM1 && ./Run.lysine
	-cd qmmm2/lysine_AM1 && ./Run.lysine_md
	-cd qmmm2/lysine_MNDO && ./Run.lysine
	-cd qmmm2/lysine_MNDO && ./Run.lysine_md
	-cd qmmm2/lysine_PM6 && ./Run.lysine
	-cd qmmm2/lysine_PM6 && ./Run.lysine_md
	-cd qmmm2/meoh-nh3_PM6-DH+ && ./Run.meoh-nh3
	-cd qmmm2/meoh-nh3_PM6-DH+ && ./Run.meoh-nh3_md
	-cd qmmm2/PM3-MMX && ./Run.water_dimer.x
	-cd qmmm2/PM3-MMX && ./Run.ammonia_dimer.x
	-cd qmmm2/PM3-MMX && ./Run.formic_acid_dimer.x
	-cd qmmm2/PM3-MMX2 && ./Run.water_dimer.x
	-cd qmmm2/PM3-MMX2 && ./Run.ammonia_dimer.x
	-cd qmmm2/PM3-MMX2 && ./Run.formic_acid_dimer.x
	-cd qmmm2/PM3-MMX2 && ./Run.methanethiol.x
	-cd qmmm2/znme2/ && ./Run.min
	-cd qmmm2/znme2/ && ./Run.md
	-cd qmmm2/cysdip/ && ./Run.mndod.min
	-cd qmmm2/cysdip/ && ./Run.mndod.md
	-cd qmmm2/cysdip/ && ./Run.mndod-mm.min
	-cd qmmm2/cysdip/ && ./Run.mndod-mm.md

test.sander.QMMM:
	-cd qmmm2/xcrd_build_test/ && ./Run.oct_nma_imaged
	-cd qmmm2/xcrd_build_test/ && ./Run.oct_nma_noimage
	-cd qmmm2/xcrd_build_test/ && ./Run.ortho_qmewald0
	-cd qmmm2/xcrd_build_test/ && ./Run.truncoct_qmewald0
	-cd qmmm2/crambin_2 && ./Run.crambin
	-cd qmmm2/crambin_2 && ./Run.crambin_md
	-cd qmmm2/crambin_2 && ./Run.crambin_md_calc_mulliken
	-cd qmmm2/crambin_2 && ./Run.crambin_md_qmgb2
	-cd qmmm2/crambin_2 && ./Run.crambin_md_qmgb2_alpb
	-cd qmmm2/crambin_2 && ./Run.crambin_md_qmgb2_aq1
	-cd qmmm2/crambin_2 && ./Run.crambin_md_qmgb2_aq2
	-cd qmmm2/crambin_2 && ./Run.crambin_md_qmgb2_saltcon
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_internal
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_dspev
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_dspevd
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_dspevx
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_dsyev
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_dsyevd
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_dsyevr
	-cd qmmm2/1NLN_test_diagonalizers && ./Run.1NLN_auto
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_min
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_aq1
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_aq2
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_mmcorr
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_qmewald2
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb2
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_aq2_NUMER_DERIV
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_lnk_dis_neg
	-cd qmmm2/1NLN_periodic_lnk_atoms && ./Run.1NLN_MD_ntb1_lnk_method2
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.notimaged
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.notimaged_md
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.notimaged_md_pme
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.notimaged_md_pme_qmewald
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.notimaged_md_pme_qmewald_lowmem
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.imaged
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.imaged_md
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.img_center
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.img_center_md
	-cd qmmm2/MG_QM_water_MM_AM1_periodic && ./Run.img_center_bigqmcut_md
	-cd qmmm2/crambin && ./Run.crambin
	-cd qmmm2/crambin && ./Run.crambin_md_hot_start
	-cd qmmm2/crambin && ./Run.crambin_md
	-cd qmmm2/tip3p_cap && ./Run.tip3p_cap
	-cd qmmm2/tip3p_cap && ./Run.tip3p_cap_md
	-cd qmmm2/tip3p_cap && ./Run.tip3p_cap_md.shake
	-cd qmmm2/tip3p_cap && ./Run.tip3p_cap_md.shake_qmshakeoff
	-cd qmmm2/MG_QM_water_MM_AM1_qmcut && ./Run.MG_QM_water_MM
	-cd qmmm2/MG_QM_water_MM_AM1_qmcut && ./Run.MG_QM_water_MM_MD
	-cd qmmm2/MG_QM_water_MM_AM1_qmcut_lowmem && ./Run.MG_QM_water_MM_MD
	-cd qmmm2/MG_QM_water_MM_AM1_qmcut_ADDMM && ./Run.MG_QM_water_MM
	-cd qmmm2/MG_QM_water_MM_AM1_qmcut_ADDMM && ./Run.MG_QM_water_MM_MD
	-cd qmmm2/variable_solvent && ./Run.mg_periodic_vsolv_wat
	-cd qmmm2/variable_solvent && ./Run.1NLN_MD_ntb1_vsolv
	-cd qmmm2/PM3-MMX && ./Run.crambin_md_pm3mmx
	-cd qmmm2/PM3-MMX2 && ./Run.crambin_md_pm3mmx2
	-cd qmmm2/crambin_PM6 && ./Run.crambin
	-cd qmmm2/crambin_PM6 && ./Run.crambin_md
	-cd qmmm2/PM3-MAIS && ./Run.zundel
	-cd qmmm2/PM3-MAIS && ./Run.hcl_wat
	-cd qmmm2/MechEm_nma-spcfwbox && ./Run.MechEm_PDDGPM3
	-cd qmmm2/qmmm_switching_func && ./Run.nma_md_switching

test.sander.QMMM.MPI:
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; $(MAKE) -k test.sander.QMMM

test.sander.DFTB.nopar:
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma_md
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma_md_qmgb0
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma_md_qmgb1
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma_md_qmgb2
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma_md_qmgb2_saltcon
	-cd qmmm_DFTB/pure_QM_MD_DFTB && ./Run.nma_disp
	-cd qmmm_DFTB/pure_QM_MD_shake_DFTB && ./Run.nma_md
	-cd qmmm_DFTB/pure_QM_MD_tight_p_DFTB && ./Run.nma
	-cd qmmm_DFTB/pure_QM_MD_tight_p_DFTB && ./Run.nma_md
	-cd qmmm_DFTB/lysine_DFTB && ./Run.lysine
	-cd qmmm_DFTB/lysine_DFTB && ./Run.lysine_md
	-cd qmmm_DFTB/lysine_DFTB && ./Run.lysine_md_lnk_dis_neg
	-cd qmmm_DFTB/lysine_DFTB && ./Run.lysine_md_lnk_method2
	-cd qmmm_DFTB/lysine_DFTB && ./Run.lysine_md_qmmm_int0
	-cd qmmm_DFTB/lysine_DFTB_qmgb0 && ./Run.lysine
	-cd qmmm_DFTB/lysine_DFTB_qmgb0 && ./Run.lysine_md
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine_md
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine_md_igb2
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine_md_igb5
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine_md_igb7
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine_md_adjust_q1
	-cd qmmm_DFTB/lysine_DFTB_qmgb2 && ./Run.lysine_md_adjust_q2
	-cd qmmm_DFTB/lysine_DFTB_qmgb2_alpb && ./Run.lysine_md
	-cd qmmm_DFTB/lysine_DFTB_calc_mulliken && ./Run.lysine
	-cd qmmm_DFTB/lysine_DFTB_calc_mulliken && ./Run.lysine_md
	-cd qmmm_DFTB/3rd-order/ala8  && ./Run.ala8_PA
	-cd qmmm_DFTB/3rd-order/ala8  && ./Run.ala8_PR
	-cd qmmm_DFTB/3rd-order/water && ./Run.water
	-cd qmmm_DFTB/3rd-order/water && ./Run.water_read_mdin
	-cd qmmm_DFTB/3rd-order/water && ./Run.water_read_file

test.sander.DFTB:
	-cd qmmm_DFTB/crambin_DFTB && ./Run.crambin
	-cd qmmm_DFTB/crambin_DFTB && ./Run.crambin_md_hot_start
	-cd qmmm_DFTB/crambin_DFTB && ./Run.crambin_md
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md_calc_mulliken
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md_qmgb2
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md_qmgb2_alpb
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md_qmgb2_aq1
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md_qmgb2_aq2
	-cd qmmm_DFTB/crambin_2_DFTB && ./Run.crambin_md_qmgb2_saltcon
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_ewald_ntb1
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_ewald_ntb1_link_atoms
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_ewald_ntb1_qmewald2
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_ewald_ntb2
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb1
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb2
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb2_link_atoms
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb2_link_atoms_aq1
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb2_link_atoms_aq2
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb2_qmewald2
	-cd qmmm_DFTB/aladip_tip3p_ewaldpme && ./Run.aladip_pme_ntb2_telec
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_qmcut && ./Run.MG_QM_water_MM
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_qmcut && ./Run.MG_QM_water_MM_MD
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.notimaged
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.notimaged_md
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.notimaged_md_pme
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.notimaged_md_pme_qmewald
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.imaged
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.imaged_md
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.img_center
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.img_center_md
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_periodic && ./Run.img_center_bigqmcut_md
	-cd qmmm_DFTB/MG_QM_water_MM_DFTB_qmcut_lowmem && ./Run.MG_QM_water_MM_MD
	-cd qmmm_DFTB/tip3p_cap_DFTB && ./Run.tip3p_cap
	-cd qmmm_DFTB/tip3p_cap_DFTB && ./Run.tip3p_cap_md
	-cd qmmm_DFTB/tip3p_cap_DFTB && ./Run.tip3p_cap_md.shake
	-cd qmmm_DFTB/tip3p_cap_DFTB && ./Run.tip3p_cap_md.shake_qmshakeoff
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_min
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1_lnk_dis_neg
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1_lnk_method2
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1_aq1
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1_aq1_disp
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1_aq2
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb1_qmewald2
	-cd qmmm_DFTB/1NLN_periodic_lnk_atoms_DFTB && ./Run.1NLN_MD_ntb2
	-cd qmmm_DFTB/3rd-order/ala8 && ./Run.ala8_PA
	-cd qmmm_DFTB/3rd-order/ala8 && ./Run.ala8_PR
	-cd qmmm2/MechEm_nma-spcfwbox && ./Run.MechEm_DFTB

test.sander.DFTB.MPI:
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; $(MAKE) -k test.sander.DFTB

# following are parallel only
test.sander.ADQMMM.MPI:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI && cd qmmm2/adqmmm_na+28water && ./Run.adqmmm
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI && cd qmmm2/adqmmm_na+28water && ./Run.adqmmm-calc_wbk2
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI && cd qmmm2/adqmmm_na+28water && ./Run.adqmmm-fixedR
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI && cd qmmm2/adqmmm_h2o-box && ./Run.adqmmm
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI && cd qmmm2/adqmmm_h2o-box && ./Run.adqmmm-fixedR-calc_wbk2

test.sander.ABFQMMM:
	-cd abfqmmm/abfqmmm_water_sp && ./Run.abfqmmm_water_sp
	-cd abfqmmm/abfqmmm_water_md && ./Run.abfqmmm_water_md
	-cd abfqmmm/abfqmmm_dmpoh_md && ./Run.abfqmmm_dmpoh_md
	-cd abfqmmm/abfqmmm_lysozyme_md && ./Run.abfqmmm_lysozyme_md
	-cd abfqmmm/adqmmm_water_md && ./Run.adqmmm_water_md
	-cd abfqmmm/hotspot_water_md && ./Run.hotspot_water_md

test.sander.ABFQMMM.MPI:
	 export TESTsander=$(AMBERHOME)/bin/sander.MPI; $(MAKE) -k test.sander.ABFQMMM

test.sander.EXTERN:
	-cd qmmm_EXTERN/pure_QM_MD_ADF && ./Run.h2o.dftb
	-cd qmmm_EXTERN/pure_QM_MD_ADF && ./Run.h2o.dftb_md
	-cd qmmm_EXTERN/pure_QM_MD_ADF && ./Run.h2o.lda-sz
	-cd qmmm_EXTERN/pure_QM_MD_ADF && ./Run.h2o.lda-sz_md
	-cd qmmm_EXTERN/pure_QM_MD_ADF && ./Run.h2o.template.lda-sz
	-cd qmmm_EXTERN/pure_QM_MD_GAMESS && ./Run.h2o.bp86_6-31g
	-cd qmmm_EXTERN/pure_QM_MD_GAMESS && ./Run.h2o.bp86_6-31g_md
	-cd qmmm_EXTERN/pure_QM_MD_GAMESS && ./Run.h2o.template.bp86_6-31g
	-cd qmmm_EXTERN/pure_QM_MD_TeraChem && ./Run.h2o.blyp_6-31g
	-cd qmmm_EXTERN/pure_QM_MD_TeraChem && ./Run.h2o.blyp_6-31g_md
	-cd qmmm_EXTERN/pure_QM_MD_TeraChem && ./Run.h2o.template.blyp_6-31g
	-cd qmmm_EXTERN/pure_QM_MD_GAUSSIAN && ./Run.h2o.blyp_sto-3g
	-cd qmmm_EXTERN/pure_QM_MD_GAUSSIAN && ./Run.h2o.blyp_sto-3g_md
	-cd qmmm_EXTERN/pure_QM_MD_GAUSSIAN && ./Run.h2o.template.blyp_sto-3g
	-cd qmmm_EXTERN/pure_QM_MD_Orca && ./Run.h2o.blyp_6-31g
	-cd qmmm_EXTERN/pure_QM_MD_Orca && ./Run.h2o.blyp_6-31g_md
	-cd qmmm_EXTERN/pure_QM_MD_Orca && ./Run.h2o.template.blyp_6-31g
	-cd qmmm_EXTERN/pure_QM_MD_NWChem && ./Run.h2o.blyp_sto-3g
	-cd qmmm_EXTERN/pure_QM_MD_NWChem && ./Run.h2o.blyp_sto-3g_md
	-cd qmmm_EXTERN/pure_QM_MD_Q-Chem && ./Run.h2o.hf_sto-3g
	-cd qmmm_EXTERN/pure_QM_MD_Q-Chem && ./Run.h2o.template.blyp_sto-3g
	-cd qmmm_EXTERN/pure_QM_MD_Q-Chem && ./Run.h2o.blyp_6-31g_md
	-cd qmmm_EXTERN/QMMM_MD_ADF && ./Run.aladip.blyp_dz.MechEm
	-cd qmmm_EXTERN/QMMM_MD_ADF && ./Run.nma-spcfw-15.blyp_dz.MechEm
	-cd qmmm_EXTERN/QMMM_MD_Orca && ./Run.aladip.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_Orca && ./Run.aladip.hf_sto-3g.MechEm
	-cd qmmm_EXTERN/QMMM_MD_Orca && ./Run.nma-spcfw-15.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_Orca && ./Run.nma-spcfw-15.hf_sto-3g.MechEm
	-cd qmmm_EXTERN/QMMM_MD_TeraChem && ./Run.aladip.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_TeraChem && ./Run.aladip.hf_sto-3g.MechEm
	-cd qmmm_EXTERN/QMMM_MD_TeraChem && ./Run.nma-spcfw-15.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_TeraChem && ./Run.nma-spcfw-15.hf_sto-3g.MechEm
	-cd qmmm_EXTERN/QMMM_MD_GAUSSIAN && ./Run.aladip.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_GAUSSIAN && ./Run.nma-spcfw-15.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_NWChem && ./Run.aladip.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_Q-Chem && ./Run.aladip.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_Q-Chem && ./Run.aladip.hf_sto-3g.MechEm
	-cd qmmm_EXTERN/QMMM_MD_Q-Chem && ./Run.nma-spcfw-15.hf_sto-3g
	-cd qmmm_EXTERN/QMMM_MD_Q-Chem && ./Run.nma-spcfw-15.hf_sto-3g.MechEm
	-cd qmmm_EXTERN/QMMM_Lio/water && ./Run.water
	-cd qmmm_EXTERN/QMMM_Lio/water && ./Run.td.water

test.sander.EXTERN.MPI:
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_ADF_PIMD && ./Run.h2o.blyp_sz.pimd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_ADF_REMD && ./Run.h2o.blyp_sz.remd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_GAMESS_PIMD && ./Run.h2o.hf_sto-3g.pimd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_GAMESS_REMD && ./Run.h2o.hf_sto-3g.remd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_GAUSSIAN_PIMD && ./Run.h2o.hf_sto-3g.pimd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_GAUSSIAN_REMD && ./Run.h2o.hf_sto-3g.remd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_Q-Chem_PIMD && ./Run.h2o.hf_sto-3g.pimd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_TeraChem_PIMD && ./Run.h2o.hf_sto-3g.pimd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_TeraChem_PIMD && ./Run.h2o.hf_sto-3g.pimd mpi
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_TeraChem_REMD && ./Run.h2o.hf_sto-3g.remd
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/pure_QM_MD_TeraChem_REMD && ./Run.h2o.hf_sto-3g.remd mpi
	export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd qmmm_EXTERN/adqmmm_na+28water_TeraChem && ./Run.adqmmm

test.sander.PIMD.partial:
	-cd PIMD/part_pimd_water && ./Run.pimd
	-cd PIMD/part_nmpimd_water && ./Run.nmpimd
	-cd PIMD/part_pimd_helium && ./Run.pimd
	-cd PIMD/part_nmpimd_helium && ./Run.nmpimd
	-cd PIMD/part_nmpimd_ntp && ./Run.nmpimd
	-cd PIMD/part_pimd_spcfw && ./Run.pimd
	-cd PIMD/part_cmd_water/equilib && ./Run.cmdyn
	-cd PIMD/part_cmd_water/start && ./Run.cmdyn
	-cd PIMD/part_cmd_water/restart && ./Run.cmdyn
	-cd PIMD/part_rpmd_water && ./Run.rpmd
	-cd ti_mass/pent_LES_PIMD && ./Run.pentadiene

test.sander.PIMD.MPI.partial:
	export TESTsander=$(AMBERHOME)/bin/sander.LES.MPI; $(MAKE) -k test.sander.PIMD.partial

test.sander.PIMD.MPI.full:
	cd PIMD/full_cmd_water/equilib && ./Run.full_cmd
	cd PIMD/full_cmd_water/start && ./Run.full_cmd
	cd PIMD/full_cmd_water/restart && ./Run.full_cmd
	cd PIMD/full_nmpimd_water && ./Run.full_pimd
	cd PIMD/full_pimd_amoeba && ./Run.full_pimd_amoeba
	cd PIMD/full_pimd_ln_water && ./Run.full_pimd_ln
	cd PIMD/full_pimd_nhc_water && ./Run.full_pimd_nhc
	cd PIMD/full_pimd_ntp_water && ./Run.full_pimd_ntp
	cd PIMD/full_rpmd_water && ./Run.full_rpmd
	cd ti_mass/pentadiene && ./Run.pentadiene

test.sander.REM:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd rem_gb_2rep && ./Run.rem
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd rem_gb_4rep && ./Run.rem
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd rxsgld_4rep && ./Run.rxsgld
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd rem_wat && ./Run.rem
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd rem_hybrid && ./Run.rem
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd h_rem && ./Run.rem
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd multid_remd && ./Run.multirem

test.sander.allREM: test.sander.REM 
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph_remd/TempRem && ./Run.cnstph_remd
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph_remd/pHREM && ./Run.pHremd
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd cnstph_remd/Explicit_pHREM && ./Run.pHremd

test.sander.OIN:
	-export TESTsander=$(AMBERHOME)/bin/sander; cd sander_OIN && ./Run.sander.OIN

test.sander.OIN.MPI:
	-export TESTsander=$(AMBERHOME)/bin/sander.MPI; cd sander_OIN_MPI && ./Run.sander.OIN.MPI

test.pmemd.REM:
	-export TESTsander='../../bin/pmemd.MPI' && cd rem_gb_2rep && ./Run.rem
	-export TESTsander='../../bin/pmemd.MPI' && cd rem_wat && ./Run.rem
	-export TESTsander='../../bin/pmemd.MPI' && cd rem_gb_4rep && ./Run.rem
	-export TESTsander='../../bin/pmemd.MPI' && cd h_rem && ./Run.rem
	-export TESTsander='../../bin/pmemd.MPI' && cd multid_remd && ./Run.multirem
	-export TESTsander='../../bin/pmemd.MPI'; cd rxsgld_4rep && ./Run.rxsgld
	-export TESTsander='../../../bin/pmemd.MPI'; cd cnstph_remd/pHREM && ./Run.pHremd
	-export TESTsander='../../../bin/pmemd.MPI'; cd cnstph_remd/TempRem && ./Run.cnstph_remd
	-export TESTsander='../../../bin/pmemd.MPI'; cd cnstph_remd/Explicit_pHREM && ./Run.pHremd

test.mic_offload.pmemd.REM:
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI" && cd rem_gb_2rep && ./Run.rem
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI" && cd rem_wat && ./Run.rem
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI" && cd rem_gb_4rep && ./Run.rem
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI" && cd h_rem && ./Run.rem
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI" && cd multid_remd && ./Run.multirem
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd rxsgld_4rep && ./Run.rxsgld
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd cnstph_remd/pHREM && ./Run.pHremd
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd cnstph_remd/TempRem && ./Run.cnstph_remd
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd cnstph_remd/Explicit_pHREM && ./Run.pHremd

test.serial.sander.AMD:
	-cd amd && $(MAKE) -k test 
	-cd scaledMD && $(MAKE) -k test 

test.parallel.sander.AMD:
	-export TESTsander='../../../bin/sander.MPI'; cd amd && $(MAKE) -k test
	-export TESTsander='../../bin/sander.MPI'; cd scaledMD && $(MAKE) -k test

test.amoeba:
	-cd amoeba_wat1 && ./Run.amoeba_wat1
	-cd amoeba_wat2 && ./Run.amoeba_wat2
	-cd amoeba_wat2 && ./Run.ntpverlet
	-cd amoeba_gb1 && ./Run.amoeba_gb1
	-cd amoeba_jac && ./Run.amoeba_jac
	-cd amoeba_formbox && ./Run.amoeba_formbox

test.sander.RISM:
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.mdiis0; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.mdiis1; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.resize; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.dist; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.dist.xyzv; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.hnc; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.list; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.pse1; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.min; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/ala && ./Run.ala.imin5; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph1; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph1.pse1; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph1.pse3; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph1-nacl; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph1-nacl.pse1; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph1-nacl.pse3; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph2; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph10; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph10-nacl; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/phosphate && ./Run.ph20; fi

test.sander.RISM.dev: #test.sander.RISM
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/dev/ala && ./Run.ala.dist.hnc; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/dev/ala && ./Run.ala.dist.pse3; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/dev/dt/spc_spc && ./Run.all; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/dev/dt/ala_spc && ./Run.all; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/dev/dt/Na_spc && ./Run.all; fi
	-if [ -n "$(TESTRISMSANDER)" ]; then cd rism3d/dev/dt/Na_spc-NaCl && ./Run.all; fi

test.sander.RISM.MPI:
	 # -if [ -n "$(TESTRISMSANDER)" ]; then \
	 # 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 # 	&& cd rism3d/ala && ./Run.ala; \
	 # fi
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd rism3d/ala && ./Run.ala.mdiis0; \
	 fi
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd rism3d/ala && ./Run.ala.mdiis1; \
	 fi
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd rism3d/ala && ./Run.ala.resize; \
	 fi
	 # -if [ -n "$(TESTRISMSANDER)" ]; then \
	 # 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 # 	&& cd rism3d/ala && ./Run.ala.dist; \
	 # fi
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd rism3d/ala && ./Run.ala.dist.xyzv; \
	 fi
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd rism3d/ala && ./Run.ala.hnc; \
	 fi
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd rism3d/ala && ./Run.ala.pse1; \
	 fi

test.sander.RISM.OIN:
	-if [ -n "$(TESTRISMSANDER)" ]; then cd sander_rism3d_asfe_oin && ./Run.ala; fi

test.sander.RISM.OIN.MPI:
	 -if [ -n "$(TESTRISMSANDER)" ]; then \
	 	export TESTsander='$(AMBERHOME)/bin/sander.MPI' \
	 	&& cd sander_rism3d_asfe_oin && ./Run.ala; \
	 fi

test.sander.crgreloc:
	-(cd crg_reloc/water && ./Run.water.min)
	-(cd crg_reloc/water && ./Run.water.md)

test.sander.crgreloc.MPI:
	-export TESTsander='$(AMBERHOME)/bin/sander.MPI' && cd crg_reloc/water && ./Run.water.min
	-export TESTsander='$(AMBERHOME)/bin/sander.MPI' && cd crg_reloc/water && ./Run.water.md

test.sander.APBS:
	export TESTsander='$(AMBERHOME)/bin/sander.APBS'; cd iapbs_radi/ && ./Run.ion.min


test.serial.pmemd.basic:  
	-export TESTsander='../../bin/pmemd'; cd 4096wat && ./Run.pure_wat
	-export TESTsander='../../bin/pmemd'; cd 4096wat && ./Run.pure_wat_nmr_temp_reg
	-export TESTsander='../../bin/pmemd'; cd 4096wat && ./Run.vrand
	-export TESTsander='../../bin/pmemd'; cd 4096wat && ./Run.frcdmp
	-export TESTsander='../../bin/pmemd'; cd 4096wat_oct && ./Run.pure_wat_oct
	-export TESTsander='../../bin/pmemd'; cd alp && ./Run.alp
	-export TESTsander='../../bin/pmemd'; cd cytosine && ./Run.cytosine
	-export TESTsander='../../bin/pmemd'; cd dhfr && ./Run.dhfr
	-export TESTsander='../../bin/pmemd'; cd dhfr && ./Run.dhfr.min
	-export TESTsander='../../bin/pmemd'; cd dhfr && ./Run.dhfr.noshake
	-export TESTsander='../../bin/pmemd'; cd ff14ipq && ./Run.ff14ipq
	-export TESTsander='../../bin/pmemd'; cd variable_14 && ./Run.variable_14_ntb1
	-export TESTsander='../../bin/pmemd'; cd variable_14 && ./Run.variable_14_igb1
	-export TESTsander='../../bin/pmemd'; cd trx && ./Run.trx
	-export TESTsander='../../bin/pmemd'; cd trx && ./Run.trx.cpln.pmemd
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.min
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.ln
	-export TESTsander='../../bin/pmemd'; cd circ_dna && ./Run.circdna
	-export TESTsander='../../bin/pmemd'; cd gb2_trx && ./Run.trxox.nogbsa && ./Run.trxox.nogbsa prmtop_an
	-export TESTsander='../../bin/pmemd'; cd gb7_trx && ./Run.trxox_md
	-export TESTsander='../../bin/pmemd'; cd gb8_trx && ./Run.trxox
	-export TESTsander='../../bin/pmemd'; cd gb8_trx && ./Run.trxox_md && ./Run.trxox_md prmtop_an
	-export TESTsander='../../bin/pmemd'; cd gbneck2nu_test && ./Run.1hji_md
	-export TESTsander='../../bin/pmemd'; cd alpb_trx && ./Run.trxox.nogbsa
	-export TESTsander='../../bin/pmemd'; cd gb1_cox2 && ./Run.cox2
	-export TESTsander='../../bin/pmemd'; cd gbsa_xfin && ./Run.gbsa
	-export TESTsander='../../bin/pmemd'; cd tip4p && ./Run.tip4p
	-export TESTsander='../../bin/pmemd'; cd tip4p && ./Run.tip4p_mcbar
	-export TESTsander='../../bin/pmemd'; cd tip4p && ./Run.tip4p_nve
	-export TESTsander='../../bin/pmemd'; cd tip5p && ./Run.tip5p
	-export TESTsander='../../bin/pmemd'; cd tip5p && ./Run.tip5p_nve
# constant pH test
	-export TESTsander='../../../bin/pmemd'; cd cnstph/implicit && ./Run.cnstph
	-export TESTsander='../../../bin/pmemd'; cd cnstph/explicit && ./Run.cnstph
#CHARMM TESTS
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.min
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.md
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.md
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.min
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.md
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.min
	-export TESTsander='../../../../bin/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.md
#AMD TESTS
	-export TESTsander='../../../bin/pmemd'; cd amd && $(MAKE) -k test OPT=pmemd
	-export TESTsander='../../bin/pmemd'; cd scaledMD && $(MAKE) -k test OPT=pmemd
#IPS TESTS
	-export TESTsander='../../bin/pmemd'; cd gact_ips && ./Run.ips
# Constant surface tension (csurften > 0)
	-export TESTsander='../../bin/pmemd'; cd csurften && ./Run.csurften_z-dir
# Lennard-Jones 12-6-4 potential
	-export TESTsander='../../bin/pmemd'; cd lj_12-6-4 && ./Run.12-6-4
# Semi-isotropic pressure coupling (csurften > 0 AND ntp = 3)
	-export TESTsander='../../bin/pmemd'; cd csurften && ./Run.csurften_z-dir_npt_3
#    NMR Opt Test Cases
	-export TESTsander='../../../../bin/pmemd'; cd nmropt && $(MAKE) pmemd_compat
# NetCDF Traj/Restart tests
	-export TESTsander='../../../bin/pmemd'; cd netcdf && $(MAKE) -k test OPT=$(NETCDF)
    
test.serial.pmemd.TI:
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/campTI && ./Run.campTI
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/pheMTI && ./Run.0
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/pheMTI && ./Run.1
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/pheMTI && ./Run.lambda0
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/pheMTI && ./Run.lambda1
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/sodium && ./Run.sodium
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/ti_ggcc && ./Run.test1
	-export TESTsander='../../../bin/pmemd'; cd pmemdTI/ti_ggcc && ./Run.test2
	-export TESTsander='../../../../bin/pmemd'; cd pmemdTI/softcore && ./Run_sc

#2 million atom tests are not run by default because they need >3GB Memory
test.serial.big:
	-cd 2_million_atom && ./Run.2mil.min
	-cd 2_million_atom && ./Run.2mil.md
	-cd 2_million_atom && ./Run.2mil.rst.md

test.parallel.big:
	-export TESTsander='../../bin/sander.MPI'; cd 2_million_atom && ./Run.2mil.min
	-export TESTsander='../../bin/sander.MPI'; cd 2_million_atom && ./Run.2mil.md
	-export TESTsander='../../bin/sander.MPI'; cd 2_million_atom && ./Run.2mil.rst.md

test.serial.pmemd.big:
	-export TESTsander='../../bin/pmemd'; cd 2_million_atom && ./Run.2mil.min
	-export TESTsander='../../bin/pmemd'; cd 2_million_atom && ./Run.2mil.md
	-export TESTsander='../../bin/pmemd'; cd 2_million_atom && ./Run.2mil.rst.md

test.serial.sander.sgld:
	-export TESTsander='../../bin/sander'; cd gact_ips && ./Run.ips_sgld
	-export TESTsander='../../bin/sander'; cd gact_ips && ./Run.ips_sgldfp
	-export TESTsander='../../bin/sander'; cd gact_ips && ./Run.ips_sgldg
	-export TESTsander='../../bin/sander'; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander='../../bin/sander'; cd gb_rna && ./Run.gbrna.sgld
	-export TESTsander='../../bin/sander'; cd gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander='../../bin/sander'; cd gb_rna && ./Run.gbrna.sgldg
	-export TESTsander='../../bin/sander'; cd gb_rna && ./Run.gbrna.sgmdg

test.parallel.sander.sgld:
	-export TESTsander='../../bin/sander.MPI'; cd gact_ips && ./Run.ips_sgld
	-export TESTsander='../../bin/sander.MPI'; cd gact_ips && ./Run.ips_sgldfp
	-export TESTsander='../../bin/sander.MPI'; cd gact_ips && ./Run.ips_sgldg
	-export TESTsander='../../bin/sander.MPI'; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander='../../bin/sander.MPI'; cd gb_rna && ./Run.gbrna.sgld
	-export TESTsander='../../bin/sander.MPI'; cd gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander='../../bin/sander.MPI'; cd gb_rna && ./Run.gbrna.sgldg
	-export TESTsander='../../bin/sander.MPI'; cd gb_rna && ./Run.gbrna.sgmdg

test.parallel.tmp.sgld:
	-export TESTsander='../../bin/sander.MPI'; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander='../../bin/sander.MPI'; cd gb_rna && ./Run.gbrna.sgmdg

test.serial.tmp.sgld:
	-export TESTsander='../../bin/pmemd'; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.sgmdg
#	-export TESTsander='../../bin/sander'; cd gact_ips && ./Run.ips_sgmdg
#	-export TESTsander='../../bin/sander'; cd gb_rna && ./Run.gbrna.sgmdg

test.serial.pmemd.sgld:
	-export TESTsander='../../bin/pmemd'; cd gact_ips && ./Run.ips_sgld
	-export TESTsander='../../bin/pmemd'; cd gact_ips && ./Run.ips_sgldfp
	-export TESTsander='../../bin/pmemd'; cd gact_ips && ./Run.ips_sgldg
	-export TESTsander='../../bin/pmemd'; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.sgld
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.sgldg
	-export TESTsander='../../bin/pmemd'; cd gb_rna && ./Run.gbrna.sgmdg

test.parallel.pmemd.sgld:
	-export TESTsander='../../bin/pmemd.MPI'; cd gact_ips && ./Run.ips_sgld
	-export TESTsander='../../bin/pmemd.MPI'; cd gact_ips && ./Run.ips_sgldfp
	-export TESTsander='../../bin/pmemd.MPI'; cd gact_ips && ./Run.ips_sgldg
	-export TESTsander='../../bin/pmemd.MPI'; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna.sgld
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna.sgldg
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna.sgmdg

test.mic_offload.pmemd.sgld:
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gact_ips && ./Run.ips_sgld
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gact_ips && ./Run.ips_sgldfp
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gact_ips && ./Run.ips_sgldg
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gact_ips && ./Run.ips_sgmdg
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna.sgld
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna.sgldg
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna.sgmdg

test.parallel.pmemd.big:
	-export TESTsander='../../bin/pmemd.MPI'; cd 2_million_atom && ./Run.2mil.min
	-export TESTsander='../../bin/pmemd.MPI'; cd 2_million_atom && ./Run.2mil.md
	-export TESTsander='../../bin/pmemd.MPI'; cd 2_million_atom && ./Run.2mil.rst.md

test.parallel.pmemd.basic:
	-export TESTsander='../../bin/pmemd.MPI'; cd 4096wat && ./Run.pure_wat
	-export TESTsander='../../bin/pmemd.MPI'; cd 4096wat && ./Run.pure_wat_nmr_temp_reg
	-export TESTsander='../../bin/pmemd.MPI'; cd 4096wat && ./Run.vrand
	-export TESTsander='../../bin/pmemd.MPI'; cd 4096wat && ./Run.frcdmp
	-export TESTsander='../../bin/pmemd.MPI'; cd 4096wat_oct && ./Run.pure_wat_oct
	-export TESTsander='../../bin/pmemd.MPI'; cd alp && ./Run.alp
	-export TESTsander='../../bin/pmemd.MPI'; cd cytosine && ./Run.cytosine
	-export TESTsander='../../bin/pmemd.MPI'; cd dhfr && ./Run.dhfr
	-export TESTsander='../../bin/pmemd.MPI'; cd dhfr && ./Run.dhfr.min
	-export TESTsander='../../bin/pmemd.MPI'; cd dhfr && ./Run.dhfr.noshake
	-export TESTsander='../../bin/pmemd.MPI'; cd ff14ipq && ./Run.ff14ipq
	-export TESTsander='../../bin/pmemd.MPI'; cd variable_14 && ./Run.variable_14_ntb1
	-export TESTsander='../../bin/pmemd.MPI'; cd trx && ./Run.trx
	-export TESTsander='../../bin/pmemd.MPI'; cd trx && ./Run.trx.cpln.pmemd
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna.min
	-export TESTsander='../../bin/pmemd.MPI'; cd gb_rna && ./Run.gbrna.ln
	-export TESTsander='../../bin/pmemd.MPI'; cd circ_dna && ./Run.circdna
	-export TESTsander='../../bin/pmemd.MPI'; cd gb2_trx && ./Run.trxox.nogbsa
	-export TESTsander='../../bin/pmemd.MPI'; cd gb7_trx && ./Run.trxox_md
	-export TESTsander='../../bin/pmemd.MPI'; cd gb8_trx && ./Run.trxox
	-export TESTsander='../../bin/pmemd.MPI'; cd gb8_trx && ./Run.trxox_md && ./Run.trxox_md prmtop_an
	-export TESTsander='../../bin/pmemd.MPI'; cd gbneck2nu_test && ./Run.1hji_md
	-export TESTsander='../../bin/pmemd.MPI'; cd alpb_trx && ./Run.trxox.nogbsa
	-export TESTsander='../../bin/pmemd.MPI'; cd gb1_cox2 && ./Run.cox2
	-export TESTsander='../../bin/pmemd.MPI'; cd gbsa_xfin && ./Run.gbsa
	-export TESTsander='../../bin/pmemd.MPI'; cd tip4p && ./Run.tip4p
	-export TESTsander='../../bin/pmemd.MPI'; cd tip4p && ./Run.tip4p_mcbar
	-export TESTsander='../../bin/pmemd.MPI'; cd tip4p && ./Run.tip4p_nve
	-export TESTsander='../../bin/pmemd.MPI'; cd tip5p && ./Run.tip5p
	-export TESTsander='../../bin/pmemd.MPI'; cd tip5p && ./Run.tip5p_nve
# constant pH test
	-export TESTsander='../../../bin/pmemd.MPI'; cd cnstph/implicit && ./Run.cnstph
	-export TESTsander='../../../bin/pmemd.MPI'; cd cnstph/explicit && ./Run.cnstph
#CHARMM TESTS
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.min
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.md
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.md
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.min
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.md
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.min
	-export TESTsander='../../../../bin/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.md
#AMD TESTS
	-export TESTsander='../../../bin/pmemd.MPI'; cd amd && $(MAKE) -k test OPT=pmemd
	-export TESTsander='../../bin/pmemd.MPI'; cd scaledMD && $(MAKE) -k test OPT=pmemd
#IPS TESTS
	-export TESTsander='../../bin/pmemd.MPI'; cd gact_ips && ./Run.ips
# Constant surface tension (csurften > 0)
	-export TESTsander='../../bin/pmemd.MPI'; cd csurften && ./Run.csurften_z-dir
# Lennard-Jones 12-6-4 potential
	-export TESTsander='../../bin/pmemd.MPI'; cd lj_12-6-4 && ./Run.12-6-4
# Semi-isotropic pressure coupling (csurften > 0 AND ntp = 3)
	-export TESTsander='../../bin/pmemd.MPI'; cd csurften && ./Run.csurften_z-dir_npt_3
#    NMR Opt Test Cases
	-export TESTsander='../../../../bin/pmemd.MPI'; cd nmropt && $(MAKE) pmemd_compat
# NetCDF Traj/Restart tests
	-export TESTsander='../../../bin/pmemd.MPI'; cd netcdf && $(MAKE) -k test OPT=$(NETCDF)

test.mic_offload.pmemd.basic:
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd 4096wat && ./Run.pure_wat
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd 4096wat && ./Run.pure_wat_nmr_temp_reg
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd 4096wat && ./Run.vrand
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd 4096wat && ./Run.frcdmp
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd 4096wat_oct && ./Run.pure_wat_oct
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd alp && ./Run.alp
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd cytosine && ./Run.cytosine
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd dhfr && ./Run.dhfr
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd dhfr && ./Run.dhfr.min
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd dhfr && ./Run.dhfr.noshake
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd ff14ipq && ./Run.ff14ipq
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd variable_14 && ./Run.variable_14_ntb1
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd trx && ./Run.trx
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd trx && ./Run.trx.cpln.pmemd
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna.min
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb_rna && ./Run.gbrna.ln
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd circ_dna && ./Run.circdna
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb2_trx && ./Run.trxox.nogbsa
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb7_trx && ./Run.trxox_md
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb8_trx && ./Run.trxox
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb8_trx && ./Run.trxox_md && ./Run.trxox_md prmtop_an
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd alpb_trx && ./Run.trxox.nogbsa
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gb1_cox2 && ./Run.cox2
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gbsa_xfin && ./Run.gbsa
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd tip4p && ./Run.tip4p
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd tip4p && ./Run.tip4p_mcbar
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd tip4p && ./Run.tip4p_nve
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd tip5p && ./Run.tip5p
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd tip5p && ./Run.tip5p_nve
# constant pH test
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd cnstph/implicit && ./Run.cnstph
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd cnstph/explicit && ./Run.cnstph
#CHARMM TESTS
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.min
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.md
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.md
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.min
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.md
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.min
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.md
#AMD TESTS
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd amd && $(MAKE) -k test OPT=pmemd
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd scaledMD && $(MAKE) -k test OPT=pmemd
#IPS TESTS
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd gact_ips && ./Run.ips
# Constant surface tension (csurften > 0)
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd csurften && ./Run.csurften_z-dir
# Lennard-Jones 12-6-4 potential
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd lj_12-6-4 && ./Run.12-6-4
# Semi-isotropic pressure coupling (csurften > 0 AND ntp = 3)
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd csurften && ./Run.csurften_z-dir_npt_3
#    NMR Opt Test Cases
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd nmropt && $(MAKE) pmemd_compat
# NetCDF Traj/Restart tests
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd netcdf && $(MAKE) -k test OPT=$(NETCDF)

#EMIL TI TESTS
test.parallel.emil:
	-export TESTsander='../../../bin/sander.MPI'; cd emil/emil_sander_gbsa && ./Run.emil $(EMIL)
	-export TESTsander='../../../bin/sander.MPI'; cd emil/emil_sander_tip3p && ./Run.emil $(EMIL)

test.parallel.emil.pmemd:
	-export TESTsander='../../../bin/pmemd.MPI'; cd emil/emil_pmemd_gbsa && ./Run.emil $(EMIL)
	-export TESTsander='../../../bin/pmemd.MPI'; cd emil/emil_pmemd_tip3p && ./Run.emil $(EMIL)
        

test.mic_offload.emil.pmemd:
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd emil/emil_pmemd_gbsa && ./Run.emil $(EMIL)
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd emil/emil_pmemd_tip3p && ./Run.emil $(EMIL)

test.serial.emil:
	-export TESTsander='../../../bin/sander'; unset DO_PARALLEL; cd emil/emil_sander_gbsa && ./Run.emil $(EMIL)
	-export TESTsander='../../../bin/sander'; unset DO_PARALLEL; cd emil/emil_sander_tip3p && ./Run.emil $(EMIL)

test.serial.emil.pmemd:
	-export TESTsander='../../../bin/pmemd';  unset DO_PARALLEL; cd emil/emil_pmemd_gbsa && ./Run.emil $(EMIL)
	-export TESTsander='../../../bin/pmemd';  unset DO_PARALLEL; cd emil/emil_pmemd_tip3p && ./Run.emil $(EMIL)

#TI TESTS
test.parallel.pmemd.TI:
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/campTI && ./Run.campTI
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.0
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.1
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.lambda0
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.lambda1
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/sodium && ./Run.sodium
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/ti_ggcc && ./Run.test1
	-export TESTsander='../../../bin/pmemd.MPI'; cd pmemdTI/ti_ggcc && ./Run.test2
	-export TESTsander='../../../../bin/pmemd.MPI'; cd pmemdTI/softcore && ./Run_sc

test.mic_offload.pmemd.TI:
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/campTI && ./Run.campTI
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/pheMTI && ./Run.0
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/pheMTI && ./Run.1
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/pheMTI && ./Run.lambda0
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/pheMTI && ./Run.lambda1
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/sodium && ./Run.sodium
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/ti_ggcc && ./Run.test1
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/ti_ggcc && ./Run.test2
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd pmemdTI/softcore && ./Run_sc

# API tests
test.serial.sanderapi:
	-cd sanderapi && $(MAKE) test
	-cd sanderlesapi && $(MAKE) test

#emap TESTS
test.serial.sander.emap:
	-export TESTsander='../../bin/sander'; cd emap/ && ./Run.emap

test.parallel.sander.emap:
	-export TESTsander='../../bin/sander.MPI'; cd emap/ && ./Run.emap

test.serial.pmemd.emap:
	-export TESTsander='../../bin/pmemd'; cd emap/ && ./Run.emap

test.parallel.pmemd.emap:
	-export TESTsander='../../bin/pmemd.MPI'; cd emap/ && ./Run.emap

test.mic_offload.pmemd.emap:
	-export TESTsander="$(AMBERHOME)/bin/pmemd.mic_offload.MPI"; cd emap/ && ./Run.emap

test.pmemd.cuda:
	cd cuda && $(MAKE) -k test.pmemd.cuda PREC_MODEL=$(PREC_MODEL) 

test.pmemd.cuda.MPI:
	cd cuda && $(MAKE) -k test.pmemd.cuda.MPI PREC_MODEL=$(PREC_MODEL) 

test.pmemd.amoeba:
	-cd amoeba_wat1 && ./Run.amoeba_wat1.pmemd
	-cd amoeba_wat2 && ./Run.amoeba_wat2.pmemd
	-cd amoeba_wat2 && ./Run.ntpverlet.pmemd
	-cd amoeba_gb1 && ./Run.amoeba_gb1.pmemd
	-cd amoeba_jac && ./Run.amoeba_jac.pmemd
	-cd amoeba_formbox && ./Run.amoeba_formbox.pmemd
	-cd amoeba_softcore && ./Run.amoeba_softcore.pmemd

test.pmemd.amoeba.MPI:
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_wat1 && ./Run.amoeba_wat1.pmemd
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_wat2 && ./Run.amoeba_wat2.pmemd
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_wat2 && ./Run.ntpverlet.pmemd
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_gb1 && ./Run.amoeba_gb1.pmemd
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_jac && ./Run.amoeba_jac.pmemd
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_formbox && ./Run.amoeba_formbox.pmemd
	-export TESTsander='../../bin/pmemd.amoeba.MPI'; cd amoeba_softcore && ./Run.amoeba_softcore.pmemd

test.ncsu.serial:
	cd ncsu && ./run-serial.sh

test.ncsu.parallel:
	cd ncsu && ./run-parallel.sh

test.lmod:
	-cd vancomycin_lmod && ./Run.vancomycin_lmod
	-cd trpcage_lmod && ./Run.lmod
	-cd vancomycin_lmod_dock && ./Run.vancomycin_lmod_dock

test.nmode:
	-cd nmode && ./Run.ala-dipeptide
	-cd nmode && ./Run.ala-dipeptide_rgroup
	-cd newton_raph && ./Run.newton_raph

test.pbsa:
	-cd sander_pbsa_dmp    && ./Run.dmp.min
	-cd sander_pbsa_ipb2   && ./Run.110D.min
	-cd sander_pbsa_delphi && ./Run.1az6.min
	-cd sander_pbsa_lpb    && ./Run.lsolver.min
	-cd sander_pbsa_radi   && ./Run.ion.min
	-cd sander_pbsa_tsr    && ./Run.tsrb.min
	-cd sander_pbsa_frc    && ./test
	-cd sander_pbsa_decres && ./Run.pbsa_decres
	-cd sander_pbsa_decpw  && ./Run.pbsa_decpw
	-cd sander_pbsa_ligand && ./Run
	-cd sander_pbsa_imin5  && ./Run.min

test.netcdf:
	-(cd ../src/netcdf/src && $(MAKE) check)
	-cd trajene_netcdf && ./Run.trajene $(NETCDF)
	-cd netcdf && $(MAKE) -k test OPT=$(NETCDF)

test.neb:
	-cd neb-testcases/neb_gb_partial && ./Run.neb_gb_partial
	-cd neb-testcases/neb_gb_full && ./Run.neb_gb_full
	-cd neb-testcases/neb_explicit && ./Run.neb_explicit

#===========================================================================

finished.serial:
	@echo ""
	@echo "Finished serial test suite for Amber 14 at `date`."
	@echo ""

finished.parallel:
	@echo ""
	@echo "Finished parallel test suite for Amber 14 at `date`."
	@echo "Some tests require 4 threads to run, while some will not"
	@echo "run with more than 2.  Please run further parallel tests with the"
	@echo "appropriate number of processors. See $(AMBERHOME)/test/README."
	@echo ""

finished.mic_offload:
	@echo ""
	@echo "Finished mic offload test suite for Amber 14 at `date`."
	@echo "Some tests require 4 threads to run, while some will not"
	@echo "run with more than 2.  Please run further parallel tests with the"
	@echo "appropriate number of processors. See $(AMBERHOME)/test/README."
	@echo ""

is_amberhome_defined:
	@(if [ -z "$(AMBERHOME)" ] ; then \
	    echo "Error: AMBERHOME should be defined or else some tests will fail !" ;\
	    exit 2 ;\
	fi ;\
	)

clean:
	-./test_amber_clean.sh

check:
	(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

#===========================================================================
