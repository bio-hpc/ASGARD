#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
# Makefile for Amber 14   SANDER 
#
include ../config.h

# PLUMED-related include file
ifdef PLUMED_INCLUDE_FILE
include $(PLUMED_INCLUDE_FILE)
endif

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI

SHELL=/bin/sh
VPATH=

# Sander Fortran source files are free format
# Add also search paths for .mod files
LOCALFLAGS = $(FREEFORMAT_FLAG) -I../pbsa -I../sqm -I../rism -I../../../include


EMILOBJ= \
    emil.o mdin_emil_dat.o

NCSU_OBJECTS = \
    ncsu-abmd-hooks.o ncsu-bbmd-ctxt.o ncsu-bbmd-hooks.o ncsu-cftree.o \
    ncsu-colvar.o ncsu-colvar-math.o ncsu-colvar-type.o ncsu-colvar-utils.o \
    ncsu-constants.o ncsu-cv-ANGLE.o ncsu-cv-COM_ANGLE.o \
    ncsu-cv-COM_DISTANCE.o ncsu-cv-COM_TORSION.o ncsu-cv-COS_OF_DIHEDRAL.o \
    ncsu-cv-DISTANCE.o ncsu-cv-HANDEDNESS.o ncsu-cv-LCOD.o \
    ncsu-cv-MULTI_RMSD.o ncsu-cv-N_OF_BONDS.o ncsu-cv-N_OF_STRUCTURES.o \
    ncsu-cv-R_OF_GYRATION.o ncsu-cv-TORSION.o ncsu-lexer.o ncsu-parser.o \
    ncsu-pmd-hooks.o ncsu-rmsd.o ncsu-sander-hooks.o ncsu-sander-proxy.o \
	ncsu-cv-PCA.o ncsu-read-pca.o \
    ncsu-smd-hooks.o ncsu-umbrella.o ncsu-utils.o ncsu-value.o mt19937.o

EVB_BASE = evb_vars.o evb_mathf.o evb_init.o evb_bcast.o \
    evb_dealloc.o check_input.o evb_io.o evb_keywrd.o \
    exchange_gauss.o exchange_warshel.o react_flux_init.o \
    react_path.o react_flux.o kappa_keck.o \
    schlegel_dg.o schlegel_full.o schlegel_diis.o schlegel_gmres.o \
    schlegel_gbasis.o schlegel_pgbasis.o \
    bnd_strtch.o ang_bend.o dihed_tors.o torinternal.o \
    evb_gfchk.o bmat_inv.o bmat_grad.o \
    unitv.o crossv.o cart2internal.o wdc_bmat.o select_coord.o \
    schlegel_uff.o schlegel_dihed.o schlegel_angle.o \
    schlegel_poh_uff.o schlegel_poh_irc.o schlegel_hcn_cart.o \
    schlegel_hcn_int_irc.o schlegel_hcn_cart_irc.o evb_outoptions.o

EVB_DEBUG = evb_2stdebug.o egap_umb_2stdebug.o dg_debug.o \
    evb_io_debug.o bmat_debug.o bmat_inv_debug.o evb_io_debug2.o \
    morse_anal2num.o xwarshel_anal2num.o xgauss_anal2num.o \
    vdw_anal2num.o rcdbonds_anal2num.o rcbond_anal2num.o \
    dg_grad_anal2num.o

EVBOBJ =  $(EVB_BASE) $(EVB_DEBUG) evb_input.o evb_alloc.o evb_matrix.o \
    morsify.o mod_vdw.o evb_ntrfc.o evb_force.o out_evb.o \
    dbonds_anal2num.o bond_anal2num.o evb_umb.o evb_umb_primitive.o

EVBPIMD =  $(EVB_BASE) $(EVB_DEBUG) evb_input.LES.o evb_alloc.LES.o evb_matrix.LES.o \
    morsify.LES.o mod_vdw.LES.o evb_ntrfc.LES.o evb_force.LES.o \
    evb_pimd_init.LES.o out_evb.LES.o dbonds_anal2num.o bond_anal2num.o \
    evb_umb.LES.o evb_umb_primitive.LES.o

AMOEBAOBJ= \
    amoeba_valence.o amoeba_multipoles.o amoeba_recip.o amoeba_interface.o\
    amoeba_direct.o amoeba_mdin.o amoeba_adjust.o amoeba_self.o\
    amoeba_vdw.o amoeba_induced.o amoeba_runmd.o 

APBSOBJ= \
    file_io_dat.APBS.o apbs_vars.APBS.o apbs.APBS.o \
    constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    multisander.o sander.APBS.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.APBS.o runmin.o rdparm.o \
    mdread.APBS.o locmem.o runmd.o getcor.o degcnt.o decnvh.o \
    fastwt.o parallel.o shake.o ene.o mdwrit.o minrit.o \
    set.o dynlib.APBS.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_bspline.o\
    ew_fft.o nonbond_list.o short_ene.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o amopen.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.o thermo_int.o \
    matinv.o assert.o erfcfun.o \
    constantph.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(AMOEBAOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o sander_lib.o amd.o scaledMD.o barostats.o random.o $(EMILOBJ)

QMOBJ= \
    abfqmmm_module.o qm_mm.o qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o qm2_extern_gau_module.o \
    qm2_extern_orc_module.o qm2_extern_nw_module.o qm2_extern_qc_module.o qm2_extern_genmpi_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o 

ifdef LIOLIBS
QMOBJ += qm2_extern_lio_module.o
endif

QMAPIOBJ = \
    abfqmmm_module.o qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o qm2_extern_gau_module.o \
    qm2_extern_orc_module.o qm2_extern_nw_module.o qm2_extern_qc_module.o qm2_extern_genmpi_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o 

ifdef LIOLIBS
QMAPIOBJ += qm2_extern_lio_module.o
endif

SEBOMDOBJ = \
    sebomd_module.o sebomd_arrays.o

SEBOMDLIB = ../sebomd/sebomd.a

MMOBJ=  file_io_dat.o constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    multisander.o sander.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o rdparm.o sander_lib.o \
    mdread.o locmem.o runmd.o getcor.o degcnt.o decnvh.o \
    fastwt.o parallel.o shake.o ene.o mdwrit.o minrit.o \
    set.o dynlib.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_bspline.o\
    ew_fft.o nonbond_list.o short_ene.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o amopen.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.o thermo_int.o \
    matinv.o assert.o erfcfun.o \
    constantph.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o\
    $(EMILOBJ) $(AMOEBAOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o amd.o scaledMD.o barostats.o random.o les.o

APIOBJ=  file_io_dat.o constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o rdparm.API.o cpparm.o sander_lib.o \
    mdread.o locmem.o runmd.o getcor.o degcnt.o decnvh.o \
    fastwt.API.o parallel.o shake.o ene.o mdwrit.o minrit.o \
    set.API.o dynlib.API.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.API.o ew_box.API.o ew_bspline.o\
    ew_fft.o nonbond_list.API.o short_ene.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o amopen.o \
    debug.API.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.API.o extra_pts.API.o thermo_int.o les.o \
    matinv.o assert.o erfcfun.o \
    constantph.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o\
    $(EMILOBJ) $(AMOEBAOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.API.o crg_reloc.o amd.o scaledMD.o barostats.o interface.o random.o qm_mm.API.o prmtop_type.o

LESAPIOBJ= memory_module.o stack.o file_io_dat.o les.o \
    constants.o state.o nose_hoover.o nose_hoover_init_LES.o \
    trace.o rmsgrd.o rand2.o lmodC.o xminC.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o relax_mat.o nmr.LES.o multitmd.o \
    cshf.o nmrcal.o mtmdcall.o pearsn.o cpparm.LES.o \
    printe.o runmin.o nonbond_list.LESAPI.o force.LES.o rdparm.LESAPI.o \
    mdread.LES.o locmem.LES.o runmd.LES.o degcnt.LES.o decnvh.o \
    fastwt.API.o parallel.LES.o shake.o ene.LES.o \
    dynlib.LESAPI.o nmlsrc.o \
    ew_force.LES.o ew_setup.LESAPI.o ew_bspline.o\
    ew_fft.o short_ene.LES.o ew_recip.LES.o pcshift.o \
    align.o csa.o rfree.o rgroup.o \
    amopen.o sander_lib.o random.o \
    debug.API.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.API.o extra_pts.LESAPI.o thermo_int.o matinv.o assert.o \
    erfcfun.o \
    constantph.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(AMOEBAOBJ) $(EMILOBJ) $(SANDER_RISM_INTERFACE) \
    $(SEBOMDOBJ) \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.API.o crg_reloc.o barostats.o prmtop_type.o \
    AmberNetcdf.o bintraj.LES.o binrestart.LES.o getcor.LES.o mdfil.LES.o mdwrit.LES.o \
    lmod.LES.o trajene.o ew_box.LESAPI.o minrit.LES.o set.LESAPI.o qm_mm.LESAPI.o \
    amd.o scaledMD.o interface.LESAPI.o

FULLPIMDOBJ = pimd_vars.o pimd_force.o pimd_init.o cmd_vars.o cmd_matrix.o 

PARTPIMDOBJ = pimd_vars.o pimd_force.LES.o pimd_init.LES.o cmd_vars.o cmd_matrix.LES.o

PUPILOBJ= file_io_dat.o constants.o state.o memory_module.o stack.o \
	nose_hoover.o nose_hoover_init.o \
	trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
	pupildata.o icosasurf.o egb.o findmask.o \
	relax_mat.o nmr.o multitmd.o \
	multisander.o sander.PUPIL.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
	printe.PUPIL.o runmin.o force.PUPIL.o rdparm.o  \
	mdread.PUPIL.o locmem.o runmd.o getcor.o degcnt.o decnvh.o \
	fastwt.o parallel.o shake.o ene.o mdwrit.o minrit.o \
	set.o dynlib.PUPIL.o mdfil.PUPIL.o nmlsrc.o \
	ew_force.o ew_setup.o ew_box.o ew_bspline.o\
	ew_fft.o nonbond_list.o short_ene.o ew_recip.o pcshift.o align.o csa.o \
	rfree.o rgroup.o amopen.o sander_lib.o \
	debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
	mexit.PUPIL.o new_time.o extra_pts.o thermo_int.o \
	matinv.o assert.o erfcfun.o \
	constantph.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
	$(EMILOBJ) $(AMOEBAOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
	spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o \
	putvalues.o fixport.o  dssp.o charmm.o crg_reloc.o amd.o scaledMD.o \
   barostats.o random.o

LSCIVROBJ = lscivr_vars.o lsc_init.o lsc_xp.o

LESOBJ= memory_module.o stack.o file_io_dat.o \
    constants.o state.o nose_hoover.o nose_hoover_init_LES.o \
    trace.o rmsgrd.o rand2.o lmodC.o xminC.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o \
    relax_mat.o nmr.LES.o multitmd.o \
    multisander.LES.o sander.LES.o  cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o nonbond_list.LES.o force.LES.o rdparm.LES.o \
    mdread.LES.o locmem.LES.o runmd.LES.o degcnt.LES.o decnvh.o \
    fastwt.o parallel.LES.o shake.o ene.LES.o \
    dynlib.LES.o nmlsrc.o \
    ew_force.LES.o ew_setup.LES.o ew_bspline.o\
    ew_fft.o short_ene.LES.o ew_recip.LES.o pcshift.o \
    align.o csa.o rfree.o rgroup.o \
    amopen.o sander_lib.o random.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.LES.o thermo_int.o matinv.o assert.o \
    erfcfun.o \
    constantph.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(AMOEBAOBJ) $(EMILOBJ) $(SANDER_RISM_INTERFACE) \
    $(SEBOMDOBJ) les.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o barostats.o \
    qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o qm2_extern_gau_module.o \
    qm2_extern_orc_module.o qm2_extern_nw_module.o qm2_extern_qc_module.o qm2_extern_genmpi_module.o \
    qmmm_adaptive_module.o abfqmmm_module.o prn_qmmm_dipole.o \
    AmberNetcdf.o bintraj.LES.o binrestart.LES.o getcor.LES.o mdfil.LES.o mdwrit.LES.o \
    lmod.LES.o trajene.o ew_box.LES.o minrit.LES.o set.LES.o qm_mm.LES.o amd.o scaledMD.o
ifdef LIOLIBS
LESOBJ += qm2_extern_lio_module.o
endif

XRAY_OBJS = \
   xray_fourier.o xray_globals.o xray_interface.o \
   xray_real_space.o xray_reciprocal_space.o xray_utils.o

.PHONY: configured configured-serial configured_parallel serial parallel

install: $(INSTALLTYPE)

SERIALPROGS = $(BINDIR)/sander$(SFX) $(BINDIR)/sander.LES$(SFX) $(BINDIR)/ambmask$(SFX)

.NOTPARALLEL: $(SERIALPROGS)

serial: configured_serial $(SERIALPROGS) $(LIBDIR)/libsander$(SHARED_SUFFIX) \
        $(LIBDIR)/libsanderles$(SHARED_SUFFIX) $(INCDIR)/sander.h 

PARALLELPROGS = $(BINDIR)/sander.MPI$(SFX) $(BINDIR)/sander.LES.MPI$(SFX)

.NOTPARALLEL: $(PARALLELPROGS)

parallel: configured_parallel $(PARALLELPROGS)

all_serial_programs: $(SERIALPROGS) 
	$(MAKE) $(BINDIR)/sander.PUPIL$(SFX)
	$(MAKE) $(BINDIR)/sander.APBS$(SFX)

#---------------------------------------------------------------------------
$(BINDIR)/sander$(SFX): libsqm $(MMOBJ) $(QMOBJ) $(FULLPIMDOBJ) $(EMIL) \
              $(SEBOMDOBJ) sebomd \
              $(LSCIVROBJ) force.o syslib \
              ../lib/nxtsec.o netlib configured_serial \
              $(NCSU_OBJECTS) $(XRAY_OBJS) libpbsa librism $(PLUMED_DEPENDENCIES)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(MMOBJ) $(QMOBJ) $(SEBOMDOBJ) \
	      $(FULLPIMDOBJ) $(LSCIVROBJ) force.o -L$(LIBDIR) -lsqm \
	      $(NCSU_OBJECTS) $(XRAY_OBJS) \
	      -lFpbsa ../lib/nxtsec.o $(EMILLIB) \
	      $(SEBOMDLIB) \
	      ../lib/sys.a $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
	      $(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.MPI$(SFX): libsqm $(MMOBJ) $(QMOBJ) $(FULLPIMDOBJ) $(EMIL) \
              $(LSCIVROBJ) $(EVBOBJ) force.o \
              $(SEBOMDOBJ) sebomd $(XRAY_OBJS) \
              syslib ../lib/nxtsec.o netlib libpbsa librism_mpi configured_parallel $(NCSU_OBJECTS) \
              $(PLUMED_DEPENDENCIES)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(MMOBJ) $(QMOBJ) $(SEBOMDOBJ) \
              $(FULLPIMDOBJ) $(LSCIVROBJ) $(EVBOBJ) force.o \
	      -L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	      $(SEBOMDLIB) $(XRAY_OBJS) \
	      ../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
	      $(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)


#---------------------------------------------------------------------------
$(BINDIR)/sander.PUPIL$(SFX): libsqm $(PUPILOBJ) $(QMOBJ) $(FULLPIMDOBJ) $(EMIL) \
		$(LSCIVROBJ)  syslib $(XRAY_OBJS) \
                $(SEBOMDOBJ) sebomd \
		../lib/nxtsec.o netlib libpbsa librism configured_serial $(NCSU_OBJECTS) \
                $(PLUMED_DEPENDENCIES)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(PUPILOBJ) $(QMOBJ) $(SEBOMDOBJ) $(FULLPIMDOBJ) \
		$(LSCIVROBJ) -L$(LIBDIR) -lsqm -lFpbsa \
	        $(SEBOMDLIB) $(XRAY_OBJS) \
		../lib/nxtsec.o $(EMILLIB) ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIBF) \
	    	$(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
		$(PUPILLIBS) $(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)


#---------------------------------------------------------------------------
$(BINDIR)/sander.LES$(SFX): libsqm $(LESOBJ) $(PARTPIMDOBJ) syslib \
        ../lib/nxtsec.o netlib \
		$(LSCIVROBJ) $(NCSU_OBJECTS) $(XRAY_OBJS) configured_serial \
		libpbsa librism $(EMIL) sebomd $(PLUMED_DEPENDENCIES)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(LESOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		$(XRAY_OBJS) -L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	        $(SEBOMDLIB) \
		../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIBF) \
	      	$(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
		$(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.LES.MPI$(SFX): libsqm $(LESOBJ) $(EVBPIMD) \
        $(PARTPIMDOBJ) syslib ../lib/nxtsec.o \
	    netlib libpbsa librism_mpi $(EMIL) $(XRAY_OBJS) \
		$(LSCIVROBJ) $(NCSU_OBJECTS) sebomd configured_parallel $(PLUMED_DEPENDENCIES)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(LESOBJ) $(EVBPIMD) \
        $(PARTPIMDOBJ) $(LSCIVROBJ) $(XRAY_OBJS) \
		-L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	        $(SEBOMDLIB) \
		../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIBF) \
	      	$(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.APBS$(SFX): libsqm $(APBSOBJ) $(QMOBJ) \
              $(LSCIVROBJ) $(FULLPIMDOBJ) $(NCSU_OBJECTS) \
              force.APBS.o syslib librism $(EMIL) \
              $(SEBOMDOBJ) sebomd $(XRAY_OBJS) \
              ../lib/nxtsec.o netlib configured_serial $(PLUMED_DEPENDENCIES)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(APBSOBJ) $(QMOBJ) $(SEBOMDOBJ) \
              $(FULLPIMDOBJ) $(LSCIVROBJ) force.APBS.o \
              -L$(APBS_LIBDIR) $(APBS_LIBS) \
              -L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	        $(SEBOMDLIB) $(XRAY_OBJS) \
              ../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
              $(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/ambmask$(SFX): 	ambmask.o findmask.o amopen.o parms.o \
       memory_module.o mexit.o ../lib/nxtsec.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ ambmask.o findmask.o \
		amopen.o mexit.o ../lib/nxtsec.o constants.o memory_module.o parms.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS)


#-----------LIBS

$(LIBDIR)/libsander$(SHARED_SUFFIX): libsqm $(APIOBJ) $(QMAPIOBJ) $(FULLPIMDOBJ) $(EMIL) \
              $(SEBOMDOBJ) sebomd $(INCDIR)/sander_api.mod \
              $(LSCIVROBJ) force.o syslib $(PLUMED_DEPENDENCIES) \
              ../lib/nxtsec.o netlib configured_serial \
              $(NCSU_OBJECTS) $(XRAY_OBJS) libpbsa librism
	$(FC) $(MAKE_SHARED) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(APIOBJ) \
		  $(QMAPIOBJ) $(SEBOMDOBJ) $(FULLPIMDOBJ) $(LSCIVROBJ) force.o \
		  -L$(LIBDIR) -lsqm \
	      $(NCSU_OBJECTS) $(XRAY_OBJS) \
	      -lFpbsa ../lib/nxtsec.o $(EMILLIB) \
	      $(SEBOMDLIB) \
	      ../lib/sys.a $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) \
	      $(LIOLIBS) $(PLUMED_LOAD)

$(LIBDIR)/libsanderles$(SHARED_SUFFIX): libsqm $(LESAPIOBJ) $(PARTPIMDOBJ) $(QMAPIOBJ) syslib \
        ../lib/nxtsec.o netlib $(INCDIR)/sanderles_api.mod \
		$(LSCIVROBJ) $(NCSU_OBJECTS) $(XRAY_OBJS) configured_serial \
		libpbsa librism $(EMIL) sebomd $(PLUMED_DEPENDENCIES)
	$(FC) $(MAKE_SHARED) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) \
       -o $@ $(LESAPIOBJ) $(QMAPIOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		$(XRAY_OBJS) -L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	        $(SEBOMDLIB) $(PLUMED_LOAD) \
		../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIBF) \
	      	$(FLIBS_RISMSANDER) $(FLIBS_FFTW3) $(FLIBSF) $(LIOLIBS)

$(INCDIR)/sander_api.mod: interface.o
	/bin/cp sander_api.mod $(INCDIR)

$(INCDIR)/sanderles_api.mod: interface.o
	/bin/cp sanderles_api.mod $(INCDIR)

$(INCDIR)/sander.h: sander.h
	/bin/cp sander.h $(INCDIR)

EMIL:
	$(MAKE) -C ../emil install

syslib:
	cd ../lib && $(MAKE) sys.a

netlib:  
	cd ../blas && $(MAKE) $(BLAS)
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../arpack && $(MAKE) install

sebomd: $(SEBOMDOBJ)
	cd ../sebomd && $(MAKE) sebomd.a

../lib/nxtsec.o:  ../lib/nxtsec.F
	cd ../lib && $(MAKE) nxtsec.o

pmemd_clib:
	cd ../pmemd/src && $(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) \
		$(CFLAGS) $(AMBERCFLAGS) pmemd_clib.c

../lib/mexit.o:  ../lib/mexit.F
	cd ../lib && $(MAKE) mexit.o

ifndef PLUMED_INCLUDE_FILE
Plumed.o: Plumed.c
	$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
	        -o Plumed.o Plumed.c
endif

rand2.o:  ../sff/rand2.c
	$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o rand2.o ../sff/rand2.c

lmodC.o:  ../sff/lmodC.c
	$(CC) -c -DSQM $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o lmodC.o ../sff/lmodC.c

xminC.o:  ../sff/xminC.c
	$(CC) -c -DSQM $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o xminC.o ../sff/xminC.c

libsqm:	
	cd ../sqm && $(MAKE) $@

librism librism_mpi:
	if [ -n "$(RISMSANDER)" ]; then \
		cd ../rism && $(MAKE) $@; \
	fi

libpbsa: decomp.o
	cd ../pbsa && $(MAKE) libFpbsa.a

clean:
	/bin/rm -f sander.PUPIL$(SFX)
	/bin/rm -f sander.APBS$(SFX)
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *.d
	/bin/rm -f *nbflag
	/bin/rm -f *__genmod.f90  # intel compiler-generated interface modules

uninstall.serial:
	/bin/rm -f $(BINDIR)/mopac.sh
	-for file in $(SERIALPROGS) ; do \
	  /bin/rm -f $${file} ; \
	done

uninstall.parallel:
	-for file in $(PARALLELPROGS) ; do \
	  /bin/rm -f $${file} ; \
	done

uninstall: uninstall.serial uninstall.parallel

depend::
	./makedepend > depend

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) and do NOT specify -mpi." ;\
		exit 2 ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
			exit 2 ;\
		} ;\
	)
#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
