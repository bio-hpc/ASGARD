########################################################################
#                                                                      #
#                    Makefile for sebomd                               #
#                                                                      #
########################################################################

include ../config.h

SANDERDIR=../sander

PROG    = sebomd.a

MODOBJS =  \
         se_corrections_params.o \
         se_corrections_module.o \
         se_corrections_tools.o \
         se_inter.o \


SEOBJS = \
         se_atmchg.o \
         se_backinbox.o \
         se_balance.o \
         block.o \
         se_bpair.o \
         se_bsort.o \
         se_bspline.o \
         se_clustsub.o \
         se_debug.o \
         se_deeprec.o \
         se_denful.o \
         se_densub.o \
         se_diag.o \
         se_diagp.o \
         se_diat.o \
         se_dihedrpbc.o \
         se_dipole.o \
         se_dmxinter.o \
         se_doferm.o \
         se_doscf.o \
         se_edriver.o \
         se_elecsr.o \
         se_energy.o  \
         se_escf.o \
         se_esqr.o \
         se_etimer.o \
         se_fock.o \
         se_fockpme.o \
         se_fourier.o \
         se_fshift.o \
         se_gcart.o \
         se_gdiat.o \
         se_gensub.o \
         se_getpar.o \
         se_glbpnt.o \
         se_goverlp.o \
         se_grepul.o \
         se_ijpair.o \
         se_info_from_sander.o \
         se_initp.o \
         se_mosub.o \
         se_opnfil.o \
         se_overlp.o \
         se_pbcgcres.o \
         se_pbcxyz.o \
         se_pme_calcb.o \
         se_pme_calcchg.o \
         se_pme_calccm1.o \
         se_pme_calcq.o \
         se_pme_calctheta.o \
         se_pme_derec.o \
         se_pme_direct2.o \
         se_pme_qm_direct.o \
         se_pme_qm_recip.o \
         se_pme_recip.o \
         se_pme_setup.o \
         se_pmix.o \
         se_rdall.o \
         se_rdelem.o \
         se_rdnum.o \
         se_rdtail.o \
         se_repul.o \
         se_resclustsub.o \
         se_mpi_allreduce.o \
         se_mpi_vars.o \
         sebomd.o \
         se_corrections_info_from_sander.o \
         se_corrections.o \
         se_init.o \
         se_switch.o \
         se_setbox.o \
         se_setunit.o \
         se_setup.o \
         se_write_charges.o \
         se_sproc.o \
         se_whatis.o \
         se_wrtdmx.o 

.SUFFIXES:
.SUFFIXES: .o .F .c .F90

all: ${PROG}

rdtail.o: rdtail.F
	@echo "      [SEBOMD] [F77] [NOOPT] $<"
	@$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) $(RISMSANDER) -o $@ $<

.F90.o:
	@echo "      [SEBOMD] [F90] [OPT] $<"
	@$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) $(RISMSANDER) -I$(SANDERDIR) -o $@ $<

.F.o:
	@echo "      [SEBOMD] [F77] [OPT] $<"
	@$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) $(RISMSANDER) -I$(SANDERDIR) -o $@ $<

${PROG}: ${MODOBJS} ${SEOBJS}
	@echo "      [SEBOMD] [AR] $@"
	@$(AR) $(PROG) ${MODOBJS} $(SEOBJS) 
	@$(RANLIB) $(PROG)

clean:
	rm -f ${PROG} ${MODOBJS} ${SEOBJS} *.mod

depend::
	./makedepend > depend

#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE

