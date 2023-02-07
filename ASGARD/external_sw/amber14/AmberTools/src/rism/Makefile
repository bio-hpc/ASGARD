SHELL=/bin/sh
include ../config.h

# rism1d Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

# ------- rism1d information: ----------------------------------------------

OBJ=rism1d.o
SRC=rism1d.F90

.PHONY:  yes no \
	librism librism_mpi rism1d$(SFX) \
	rism3d.snglpnt$(SFX) rism3d.snglpnt.MPI$(SFX) \
	install install_mpi uninstall clean


rism1d$(SFX): $(OBJ) librism
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o rism1d$(SFX) $(OBJ) \
	        -L$(LIBDIR) -lrism $(FLIBS_FFTW3) \
	        $(FLIBSF)  $(LDFLAGS) $(AMBERLDFLAGS)  $(LIBDIR)/libxblas-amb.a

# ------- rism3d information: ----------------------------------------------

LIBOBJ = rism1d_c.o solvmdl_c.o rism1d_potential_c.o\
	rism1d_closure_c.o rism1d_kh_c.o rism1d_hnc_c.o rism1d_py_c.o \
	rism1d_mv0_c.o rism1d_psen_c.o \
	quaternion.o array_util.o rism_util.o rism_report_c.o \
	rism3d_grid_c.o rism3d_closure_c.o \
	rism3d_kh_c.o rism3d_hnc_c.o rism3d_psen_c.o \
	rism3d_c.o rism3d_potential_c.o rism3d_csv.o \
	rism3d_xyzv.o rism3d_opendx.o \
	rism3d_solv_c.o rism3d_solu_c.o \
	pubfft.o rism3d_fft.o rism_parm.o \
	mdiis_orig_c.o mdiis_blas_c.o mdiis_blas2_c.o mdiis_c.o \
	fce_c.o erfcfun.o safemem.o blend.o rism_timer_c.o constants.o \
	getopts_c.o fftw3.o #rism3d_debug_c.o

LIBOBJ_MPI = quaternion.o array_util.o rism_util.o rism_report_c.o \
	rism3d_grid_c.o rism3d_closure_c.o \
	rism3d_kh_c.o rism3d_hnc_c.o rism3d_psen_c.o \
	rism3d_c.o rism3d_potential_c.o rism3d_csv.o \
	rism3d_xyzv.o rism3d_opendx.o \
	rism3d_solv_c.o rism3d_solu_c.o \
	pubfft.o rism3d_fft.o rism_parm.o \
	mdiis_orig_c.o mdiis_blas_c.o mdiis_blas2_c.o mdiis_c.o \
	fce_c.o erfcfun.o safemem.o blend.o rism_timer_c.o constants.o \
	fftw3.o #rism3d_debug_c.o


librism: $(LIBOBJ)
	$(AR) $(LIBDIR)/$@.a $(LIBOBJ)
	$(RANLIB) $(LIBDIR)/$@.a

librism_mpi: $(LIBOBJ_MPI)
	$(AR) $(LIBDIR)/$@.a $(LIBOBJ_MPI)
	$(RANLIB) $(LIBDIR)/$@.a

# ------ Single-Point information ------------------------------------------
SPSRC=rism3d.snglpnt.nab
rism3d.snglpnt$(SFX): $(SPSRC)
	$(BINDIR)/nab $(SPSRC) -o $(BINDIR)/rism3d.snglpnt$(SFX)
rism3d.snglpnt.MPI$(SFX): $(SPSRC)
	$(BINDIR)/mpinab -v $(SPSRC) -o $(BINDIR)/rism3d.snglpnt.MPI$(SFX)
# --------------------------------------------------------------------------

yes: install

no:

install: librism rism1d$(SFX) 
	/bin/mv rism1d$(SFX)  $(BINDIR)

install_mpi: librism_mpi 

uninstall:
	-rm -f $(BINDIR)/rism1d$(SFX) \
	$(BINDIR)/rism3d.snglpnt$(SFX) \
	$(BINDIR)/rism3d.snglpnt.MPI$(SFX) \
	$(LIBDIR)/librism.a \
	$(LIBDIR)/librism_mpi.a

clean:
	-/bin/rm -f *.o  *.mod rism1d$(SFX) rism3d.thermo$(SFX) \
	rism3d.snglpnt$(SFX) rism3d.snglpnt.MPI$(SFX)

depend::
	./makedepend > depend

include depend
