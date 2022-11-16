#
#************************************************************************
# 
#************************************************************************
#
include ../config.h

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI
.SUFFIXES:
.SUFFIXES: .F90 .o

# PBSA Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

SRC=\
	dsvdc.F90  pbsa.F90 pb_read.F90 pb_write.F90 getcoor.F90 \
	runmd.F90 runmin.F90 force.F90 \
	pb_init.F90 sa_driver.F90 np_force.F90 pb_force.F90 pb_fddrv.F90 pb_iimdrv.F90 pb_nlsolver.F90 \
	pb_exmol.F90 pb_mpfrc.F90 pb_direct.F90 pb_list.F90 timer.F90 pb_fftsolv.F90 pb_fft.F90 \
	egb.F90 ene.F90 locmem.F90 myopen.F90 rdparm.F90 pb_lsolver.F90 \
	decomp.F90 rgroup.F90 rfree.F90 debug.F90 \
	coed20.F90  indexg.F90  irre32.F90  phidev.F90  project.F90  qld.F90 \
        coed6.F90   irre31.F90  miniop.F90 iimod.F90  prodis.F90  qint.F90     wint.F90 \
        transf.F90 curv.F90 GrToPr.F90 jumps.F90 problem.F90 \
        IIM.F90 regular.F90 gmres.F90 pbsa_lib.F90 \
	daxpy.F90  ddot.F90    dheqr.F90  dpigmr.F90  dscal.F90   dslui2.F90  dxlcal.F90 \
	dchkw.F90  dgmres.F90  dnrm2.F90  drlcal.F90  dsilus.F90  dslui.F90   isdgmr.F90 \
	dcopy.F90  dhels.F90   dorth.F90  ds2y.F90    dslugm.F90  dsmv.F90    qs2i1d.F90 \
	d1mach.F90  fdump.F90   i1mach.F90  j4save.F90  xercnt.F90  xerhlt.F90  xermsg.F90  xerprn.F90  xersve.F90  xgetua.F90 \
        bicg.F90 dbcg.F90 dslui4.F90 dsluti.F90 dsmtv.F90 dslubc.F90 isdbcg.F90 \
        pb_chkfrc.F90 svdcmp.F90 svbksb.F90 pythag.F90 fftw3.F90 memory_module.F90 \
	pb_augdrv.F90  gmresX.F90 interpX.F90 matvec3.F90 gen_dx_file.F90 aug_iccg.F90 membrane.F90 parms.F90 rdpqr.F90

OBJ=\
	dsvdc.o  pbsa.o pb_read.o pb_write.o getcoor.o \
	runmd.o runmin.o force.o pbsa_lib.o \
	pb_init.o sa_driver.o np_force.o pb_force.o pb_fddrv.o pb_iimdrv.o pb_nlsolver.o \
	pb_exmol.o pb_mpfrc.o pb_direct.o pb_list.o timer.o pb_fftsolv.o pb_fft.o \
	egb.o ene.o locmem.o myopen.o rdparm.o pb_lsolver.o \
	decomp.o rgroup.o rfree.o debug.o \
	coed20.o  indexg.o  irre32.o  phidev.o  project.o  qld.o \
	coed6.o   irre31.o  miniop.o iimod.o  prodis.o  qint.o     wint.o \
	transf.o curv.o GrToPr.o jumps.o problem.o \
	IIM.o regular.o gmres.o \
	daxpy.o  ddot.o    dheqr.o  dpigmr.o  dscal.o   dslui2.o  dxlcal.o \
	dchkw.o  dgmres.o  dnrm2.o  drlcal.o  dsilus.o  dslui.o   isdgmr.o \
	dcopy.o  dhels.o   dorth.o  ds2y.o    dslugm.o  dsmv.o    qs2i1d.o \
	d1mach.o  fdump.o   i1mach.o  j4save.o  xercnt.o  xerhlt.o  xermsg.o  \
        xerprn.o  xersve.o  xgetua.o \
	bicg.o dbcg.o dslui4.o dsluti.o dsmtv.o dslubc.o isdbcg.o \
        pb_chkfrc.o svdcmp.o svbksb.o pythag.o\
        pb_augdrv.o gmresX.o interpX.o matvec3.o gen_dx_file.o \
        aug_iccg.o membrane.o fftw3.o parms.o memory_module.o rdpqr.o

LIBPBSAOBJS= \
	interface.o pbsaapi.LIBPBSA.o   myopen.o    pb_read.LIBPBSA.o pb_force.LIBPBSA.o \
	pb_fddrv.o  pb_iimdrv.o sa_driver.LIBPBSA.o pb_nlsolver.o     pb_lsolver.o \
	dsvdc.o     timer.o     pb_init.LIBPBSA.o   pb_exmol.LIBPBSA.o        np_force.LIBPBSA.o \
	pb_direct.LIBPBSA.o pb_list.LIBPBSA.o   pb_write.LIBPBSA.o  pb_mpfrc.o        decomp.o \
	rgroup.o    rfree.o     coed20.o    indexg.o          irre32.o pbsa_lib.o \
	phidev.o    project.o   qld.o       coed6.o           irre31.o \
	miniop.o  iimod.o   prodis.o    qint.o      wint.o            transf.o \
	curv.o      GrToPr.o    jumps.o     problem.LIBPBSA.o         IIM.o \
        regular.o   pb_augdrv.o	    gmresX.o          matvec3.o\
        gmres.o     bicg.o      project.o   transf.o          dslubc.o \
        dslugm.o    curv.o      irre31.o    irre32.o          interpX.o\
	            wint.o      dchkw.o     ds2y.o            GrToPr.o \
	dgmres.o    dslui.o     dsluti.o    pb_fftsolv.LIBPBSA.o      pb_fft.LIBPBSA.o \
	dsmv.o      dsmtv.o     dsilus.o    jumps.o           dbcg.o \
	qld.o       phidev.o    dpigmr.o    dslui2.o          dslui4.o \
	isdbcg.o    qs2i1d.o    xermsg.o    d1mach.o          xercnt.o \
	xerhlt.o    drlcal.o    j4save.o    isdgmr.o          fdump.o \
	xerprn.o    dorth.o     dhels.o     xersve.o          dheqr.o \
	i1mach.o    xgetua.o    dxlcal.o    fftw3.LIBPBSA.o memory_module.o \
        svdcmp.o    svbksb.o    pythag.o aug_iccg.o gen_dx_file.LIBPBSA.o membrane.LIBPBSA.o \
	parms.LIBPBSA.o

SANDERPBSAOBJS= dsvdc.o sa_driver.SANDER.o pb_nlsolver.o \
        pb_lsolver.o aug_iccg.o pb_force.SANDER.o np_force.SANDER.o \
        pb_exmol.SANDER.o pb_mpfrc.SANDER.o pb_direct.SANDER.o \
        pb_list.SANDER.o pb_read.SANDER.o pb_write.SANDER.o pb_init.SANDER.o \
	pb_fddrv.o pb_iimdrv.o timer.SANDER.o pbsa_lib.o \
	coed20.o  indexg.o  irre32.o  phidev.o  project.o  qld.o \
	coed6.o   irre31.o  miniop.o iimod.o  prodis.o  qint.o     wint.o \
	transf.o curv.o GrToPr.o jumps.o problem.SANDER.o \
	IIM.o \
        regular.o \
	gmres.o  daxpy.o  ddot.o   dheqr.o  dpigmr.o dscal.o  dslui2.o \
        dxlcal.o dchkw.o  dgmres.o dnrm2.o  drlcal.o dsilus.o dslui.o  \
        isdgmr.o dcopy.o  dhels.o  dorth.o  ds2y.o   dslugm.o dsmv.o   \
        qs2i1d.o d1mach.o fdump.o  i1mach.o j4save.o xercnt.o xerhlt.o \
        xermsg.o xerprn.o xersve.o xgetua.o \
        bicg.o   dbcg.o   dslui4.o dsluti.o dsmtv.o  dslubc.o isdbcg.o \
	svdcmp.o svbksb.o pythag.o gen_dx_file.SANDER.o membrane.SANDER.o \
	memory_module.SANDER.o

install: install.$(INSTALLTYPE)

install.serial: serial
	/bin/mv pbsa$(SFX) $(BINDIR)

install.parallel: parallel
	/bin/mv pbsa.MPI$(SFX) $(BINDIR)

serial: configured_serial pbsa$(SFX)

parallel: configured_parallel pbsa.MPI$(SFX)

uninstall: 
	-/bin/rm -f $(BINDIR)/pbsa$(SFX) $(BINDIR)/pbsa.MPI$(SFX) $(LIBDIR)/libpbsa.a $(LIBDIR)/libFpbsa.a

libinstall: libinstall.$(INSTALLTYPE)

libinstall.openmp: libinstall.serial

libinstall.serial: libpbsa.a
	/bin/mv libpbsa.a $(LIBDIR)

libinstall.parallel: libpbsa.a
	/bin/mv libpbsa.a $(LIBDIR)
#---------------------------------------------------------------------------
configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\
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

#---------------------------------------------------------------------------
pbsa$(SFX): FFTW3 $(OBJ) syslib netlib c9x-complex configured_serial 
	$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o pbsa$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o $(FLIBS_FFTW3) \
	        -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

pbsa.MPI$(SFX): FFTW3.MPI $(OBJ) syslib netlib c9x-complex configured_parallel 
	$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o pbsa.MPI$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o $(FLIBS_FFTW3) \
	        -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

simplepbsa$(SFX): simplepbsa.o gopt.o libpbsa.a sfflib syslib netlib c9x-complex
	$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o simplepbsa$(SFX) simplepbsa.o gopt.o \
		libpbsa.a $(LIBDIR)/libsff.a ../lib/nxtsec.o -L$(LIBDIR) $(FLIBSF) \
		$(LDFLAGS) $(AMBERLDFLAGS)
	/bin/mv simplepbsa$(SFX) $(BINDIR)

libpbsa.a: FFTW3 $(LIBPBSAOBJS) syslib
	-rm -f libpbsa.a
	$(AR) libpbsa.a $(LIBPBSAOBJS) ../lib/random.o
	$(RANLIB) libpbsa.a

prelibFpbsa:
	@(if [ ! -f ../sander/decomp.o ]; then \
	     echo "Error: libFpbsa requires decomp.o from SANDER's folder"; \
	     exit 2;\
	  fi;\
	)
	/bin/rm -f decomp.o decomp.mod DECOMP.mod
	cp -p ../sander/decomp.o .
	cp -p ../sander/decomp.mod .

libFpbsa.a: prelibFpbsa libFpbsa.$(INSTALLTYPE)
	/bin/rm -f decomp.o decomp.mod DECOMP.mod

libFpbsa.serial: $(SANDERPBSAOBJS)
	$(AR) libFpbsa.a $(SANDERPBSAOBJS)
	$(RANLIB) libFpbsa.a
	/bin/mv libFpbsa.a $(LIBDIR)

libFpbsa.parallel: $(SANDERPBSAOBJS)
	$(AR) libFpbsa.a $(SANDERPBSAOBJS)
	$(RANLIB) libFpbsa.a
	/bin/mv libFpbsa.a $(LIBDIR)

FFTW3:
	@if [ ! -z "$(FFTWLIB)" ]; then \
	if [ "$(FLIBS_FFTW3)" = '-lfftw3xf_$(COMPILER)' ]; then\
		cd $(MKL)/interfaces/fftw3xf && $(MAKE) lib$(MKL_PROCESSOR) \
		compiler=$(COMPILER) PRECISION=MKL_DOUBLE \
		INSTALL_DIR=$(LIBDIR); \
	else \
	    if [ ! -f $(INCDIR)/fftw3.f03 ]; then\
		cd ../fftw-3.3 && $(MAKE) && $(MAKE) -j 1 install; \
	    fi \
	fi \
	fi

FFTW3.MPI:
	@if [ "$(PBSAFLAG)" = "-DFFTW" ]; then \
	if [ "$(FLIBS_FFTW3)" = '-lfftw3xf_$(COMPILER)' ]; then\
		cd $(MKL)/interfaces/fftw3xf && $(MAKE) lib$(MKL_PROCESSOR) \
		compiler=$(COMPILER) PRECISION=MKL_DOUBLE \
		INSTALL_DIR=$(LIBDIR); \
	else \
	    if [ ! -f $(INCDIR)/fftw3-mpi.f03 ]; then\
		cd ../fftw-3.3 && $(MAKE) && $(MAKE) -j 1 install; \
	    fi \
	fi \
	fi

syslib:
	cd ../lib && $(MAKE) nxtsec.o random.o

sfflib:
	cd ../sff && $(MAKE) install

netlib: 
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)
	cd ../arpack && $(MAKE) install

c9x-complex:
	@if test $(C9XCOMPLEX) != "skip"; then \
	   cd ../c9x-complex && $(MAKE) libmc.a; \
	 fi

pb_init.o: pb_init.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

pb_init.SANDER.o: pb_init.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -DSANDER -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

%.LIBPBSA.o: %.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

%.SANDER.o: %.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -DSANDER -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -I../sander \
	-o $@ $<

interface.o: interface.c
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

.F90.o:   $<
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

depend::
	./makedepend > depend

clean:
	/bin/rm -f pbsa$(SFX) libpbsa.a
	/bin/rm -f *.o
	/bin/rm -f 
	/bin/rm -f *.mod *__genmod.F90
	/bin/rm -f *.d

distclean: clean



#----------HEADER DEPENDENCIES (currently from makedepend) 
include depend
# DO NOT DELETE
