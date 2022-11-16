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
	gbnsr6.F90 pb_read.F90  gb_read.F90 pb_write.F90 getcoor.F90 \
	runmd.F90 runmin.F90 force.F90 \
	pb_init.F90 sa_driver.F90 np_force.F90 variable_module.F90 \
	pb_exmol.F90 NSR6routines.F90  pb_list.F90 timer.F90 \
	egb.F90 ene.F90 locmem.F90 myopen.F90 rdparm.F90 \
	decomp.F90 rgroup.F90 rfree.F90 debug.F90 \
        svdcmp.F90 svbksb.F90 pythag.F90  memory_module.F90 \
        gen_dx_file.F90 parms.F90 

OBJ=\
	gbnsr6.o pb_read.o  gb_read.o  pb_write.o getcoor.o \
	runmd.o runmin.o force.o\
	pb_init.o sa_driver.o np_force.o variable_module.o\
	 pb_exmol.o NSR6routines.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o \
	decomp.o rgroup.o rfree.o debug.o \
        svdcmp.o svbksb.o \
        parms.o pythag.o memory_module.o \
        gen_dx_file.o

install: install.$(INSTALLTYPE)

install.serial: serial
	/bin/mv gbnsr6$(SFX) $(BINDIR)

serial: configured_serial gbnsr6$(SFX)


uninstall: 
	-/bin/rm -f $(BINDIR)/gbnsr6$(SFX) 

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

#---------------------------------------------------------------------------
gbnsr6$(SFX): $(OBJ) syslib netlib c9x-complex configured_serial 
	$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o gbnsr6$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o  \
	        -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

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
	/bin/rm -f *.o
	/bin/rm -f 
	/bin/rm -f *.mod *__genmod.F90
	/bin/rm -f *.d

distclean: clean



#----------HEADER DEPENDENCIES (currently from makedepend) 
include depend
# DO NOT DELETE
