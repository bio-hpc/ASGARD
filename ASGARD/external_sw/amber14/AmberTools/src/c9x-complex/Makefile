include ../config.h

# Makefile for C9X complex variables library
# Edit endian-ness and computer type properly in mconf.h before compiling.
# For additional test programs see also celefunt.mak.

# Caution: the tests contain some references to the cephes math library
# or references to standard float or long double functions that might
# not be available.
# This should be enough to run at least the double precision functions:
LIBS = stubs.o -lm
# These are intended to get the cephes libraries:
#LIBS = -lml -lmd -lmf
#LIBS = ../libml.a ../libmd.a ../libmf.a

# float precision
SFILES = cmplxf.o clogf.o cgammaf.o
# double precision
DFILES = cmplx.o clog.o cgamma.o
# long double precision
#LFILES = 
LFILES = cmplxl.o clogl.o cgammal.o
#  LIBMCFILES = $(SFILES) $(DFILES) $(LFILES) stubs.o
LIBMCFILES = $(DFILES) stubs.o

all: libmc.a

libmc.a: $(LIBMCFILES)
	rm -f libmc.a
	$(AR) libmc.a $(LIBMCFILES)
	ranlib libmc.a
	mv libmc.a $(LIBDIR)

uninstall:
	-rm -f $(LIBDIR)/libmc.a

clog.o: clog.c complex.h protos.h
cmplx.o: cmplx.c complex.h protos.h
clogf.o: clogf.c complex.h protos.h
cmplxf.o: cmplxf.c complex.h protos.h
clogl.o: clogl.c complex.h protos.h
cmplxl.o: cmplxl.c complex.h protos.h
#clog-nans.o: clog-nans.c complex.h protos.h
#cmplx-nans.o: cmplx-nans.c complex.h protos.h

clean:
	rm -f *.o

skip:
	@echo "using built-in (C99) complex.h functionality"

