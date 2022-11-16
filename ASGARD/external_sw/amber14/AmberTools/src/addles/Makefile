include ../config.h

COMMON_FILES = 	LINE.h TOP.h TTOP.h \
		MISC.h UNITS.h DEBUG.h \
		T3TOP.h SIZE.h TEMP.h

LES_FILES = 	lesmain.o addspace.o readprm.o \
		writprm.o readcrd.o writcrd.o  \
		pick.o  	 rline.o   nxt.o \
		intgr.o   find.o    of.o \
		geti.o    unit.o    getc.o \
		alert.o   echo.o    get4c.o \
		getd.o    wlesprm.o lesprm.o random.o \
		les2prm.o checksz.o addspace_atm1st.o

LIBSRC=../lib/random.F90 ../lib/mexit.F ../lib/nxtsec.F

LIBOBJ=../lib/mexit.o  ../lib/nxtsec.o
LIBOBJNOPATH=mexit.o  nxtsec.o


$(BINDIR)/addles$(SFX) : $(LES_FILES) libobj
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) \
		-o $(BINDIR)/addles$(SFX) $(LES_FILES) \
		$(LIBOBJ) $(LDFLAGS) $(AMBERLDFLAGS)


#  default rules for Fortran and C compilation:

.F.o:   $<
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

.F90.o:   $<
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

#-----------LIBS

libobj:
	cd ../lib && $(MAKE) $(LIBOBJNOPATH)
		
install: $(BINDIR)/addles$(SFX)

clean:
	-rm -f *.o

uninstall:
	-rm -f $(BINDIR)/addles$(SFX)

$(LES_FILES) : $(COMMON_FILES)

random.o: random.F90 random.h
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ random.F90
