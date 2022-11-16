#************************************************************************
# Makefile for AmberTools14 MM_PBSA
#
SHELL=/bin/sh
include ../config.h

LIBSRC= ../lib/amopen.F ../lib/mexit.F

LIBOBJ= ../lib/amopen.o ../lib/mexit.o

PROGS = make_crd_hg$(SFX) mm_pbsa_nabnmode$(SFX)

all: $(PROGS)

make_crd_hg$(SFX): make_crd_hg.o $(LIBOBJ)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o make_crd_hg$(SFX) make_crd_hg.o \
		$(LIBOBJ) $(LDFLAGS) $(AMBERLDFLAGS)

#Note dependency on nab from AMBERTools here.	
mm_pbsa_nabnmode$(SFX): mm_pbsa_nabnmode.nab
	$(BINDIR)/nab $(NABFLAGS) -o mm_pbsa_nabnmode mm_pbsa_nabnmode.nab

../lib/amopen.o: ../lib/amopen.F
	cd ../lib; $(MAKE) amopen.o

../lib/mexit.o: ../lib/mexit.F
	cd ../lib; $(MAKE) mexit.o

install: all
	mv $(PROGS) $(BINDIR)
	cp mm_pbsa.pl mm_pbsa_statistics.pl $(BINDIR)

clean:
	/bin/rm -f *.o  $(PROGS) mm_pbsa_nabnmode.c

uninstall: clean
	-rm -f $(BINDIR)/mm_pbsa.pl
	-rm -f $(BINDIR)/mm_pbsa_statistics.pl
	-for file in $(PROGS) ; do \
	rm -f $(BINDIR)/$${file} ; \
	done
