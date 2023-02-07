SHELL=/bin/sh
SRCDIR=../..
include ../../config.h

makeDIST_RST$(SFX):	makeDIST_RST.o
		$(CC) $(CFLAGS) $(AMBERCFLAGS) -o makeDIST_RST$(SFX) makeDIST_RST.o \
			$(LDFLAGS) $(AMBERLDFLAGS) $(LM)

nef_to_RST$(SFX):	nef_to_RST.o $(LIBDIR)/libcifparse.a
		$(CC) $(CFLAGS) $(AMBERCFLAGS) -o nef_to_RST$(SFX) nef_to_RST.o \
			$(LDFLAGS) $(AMBERLDFLAGS) $(LM) -L$(LIBDIR) -lcifparse

makeANG_RST$(SFX):	makeANG_RST.o
		$(CC) $(CFLAGS) $(AMBERCFLAGS) -o makeANG_RST$(SFX) makeANG_RST.o \
			$(LDFLAGS) $(AMBERLDFLAGS) $(LM)

install: 	makeDIST_RST$(SFX) makeANG_RST$(SFX) nef_to_RST$(SFX)
		mv nef_to_RST$(SFX) makeDIST_RST$(SFX) makeANG_RST$(SFX) $(BINDIR)
		cp makeCHIR_RST makeDIP_RST.cyana $(BINDIR)
		cp map.NEF-AMBER map.DG-AMBER $(DATDIR)

uninstall:
		-rm -f $(BINDIR)/makeDIST_RST$(SFX)
		-rm -f $(BINDIR)/makeANG_RST$(SFX)
		-rm -f $(BINDIR)/makeCHIR_RST
		-rm -f $(BINDIR)/makeDIP_RST.cyana
		-rm -f $(DATDIR)/map.DG-AMBER

clean:
	/bin/rm -f *.o 

makeDIST_RST.o:  global.h
nef_to_RST.o:  global.h

#  default rule for C compilation:

.c.o:
	$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

$(LIBDIR)/libcifparse.a:
	cd ../../cifparse && $(MAKE) install
