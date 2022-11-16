include ../config.h

.c.o:
	$(CC) -c -Dflex $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(RISMSFF) $(NETCDFINC) -o $@ $<

OBJS = binpos.o conjgrad.o lmodC.o memutil.o nblist.o newton.o nmode.o \
	prm.o rand2.o sasad.o sff.o time.o xminC.o AmberNetcdf.o atomic_number.o $(SFF_RISM_INTERFACE)


install: libsff hcp_getpdb

libsff: $(OBJS)
	-rm -f $(LIBDIR)/libsff.a
	$(AR) $(LIBDIR)/libsff.a $(OBJS)	
	$(RANLIB) $(LIBDIR)/libsff.a

libsff_mpi: $(OBJS)
	-rm -f $(LIBDIR)/libsff_mpi.a
	$(AR) $(LIBDIR)/libsff_mpi.a $(OBJS)	
	$(RANLIB) $(LIBDIR)/libsff_mpi.a

uninstall:
	-rm -f $(BINDIR)/hcp_getpdb
	-rm -f $(LIBDIR)/libsff.a
	-rm -f $(LIBDIR)/libsff_mpi.a

clean:
	-/bin/rm -f $(OBJS) lex.mm_options.c

lex.mm_options.c: mm_options.l
	$(LEX) -t mm_options.l | sed -e 's/{stdout}/{DONTchangeSTDOUThere}/' -e 's/yy/mmo/g' -e 's/stdout/nabout/' -e 's/DONTchangeSTDOUThere/stdout/' > lex.mm_options.c

#  dependencies (hand-created for now):

binpos.o: memutil.h

conjgrad.o: memutil.h

lmodC.o: sff.h

memutil.o: sff.h

nblist.o: sff.h

newton.o: sff.h memutil.h

nmode.o: sff.h memutil.h

prm.o:   sff.h

rand2.c:  sff.h

sasad.o: sff.h

atomic_number.o: atomic_number.h

sff.o: sff.h memutil.h timer.h\
       debug.h eff_box.c dssp.c eff.c \
       gbsa.c egbr6.c sff2.c mask.c rattle.c lex.mm_options.c hcp.c hcp_gb.c timer.c

time.o:  sff.h

xminC.o:  sff.h

AmberNetcdf.o:  AmberNetcdf.h

hcp_getpdb: hcp_getpdb.c
	$(CC) -o $(BINDIR)/hcp_getpdb hcp_getpdb.c

.PHONY: ../rism/amber_rism_interface.NAB.o
../rism/amber_rism_interface.NAB.o:
	cd ../rism && \
	$(MAKE) amber_rism_interface.NAB.o


