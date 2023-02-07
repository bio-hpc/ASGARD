include ../../config.h

.c.o:
	$(CC) -c $(COPTFLAGS) $(CFLAGS) -o $@ $<

OBJS = crdmanip.o grid.o listRead.o matrix.o pdbRead.o QuatAlign.o ran2.o \
       topRead.o vector.o

install: libXtal.a
	mv libXtal.a $(LIBDIR)

libXtal.a: $(OBJS)
	-rm -f libXtal.a
	$(AR) libXtal.a $(OBJS)	
	$(RANLIB) libXtal.a

uninstall:
	-rm -f $(LIBDIR)/libXtal.a

clean:
	-/bin/rm -f $(OBJS) libXtal.a
