include ../../config.h

.c.o:
	$(CC) -c $(COPTFLAGS) $(CFLAGS) -o $@ $< -I../Clibs/

.NOTPARALLEL: libXtal

install: libXtal AddToBox$(SFX) ChBox$(SFX) \
         PropPDB$(SFX) UnitCell$(SFX)
	mv AddToBox ChBox PropPDB UnitCell ${BINDIR}

AddToBox$(SFX): AddToBox.o $(LIBDIR)/libXtal.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

ChBox$(SFX): ChBox.o $(LIBDIR)/libXtal.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

CheckMD$(SFX): CheckMD.o $(LIBDIR)/libXtal.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

PropPDB$(SFX): PropPDB.o $(LIBDIR)/libXtal.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

UnitCell$(SFX): UnitCell.o $(LIBDIR)/libXtal.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

pdb2rdf$(SFX): pdb2rdf.o $(LIBDIR)/libXtal.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(LIBDIR)/libXtal.a: libXtal

libXtal:
	cd ../Clibs && $(MAKE) install

uninstall: clean
	-rm $(BINDIR)/AddToBox $(BINDIR)/ChBox $(BINDIR)/CheckMD \
	 $(BINDIR)/PropPDB $(BINDIR)/UnitCell
	-( cd ../Clibs && $(MAKE) uninstall )

clean:
	-rm *.o
	-( cd ../Clibs && $(MAKE) clean )
