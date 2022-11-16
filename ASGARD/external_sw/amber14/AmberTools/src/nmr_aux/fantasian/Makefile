include ../../config.h

OBJ = fantasian.o ../../lib/random.o ../../lib/mexit.o
EXE = fantasian$(SFX)

../../lib/random.o:
	cd ../../lib && $(MAKE) random.o

../../lib/mexit.o:
	cd ../../lib && $(MAKE) mexit.o

$(EXE):		$(OBJ) 
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $(EXE) $(OBJ) $(LDFLAGS) $(AMBERLDFLAGS)

install:	$(EXE)
	mv $(EXE) $(BINDIR)

clean:
	/bin/rm -f *.o  $(EXE)

uninstall:
	-rm -f $(BINDIR)/$(EXE)

#  default rule for Fortran compilation:

.F.o:   $<
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<
