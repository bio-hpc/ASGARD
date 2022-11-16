include ../config.h

SOURCE=	main.c rdparm.c dispatch.c help.c utility.c second.c \
	io.c trajectory.c netcdf_ptraj.c parallel_ptraj.c evec.c torsion.c mask.c \
	rms.c display.c interface.c energy.c \
	experimental.c ptraj.c actions.c analyze.c \
        thermo.c pubfft.c cluster.c clusterLib.c

OBJECTS=main.o rdparm.o dispatch.o help.o utility.o second.o \
	io.o trajectory.o netcdf_ptraj.o parallel_ptraj.o evec.o torsion.o mask.o \
	rms.o display.o interface.o energy.o \
	experimental.o ptraj.o actions.o analyze.o \
        thermo.o pubfft.o cluster.o clusterLib.o

HEADERS= ptraj.h ptraj_local.h rdparm.h utility.h \
	help.h dispatch.h io.h trajectory.h netcdf_ptraj.h parallel_ptraj.h evec.h torsion.h mask.h \
	vector.h rms.h display.h interface.h energy.h \
	experimental.h constants.h actions.h analyze.h version.h \
	cluster.h clusterLib.h

LIBS = -L$(LIBDIR) -lpdb $(FLIBS_PTRAJ)

.c.o:
	$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(NETCDFINC) -o $@ $<

all: rdparm$(SFX) ptraj$(SFX)

parallel: ptraj.MPI$(SFX)
	/bin/mv ptraj.MPI$(SFX) $(BINDIR)

main.o: $(HEADERS)
help.o: rdparm.h help.c help.h utility.h
dispatch.o: dispatch.c dispatch.h utility.h actions.h
utility.o: utility.c utility.h version.h
io.o: io.c io.h rdparm.h utility.h
trajectory.o: trajectory.c trajectory.h rdparm.h vector.h utility.h
netcdf_ptraj.o: netcdf_ptraj.c netcdf_ptraj.h utility.h trajectory.h rdparm.h
parallel_ptraj.o: parallel_ptraj.c parallel_ptraj.h trajectory.h
evec.o: evec.c evec.h
ptraj.o: ptraj.c ptraj.h rdparm.h vector.h actions.h dispatch.h utility.h version.h
torsion.o: torsion.c rdparm.h vector.h utility.h
rms.o: rms.c rms.h vector.h utility.h
display.o: display.c display.h rdparm.h utility.h
interface.o: interface.c interface.h dispatch.h dispatch.c rdparm.h utility.h
experimental.o: experimental.c experimental.h rdparm.h utility.h
actions.o: ptraj.h actions.h utility.h
analyze.o: ptraj.h analyze.h utility.h
rdparm.o: rdparm.c $(HEADERS)
energy.o: energy.c energy.h ptraj.h actions.h utility.h
mask.o: mask.c mask.h
cluster.o: cluster.c ptraj.h version.h contributors.h
clusterLib.o: clusterLib.c ptraj.h version.h contributors.h

thermo.o:  thermo.F90
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

pubfft.o:  pubfft.F90
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

rdparm$(SFX): libs $(OBJECTS)
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o rdparm$(SFX) $(OBJECTS) $(LIBS) $(NETCDFLIB) $(LM)

ptraj$(SFX): libs netlib $(OBJECTS)
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o ptraj$(SFX) $(OBJECTS) $(LIBS) $(NETCDFLIB) $(LM)

ptraj.MPI$(SFX): libs netlib $(OBJECTS) 
	$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o ptraj.MPI$(SFX) $(OBJECTS) $(LIBS) $(PNETCDFLIB) $(LM)

libs: 
	cd pdb && $(MAKE) 
	cd ../arpack && $(MAKE)

netlib:
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)

clean:
	cd pdb && $(MAKE) clean
	/bin/rm -f $(OBJECTS) rdparm$(SFX) ptraj$(SFX)

install: all
	/bin/mv rdparm$(SFX) ptraj$(SFX) $(BINDIR)

uninstall:
	-rm -f $(BINDIR)/rdparm$(SFX)
	-rm -f $(BINDIR)/ptraj$(SFX)
	-rm -f $(BINDIR)/ptraj$(SFX).MPI
	-rm -f $(LIBDIR)/libpdb.a

