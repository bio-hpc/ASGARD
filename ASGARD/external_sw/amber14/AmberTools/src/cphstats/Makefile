include ../config.h

PROGNAME = $(BINDIR)/cphstats$(SFX)

.NOTPARALLEL:clean install all

OBJS = main.o cpin.o parse_cpin.o string_manip.o cloptions.o test.o cpout.o \
       prottraj.o utilities.o conprob.o

all: $(PROGNAME)

$(PROGNAME): $(OBJS)
	$(CXX) -o $(PROGNAME) $(OBJS) $(LDFLAGS) $(FOPTFLAGS) $(AMBERLDFLAGS) $(ZLIB) $(NOFORTRANMAIN) 

install: all

uninstall:
	/bin/rm -f $(PROGNAME)

.SUFFIXES: .F90 .cpp

depends::
	./makedepends

.F90.o:
	$(FC) $(FFLAGS) $(FOPTFLAGS) -c $*.F90

.cpp.o:
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(CXXOPTFLAGS) -c $*.cpp

clean:
	/bin/rm -f $(OBJS)

include depends
