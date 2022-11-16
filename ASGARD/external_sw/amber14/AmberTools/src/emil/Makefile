include ../config.h
include emilfiles

OBJECTS=$(EMIL_CPPSOURCES:.C=.o)

all: install
yes: install
install: $(LIBDIR)/libemil.a

# Only build the static library: the executable
# which used to serve as a test harness for the lib
# hasn't been needed for a while. 
$(LIBDIR)/libemil.a: $(OBJECTS) ../config.h
	$(AR) $(LIBDIR)/libemil.a $(OBJECTS)
	$(RANLIB) $(LIBDIR)/libemil.a

# All objects depend on all headers.
$(OBJECTS): $(EMIL_CPPHEADERS) ../config.h

.C.o:  
	$(CXX) -c $(CXXOPTFLAGS) $(EMIL_MPIFLAGS) -o $@ $<


#.C.o:
#	$(CXX) -Wall -c -g $(EMIL_MPIFLAGS) -o $@ $<


.PHONY:clean
clean:
	rm -f *.o *~ \#* 
