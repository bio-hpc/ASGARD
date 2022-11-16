include ../config.h

install: all
	/bin/mv saxs${SFX} saxs_md${SFX} $(BINDIR)

all: saxs$(SFX) saxs_md$(SFX)

saxs$(SFX): saxs.o sphere_lebedev_rule.o
	$(CXX) $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) \
	   -o $@ saxs.o sphere_lebedev_rule.o

saxs_md$(SFX): saxs_md.o sphere_lebedev_rule.o
	$(CXX) $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) \
	   -o $@ saxs_md.o sphere_lebedev_rule.o

.cpp.o:
	$(CXX) -c $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) -I$(INCDIR) -o $@ $<
	
sphere_lebedev_rule.cpp:  sphere_lebedev_rule.hpp

clean:
	/bin/rm -f *.o

uninstall:
	/bin/rm -f $(BINDIR)/saxs$(SFX) $(BINDIR)/saxs_md$(SFX)

