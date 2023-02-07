include ../../config.h

TARGET    = phenix_amber_interface.so

LIBFLAGS  = -lnetcdf -lfftw3 -lmdgx -lboost_python -lpython$(PYTHON_VER)
#  last two libraries above are in /usr/lib on my Ubuntu 12.04 distro
#  to get boost:
#  sudo apt-get install libboost-all-dev

OBJS = getmdgxfrc.o phenix_amber_interface.o

yes:
	@echo "NOTE:  Build Amber_Phenix manually by typing either:"
	@echo "       'make withphenixsource' or 'make withoutphenixsource'."

no:
	@echo "NOTE:  Skipping the building of Amber_Phenix."
	@echo "       (To build Amber_Phenix rerun configure with '-phenix'.)"

withphenixsource: INCDIR += $(PHENIX_REPO)
withphenixsource: INCDIR += $(PHENIX_INC)
withphenixsource: CXXFLAGS += -DAMBERPHENIX
withphenixsource: $(TARGET)
withoutphenixsource: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -shared -o $(TARGET) $(OBJS) -L$(LIBDIR) $(LIBFLAGS)

phenix_amber_interface.o: phenix_amber_interface.cpp
	$(CXX) -fPIC $(CXXFLAGS) -I$(INCLUDE_PY) -I$(INCDIR) -c $<

getmdgxfrc.o: getmdgxfrc.c $(INCDIR)/mdgx.h
	$(CC) -fPIC $(CFLAGS) -I$(INCDIR) -c $<

phenix: withphenixsource
nophenix: withoutphenixsource

clean:
	/bin/rm -f $(OBJS) $(TARGET)

.PHONY: clean no nophenix phenix withoutphenixsource withphenixsource yes

