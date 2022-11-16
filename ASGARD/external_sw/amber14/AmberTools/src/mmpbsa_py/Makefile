# Makefile for Python reimplementation of MM/PB(GB)SA written by Jason Swails
include ../config.h

SHELL = /bin/sh

install: $(BINDIR)/mmpbsa_py_nabnmode$(SFX) $(BINDIR)/mmpbsa_py_energy$(SFX)
	$(PYTHON) setup.py install --prefix=$(AMBERHOME) --install-scripts=$(BINDIR)

$(BINDIR)/mmpbsa_py_nabnmode$(SFX): mmpbsa_entropy.nab
	$(BINDIR)/nab -o $(BINDIR)/mmpbsa_py_nabnmode$(SFX) mmpbsa_entropy.nab

$(BINDIR)/mmpbsa_py_energy$(SFX): mmpbsa_energy.nab
	$(BINDIR)/nab -o $(BINDIR)/mmpbsa_py_energy$(SFX) mmpbsa_energy.nab

serial: install

parallel:
	-$(PYTHON) setup.py install --par --prefix=$(AMBERHOME) --install-scripts=$(BINDIR)

docs:
	doxygen mmpbsa.Doxyfile

clean:
	-(/bin/rm -fr build/ mmpbsa_entropy.c mmpbsa_energy.c)

uninstall:
	/bin/rm -fr $(BINDIR)/mmpbsa_py_nabnmode $(BINDIR)/mmpbsa_py_energy
	/bin/rm -fr $(BINDIR)/MMPBSA.py $(BINDIR)/MMPBSA.py.MPI $(BINDIR)/ante-MMPBSA.py
