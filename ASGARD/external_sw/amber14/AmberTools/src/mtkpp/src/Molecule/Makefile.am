DIAGNOSTICS_LD=-L$(abs_top_builddir)/src/Diagnostics/$(SRCDIR)
DIAGNOSTICS_LIB=-lDiagnostics
GRAPH_LD=-L$(abs_top_builddir)/src/Graph/$(SRCDIR)
GRAPH_LIB=-lGraph
LOG_LD=-L$(abs_top_builddir)/src/Log/$(SRCDIR)
LOG_LIB=-lLog

INCLUDES = -I@boost_inc@ -I$(srcdir)/.. -I$(top_srcdir)/src/eigen3b2

lib_LTLIBRARIES = libMolecule.la

libMolecule_la_SOURCES = \
collection.cpp \
molecule.cpp \
submolecule.cpp \
atom.cpp \
element.cpp \
connections.cpp \
ring.cpp \
conformer.cpp \
parameters.cpp \
atomType.cpp \
atomTyper.cpp \
stdLibrary.cpp \
stdGroup.cpp \
stdFrag.cpp \
protonate.cpp \
proProtonate.cpp \
ligProtonate.cpp \
watProtonate.cpp \
fingerPrint.cpp \
functionalize.cpp \
hydrophobize.cpp \
hybridize.cpp \
pharmacophore.cpp \
superimpose.cpp \
selection.cpp \
metalCenter.cpp \
seqAlign.cpp \
complex.cpp

libMoleculeincludedir = $(includedir)/Molecule

libMoleculeinclude_HEADERS = \
collection.h \
molecule.h \
submolecule.h \
atom.h \
element.h \
connections.h \
bond.h \
angle.h \
improper.h \
torsion.h \
ring.h \
conformer.h \
atomTyper.h \
atomType.h \
parameters.h \
stdLibrary.h \
stdGroup.h \
stdFrag.h \
protonate.h \
proProtonate.h \
ligProtonate.h \
watProtonate.h \
fingerPrint.h \
functionalize.h \
hydrophobize.h \
hybridize.h \
pharmacophore.h \
superimpose.h \
selection.h \
metalCenter.h \
seqAlign.h \
complex.h \
utility.h \
version.h

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in

libMolecule_la_LDFLAGS = $(GRAPH_LD) $(LOG_LD) $(DIAGNOSTICS_LD)
libMolecule_la_LIBADD = $(GRAPH_LIB) $(LOG_LIB) $(DIAGNOSTICS_LIB) @FLIBS@
