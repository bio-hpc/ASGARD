DIAGNOSTICS_LD=-L$(abs_top_builddir)/src/Diagnostics/$(SRCDIR)
DIAGNOSTICS_LIB=-lDiagnostics
LOG_LD=-L$(abs_top_builddir)/src/Log/$(SRCDIR)
LOG_LIB=-lLog

INCLUDES = -I@boost_inc@ -I$(srcdir)/..

lib_LTLIBRARIES = libMM.la

libMM_la_SOURCES = \
mmPotential.cpp \
amber.cpp \
amberBond.cpp \
amberAngle.cpp \
amberTorsion.cpp \
amberImproper.cpp \
amberNonBonded.cpp

libMMincludedir = $(includedir)/MM

libMMinclude_HEADERS = \
mmPotential.h \
amber.h \
amberBond.h \
amberAngle.h \
amberTorsion.h \
amberImproper.h \
amberNonBonded.h 

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in

if DARWIN_ACCEL
libMM_la_LDFLAGS = -framework Accelerate $(LOG_LD)
else
libMM_la_LDFLAGS = $(LOG_LD)
endif

libMM_la_LIBADD = $(LOG_LIB)
