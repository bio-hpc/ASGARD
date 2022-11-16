DIAGNOSTICS_LD=-L$(abs_top_builddir)/src/Diagnostics/$(SRCDIR)
DIAGNOSTICS_LIB=-lDiagnostics
INCLUDES = -I$(srcdir)/..

lib_LTLIBRARIES = libMinimizers.la

libMinimizers_la_SOURCES = \
lbfgs.cpp

libMinimizersincludedir = $(includedir)/Minimizers

libMinimizersinclude_HEADERS = \
lbfgs.h

if DARWIN_ACCEL
libMinimizers_la_LDFLAGS = -framework Accelerate
else
libMinimizers_la_LDFLAGS =
endif

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in
