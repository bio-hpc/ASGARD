DIAGNOSTICS_LD=-L$(abs_top_builddir)/src/Diagnostics/$(SRCDIR)
DIAGNOSTICS_LIB=-lDiagnostics

INCLUDES = -I@boost_inc@  -I$(srcdir)/.. -I$(top_srcdir)/src/eigen3b2

lib_LTLIBRARIES = libStatistics.la

libStatistics_la_SOURCES = \
sheet.cpp \
BaseStats.cpp \
ols.cpp \
pca.cpp \
pls.cpp


libStatisticsincludedir = $(includedir)/Statistics

libStatisticsinclude_HEADERS = \
table.h \
sheet.h \
BaseStats.h \
ols.h \
pca.h \
pls.h

if DARWIN_ACCEL
libStatistics_la_LDFLAGS = -framework Accelerate
else
libStatistics_la_LDFLAGS =
endif

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in

libStatistics_la_LIBADD = $(FLIBS)
