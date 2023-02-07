#if DARWIN_ACCEL
#libGA_la_LDFLAGS = -framework Accelerate
#else
#LAPACK_LIB = -L/usr/lib64 -llapack -lblas -lg2c
#endif
#libGA_la_LIBADD = $(LAPACK_LIB)

XERCES_LIB=-L@xerces_path@/lib/ -lxerces-c
# old # INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(srcdir)/../Parsers -I$(srcdir)/../Utils

INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(srcdir)/..

lib_LTLIBRARIES = libGA.la

libGA_la_SOURCES = \
gaWorld.cpp \
gaRegion.cpp \
gaPopulation.cpp \
gaIndividual.cpp \
gaChromosome.cpp \
gaGene.cpp \
gaOperators.cpp \
gaAverage.cpp \
gaCrossOver.cpp \
gaMutate.cpp \
gaGaussian.cpp \
gaSelection.cpp \
gaOutput.cpp

libGAincludedir = $(includedir)/GA

libGAinclude_HEADERS = \
gaWorld.h \
gaRegion.h \
gaPopulation.h \
gaIndividual.h \
gaChromosome.h \
gaGene.h \
gaOperators.h \
gaAverage.h \
gaCrossOver.h \
gaMutate.h \
gaGaussian.h \
gaSelection.h \
gaOutput.h \
gaVersion.h

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in
