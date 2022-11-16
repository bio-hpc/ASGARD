# INCLUDES = -I$(srcdir)/../Utils

lib_LTLIBRARIES = libGraph.la

INCLUDES = -I@boost_inc@  -I$(srcdir)/..

libGraph_la_SOURCES = \
vertex.cpp \
graph.cpp

libGraphincludedir = $(includedir)/Graph
#$(includedir)/Diagnostics

libGraphinclude_HEADERS = \
edge.h \
vertex.h \
graph.h

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in
