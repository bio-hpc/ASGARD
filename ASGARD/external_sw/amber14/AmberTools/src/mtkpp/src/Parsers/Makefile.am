DIAGNOSTICS_LD=-L$(abs_top_builddir)/src/Diagnostics/$(SRCDIR)
DIAGNOSTICS_LIB=-lDiagnostics
MOLECULE_LD=-L$(abs_top_builddir)/src/Molecule/$(SRCDIR)
MOLECULE_LIB=-lMolecule
GA_LD=-L$(abs_top_builddir)/src/GA/$(SRCDIR)
GA_LIB=-lGA
LOG_LD=-L$(abs_top_builddir)/src/Log/$(SRCDIR)
LOG_LIB=-lLog
STATS_LD=-L$(abs_top_builddir)/src/Statistics/$(SRCDIR)
STATS_LIB=-lStatistics

TINYXML_LD=-L$(abs_top_builddir)/src/tinyxml/$(SRCDIR)
TINYXML_LIB=-ltinyxml

if USE_XERCES
XERCES_LD=-L@xerces_path@/lib/
XERCES_LIB=-lxerces-c
endif

lib_LTLIBRARIES = libParsers.la

if USE_QT
moc_baseParser.cpp: baseParser.h
	$(QTDIR)/bin/moc $(QT_DEFINES) $(QT_INCPATH) $(srcdir)/baseParser.h -o moc_baseParser.cpp

baseParser_moc = \
moc_baseParser.cpp

moc_elementParser.cpp: elementParser.h
	$(QTDIR)/bin/moc $(QT_DEFINES) $(QT_INCPATH) $(srcdir)/elementParser.h -o moc_elementParser.cpp

elementParser_moc = \
moc_elementParser.cpp

moc_stdLibParser.cpp: stdLibParser.h
	$(QTDIR)/bin/moc $(QT_DEFINES) $(QT_INCPATH) $(srcdir)/stdLibParser.h -o moc_stdLibParser.cpp

stdLibParser_moc = \
moc_stdLibParser.cpp

moc_paramParser.cpp: paramParser.h
	$(QTDIR)/bin/moc $(QT_DEFINES) $(QT_INCPATH) $(srcdir)/paramParser.h -o moc_paramParser.cpp

paramParser_moc = \
moc_paramParser.cpp

moc_dMParser.cpp: dMParser.h
	$(QTDIR)/bin/moc $(QT_DEFINES) $(QT_INCPATH) $(srcdir)/dMParser.h -o moc_dMParser.cpp

dMParser_moc = \
moc_dMParser.cpp

moc_mtkppParser.cpp: mtkppParser.h
	$(QTDIR)/bin/moc $(QT_DEFINES) $(QT_INCPATH) $(srcdir)/mtkppParser.h -o moc_mtkppParser.cpp

mtkppParser_moc = \
moc_mtkppParser.cpp

INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(srcdir)/.. $(QT_INCPATH) $(QT_DEFINES)
tmpQT_LIBS = QT_LIBS

else

INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(srcdir)/.. -I$(top_srcdir)/src/eigen3b2
endif

# INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(srcdir)/../tinyxml -I$(srcdir)/../Molecule 
# -I$(srcdir)/../GA -I$(srcdir)/../Statistics -I$(srcdir)/../Utils -I$(srcdir)/../Log $(QT_INCPATH) $(QT_DEFINES)

libParsers_la_SOURCES = \
parsingException.cpp \
baseParser.cpp \
$(baseParser_moc) \
stdLibParser.cpp \
$(stdLibParser_moc) \
commLineOptions.cpp \
molParser.cpp \
mol2Parser.cpp \
pdbParser.cpp \
pamParser.cpp \
dcParser.cpp \
sdfParser.cpp \
xyzParser.cpp \
gaussianParser.cpp \
zmatParser.cpp \
prepParser.cpp \
acParser.cpp \
frcmodParser.cpp \
amberParser.cpp \
atomTypeParser.cpp \
gaParser.cpp \
StringManip.cpp \
elementParser.cpp \
$(elementParser_moc) \
paramParser.cpp \
$(paramParser_moc) \
dMParser.cpp \
$(dMParser_moc) \
mtkppParser.cpp \
$(mtkppParser_moc)

libParsersincludedir = $(includedir)/Parsers

libParsersinclude_HEADERS = \
parsingException.h \
stdLibParser.h \
baseParser.h \
mtkppParser.h \
commLineOptions.h \
inputParser.h \
molParser.h \
mol2Parser.h \
pdbParser.h \
pamParser.h \
dcParser.h \
sdfParser.h \
xyzParser.h \
elementParser.h \
paramParser.h \
atomTypeParser.h \
gaParser.h \
StringManip.h \
domErrorHandler.h \
gaussianParser.h \
zmatParser.h \
prepParser.h \
acParser.h \
frcmodParser.h \
dMParser.h \
amberParser.h \
xmlConvertors.h

DISTCLEANFILES = Makefile .deps
MAINTAINERCLEANFILES = Makefile.in

libParsers_la_LDFLAGS = $(MOLECULE_LD) $(STATS_LD) $(LOG_LD) $(GA_LD) $(TINYXML_LD) $(XERCES_LD)
libParsers_la_LIBADD = $(MOLECULE_LIB) $(STATS_LIB) $(LOG_LIB) $(GA_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(QT_LIBS)
