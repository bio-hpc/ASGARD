if IS_DARWIN

if USE_QT
QT_LF=-framework QtCore -framework QtXml
endif

fileFormatsTest_LDFLAGS = -bind_at_load $(QT_LF) $(fileFormatsTest_ldflags_tmp)
ringTest_LDFLAGS = -bind_at_load $(QT_LF) $(ringTest_ldflags_tmp)
mmTest_LDFLAGS = -bind_at_load $(QT_LF) $(mmTest_ldflags_tmp)
hybridizeTest_LDFLAGS = -bind_at_load $(QT_LF) $(hybridizeTest_ldflags_tmp)
linearAlgebraTest_LDFLAGS = -bind_at_load $(QT_LF) $(linearAlgebraTest_ldflags_tmp)
statsTest_LDFLAGS = -bind_at_load $(QT_LF) $(linearAlgebraTest_ldflags_tmp)

else

fileFormatsTest_LDFLAGS = $(fileFormatsTest_ldflags_tmp)
ringTest_LDFLAGS = $(ringTest_ldflags_tmp)
mmTest_LDFLAGS = $(mmTest_ldflags_tmp)
hybridizeTest_LDFLAGS = $(hybridizeTest_ldflags_tmp)
linearAlgebraTest_LDFLAGS = $(linearAlgebraTest_ldflags_tmp)
statsTest_LDFLAGS = $(linearAlgebraTest_ldflags_tmp)

endif

if USE_LOG4CXX
LOG4CXX_LD=-L@log4cxx_inc@/lib
LOG4CXX_LIB=-llog4cxx
endif

if USE_XERCES
XERCES_LD=-L@xerces_path@/lib
XERCES_LIB=-lxerces-c
endif

if USE_QT
tmpQT_INCPATH = $(QT_INCPATH)
tmpQT_LIBS = $(QT_LIBS)
endif

if IS_LINUX
  LAPACK_LD=-L/usr/lib64
  LAPACK_LIB=-llapack -lblas -lg2c
endif

MOLECULE_LD=-L$(abs_top_builddir)/src/Molecule
MOLECULE_LIB=-lMolecule
DIAGNOSTICS_LD=-L$(abs_top_builddir)/src/Diagnostics
DIAGNOSTICS_LIB=-lDiagnostics
LOG_LD=-L$(abs_top_builddir)/src/Log
LOG_LIB=-lLog
GRAPH_LD=-L$(abs_top_builddir)/src/Graph
GRAPH_LIB=-lGraph
GA_LD=-L$(abs_top_builddir)/src/GA
GA_LIB=-lGA
MM_LD=-L$(abs_top_builddir)/src/MM
MM_LIB=-lMM
MINIMIZERS_LD=-L$(abs_top_builddir)/src/Minimizers
MINIMIZERS_LIB=-lMinimizers
PARSERS_LD=-L$(abs_top_builddir)/src/Parsers
PARSERS_LIB=-lParsers
STATS_LD=-L$(abs_top_builddir)/src/Statistics
STATS_LIB=-lStatistics
TINYXML_LD=-L$(abs_top_builddir)/src/tinyxml
TINYXML_LIB=-ltinyxml

# INCLUDES = -I@boost_inc@ -I@xerces_path@/include  -I$(top_srcdir)/src -I$(top_srcdir)/include -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/Parsers -I$(top_srcdir)/src/Molecule -I$(top_srcdir)/src/Statistics -I$(top_srcdir)/src/Utils $(tmpQT_INCPATH)

INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(top_srcdir)/src -I$(top_srcdir)/include $(tmpQT_INCPATH) -I$(top_srcdir)/src/eigen3b2

EXTRA_DIST = $(ref_sources)

AUTOMAKE_OPTIONS = -Wno-portability

noinst_PROGRAMS = \
       fileFormatsTest \
       linearAlgebraTest \
       statsTest \
       ringTest \
       mmTest \
       hybridizeTest

check_PROGRAMS = \
       fileFormatsTest \
       linearAlgebraTest \
       statsTest \
       ringTest \
       mmTest \
       hybridizeTest

#### FILEFORMATSTEST ####
fileFormatsTest_SOURCES = fileFormatsTest.cpp
fileFormatsTest_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(STATS_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
fileFormatsTest_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(STATS_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(LOG4CXX_LIB) ${FLIBS}

#### LINEARALGEBRATEST ####
linearAlgebraTest_SOURCES = linearAlgebraTest.cpp
linearAlgebraTest_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(STATS_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
linearAlgebraTest_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(STATS_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) ${FLIBS}

#### STATSTEST ####
statsTest_SOURCES = statsTest.cpp
statsTest_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(STATS_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
statsTest_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(STATS_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) ${FLIBS}

#### ringTest ####
ringTest_SOURCES = ringTest.cpp
ringTest_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(STATS_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
ringTest_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(STATS_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) ${FLIBS}

#### MMTEST ####
mmTest_SOURCES = mmTest.cpp
mmTest_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(STATS_LD) $(MM_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
mmTest_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(STATS_LIB) $(MM_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) ${FLIBS}

#### HYBRIDIZETEST ####
hybridizeTest_SOURCES = hybridizeTest.cpp
hybridizeTest_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(STATS_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
hybridizeTest_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(STATS_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) ${FLIBS}

noext_tests = $(patsubst %$(EXEEXT), %, $(check_PROGRAMS))
ref_sources = $(patsubst %$(EXEEXT), %.ref, $(check_PROGRAMS))

create: $(check_PROGRAMS)

# check-local: $(check_PROGRAMS) $(ref_sources)
# 	@for i in $(noext_tests); do \
#		echo ------------------------------------------------------; \
#		rm -rf $(top_builddir)/tests/$$i.tmp; \
#		mkdir $(top_builddir)/tests/$$i.tmp; \
#		./$$i$(EXEEXT) -l $$i.log -o $$i.out -b $(top_builddir)/tests/$$i.tmp/; \
#		diff $$i.out $(top_srcdir)/tests/$$i.ref > $$i.diff; \
#		if test -s $$i.diff; then \
#			echo "Test \`$$i' FAILED."; \
#		else \
#			echo "Test \`$$i' PASSED."; \
#		fi; \
#	done;
#	@echo ------------------------------------------------------

#clean-local:
#	rm -rf *.tmp *.diff *.log *.out hybridizeData/*.sdf temp*.xml
