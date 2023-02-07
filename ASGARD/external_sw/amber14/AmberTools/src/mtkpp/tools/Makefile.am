##### FLAGS #####
if IS_DARWIN
if USE_QT
QT_LF=-framework QtCore -framework QtXml
endif
MTKppConstants_LDFLAGS = -bind_at_load $(QT_LF) $(MTKppConstants_ldflags_tmp)
# tester_LDFLAGS = -bind_at_load $(QT_LF) $(tester_ldflags_tmp)
else
MTKppConstants_LDFLAGS = $(MTKppConstants_ldflags_tmp)
# tester_LDFLAGS = $(tester_ldflags_tmp)
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

DIAGNOSTICS_LD = -L$(abs_top_builddir)/src/Diagnostics
DIAGNOSTICS_LIB = -lDiagnostics
MOLECULE_LD = -L$(top_builddir)/src/Molecule
MOLECULE_LIB = -lMolecule
LOG_LD = -L$(top_builddir)/src/Log
LOG_LIB = -lLog
GRAPH_LD = -L$(top_builddir)/src/Graph
GRAPH_LIB = -lGraph
GA_LD = -L$(top_builddir)/src/GA
GA_LIB = -lGA
MM_LD = -L$(top_builddir)/src/MM
MM_LIB = -lMM
MINIMIZERS_LD = -L$(top_builddir)/src/Minimizers
MINIMIZERS_LIB = -lMinimizers
PARSERS_LD = -L$(top_builddir)/src/Parsers
PARSERS_LIB = -lParsers
STATS_LD = -L$(top_builddir)/src/Statistics
STATS_LIB = -lStatistics
TINYXML_LD = -L$(top_builddir)/src/tinyxml
TINYXML_LIB = -ltinyxml

# INCLUDES = -I@boost_inc@ -I@xerces_path@/include  -I$(top_srcdir)/src -I$(top_srcdir)/include -I$(top_srcdir)/src/tinyxml -I$(top_srcdir)/src/Parsers -I$(top_srcdir)/src/Statistics -I$(top_srcdir)/src/Utils $(QT_INCPATH)

INCLUDES = -I@boost_inc@ -I@xerces_path@/include -I$(top_srcdir)/src -I$(top_srcdir)/include $(tmpQT_INCPATH) -I$(top_srcdir)/src/eigen3b2

##### ADD PROGRAM NAME HERE ####
bin_PROGRAMS = MTKppConstants \
               capActiveSite \
               frcmod2xml \
               func \
               hybrid \
               mmE \
               prep2xml \
               protonator \
               sequenceAligner \
               stats \
               superimposer \
               stdLib2Sdf \
               pdbSearcher \
               MCPB
#               tester

#### MTKPPCONSTANTS ####
MTKppConstants_SOURCES = MTKppConstants.cpp
MTKppConstants_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(tmpQT_LIBS) $(LOG4CXX_LD)
MTKppConstants_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB)

#### TESTER ####
# tester_SOURCES = tester.cpp
# tester_ldflags_tmp = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(tmpQT_LIBS) $(LOG4CXX_LD)
# tester_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB)  $(tmpQT_LIBS) $(LOG4CXX_LIB)

#### CAPACTIVESITE ####
capActiveSite_SOURCES = capActiveSite.cpp
capActiveSite_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
capActiveSite_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### FRCMOD2XML ####
frcmod2xml_SOURCES = frcmod2xml.cpp
frcmod2xml_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
frcmod2xml_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### FUNC ####
func_SOURCES = func.cpp
func_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
func_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### HYBRID ####
hybrid_SOURCES = hybrid.cpp
hybrid_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
hybrid_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### MME ####
mmE_SOURCES = mmE.cpp
mmE_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(MM_LD) $(MINIMIZERS_LD) $(LOG_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
mmE_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(MM_LIB) $(MINIMIZERS_LIB) $(LOG_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### PREP2XML ####
prep2xml_SOURCES = prep2xml.cpp
prep2xml_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
prep2xml_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### PROTONATOR ####
protonator_SOURCES = protonator.cpp
protonator_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
protonator_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### SEQUENCEALIGNER ####
sequenceAligner_SOURCES = sequenceAligner.cpp
sequenceAligner_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(STATS_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
sequenceAligner_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(STATS_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(LOG4CXX_LIB) $(tmpQT_LIBS) ${FLIBS}

#### STATS ####
stats_SOURCES = stats.cpp
stats_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(STATS_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
stats_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(STATS_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(LOG4CXX_LIB) $(tmpQT_LIBS) ${FLIBS}

#### SUPERIMPOSER ####
superimposer_SOURCES = superimposer.cpp
superimposer_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
superimposer_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### STDLIB2SDF ####
stdLib2Sdf_SOURCES = stdLib2Sdf.cpp
stdLib2Sdf_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(LOG_LD) $(GRAPH_LD) $(GA_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
stdLib2Sdf_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(LOG_LIB) $(GRAPH_LIB) $(GA_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### PDBSEARCHER ####
pdbSearcher_SOURCES = pdbSearcher.cpp
pdbSearcher_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(MM_LD) $(MINIMIZERS_LD) $(LOG_LD) $(GRAPH_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
pdbSearcher_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(MM_LIB) $(MINIMIZERS_LIB) $(LOG_LIB) $(GRAPH_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

#### MCPB ####
MCPB_SOURCES = MCPB.cpp
MCPB_LDFLAGS = $(DIAGNOSTICS_LD) $(MOLECULE_LD) $(MM_LD) $(MINIMIZERS_LD) $(LOG_LD) $(GRAPH_LD) $(STATS_LD) $(PARSERS_LD) $(TINYXML_LD) $(XERCES_LD) $(LOG4CXX_LD)
MCPB_LDADD = $(DIAGNOSTICS_LIB) $(MOLECULE_LIB) $(MM_LIB) $(MINIMIZERS_LIB) $(LOG_LIB) $(GRAPH_LIB) $(STATS_LIB) $(PARSERS_LIB) $(TINYXML_LIB) $(XERCES_LIB) $(tmpQT_LIBS) $(LOG4CXX_LIB) ${FLIBS}

##### CLEAN UP #####
clean-local:
	rm -f *.tmp
