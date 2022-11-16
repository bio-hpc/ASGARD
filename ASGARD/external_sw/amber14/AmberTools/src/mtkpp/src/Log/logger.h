#ifndef _LOGGER_H_
#define _LOGGER_H_

#ifdef USE_LOG4CXX
/** @file Declarations for logging functions and macros. */

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include "log4cxx/propertyconfigurator.h"
#include <log4cxx/consoleappender.h>
#include <log4cxx/patternlayout.h>

/// The default logger name. Use this if none of the special cases below apply.
#define LOGGER_NAME "mtkpp"

/// Logger name for Parser code.
#define LOGGER_PARSERS "mtkpp.parsers"

/// Logger name for Molecule code.
#define LOGGER_MOLECULE "mtkpp.backbone"

/// Logger name for GA code.
#define LOGGER_GA "mtkpp.ga"

/// Logger name for GRAPH code.
#define LOGGER_GRAPH "mtkpp.graph"

/// Logger name for MINIMIZERS code.
#define LOGGER_MINIMIZERS "mtkpp.minimizers"

/// Logger name for MM code.
#define LOGGER_MM "mtkpp.mm"

/// Logger name for STATISTICS code.
#define LOGGER_STATISTICS "mtkpp.statistics"

#endif // USE_LOG4CXX

#endif // _LOGGER_H_
