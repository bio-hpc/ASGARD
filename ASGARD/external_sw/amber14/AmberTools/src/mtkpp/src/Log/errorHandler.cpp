/*! 
   \file errorHandler.cpp
   \brief Error Handling within MTK++
   \author Martin Peters

   $Date: 2010/03/29 20:26:54 $
   $Revision: 1.3 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   This file is part of MTK++.

   MTK++ is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   MTK++ is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lessser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ----------------------------------------------------------------------------
*/

#include <sstream>

#include "errorHandler.h"
#include "error.h"

#ifdef USE_LOG4CXX

#include "logger.h"
/*
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
*/
#endif

namespace MTKpp
{

//! Error Logger
errorHandler errorLogger;

// ============================================================
// Function : errorHandler()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
errorHandler::errorHandler()
{
    this->level = 4;
    this->stream = &std::clog;
    this->log = true;
    this->maxMessages = 100;
    this->nMessages = 0;

#ifdef USE_LOG4CXX
    // Configure logging
    std::string MTKppDIR = ".";
    std::string logConfigureDir = ".";
    if (getenv("MTKppDIR")) {
      MTKppDIR = getenv("MTKppDIR");
      logConfigureDir = MTKppDIR + "/data/etc";
    }

/*
    boost::filesystem::path logConfigurePath("./mtkpp-log.config");
    if (boost::filesystem::exists(logConfigurePath)) {
      logConfigureDir = ".";
    }
*/

    //log4cxx::LayoutPtr layout(new log4cxx::PatternLayout("%d  %F:%L%n  %-5p %m%n"));
    //log4cxx::AppenderPtr appender(new log4cxx::ConsoleAppender(layout));
    //log4cxx::BasicConfigurator::configure(appender);
    log4cxx::PropertyConfigurator::configure(logConfigureDir + "/mtkpp-log.config");
    //log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());

    // Get logger object
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger(LOGGER_NAME));

    LOG4CXX_INFO(logger, "MTK Starting");
#endif
}

// ============================================================
// Function : ~errorHandler()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
errorHandler::~errorHandler() {}

// ============================================================
// Function : throwError()
// ------------------------------------------------------------
// 
// ============================================================
void errorHandler::throwError(error e)
{
    this->messages.push_back(e);
    if (this->maxMessages != 0 && (this->messages.size() > this->maxMessages)) {
      this->messages.pop_front();
    }

    if (this->log && (e.getType() <= this->level)) {
#ifdef USE_LOG4CXX
      // Get logger object
      log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger(LOGGER_NAME));

      std::stringstream ss;
      ss << e;
      switch(e.getType()) {
        case 0:
          LOG4CXX_INFO(logger, ss.str());
          break;
        case 1:
          LOG4CXX_INFO(logger, ss.str());
          break;
        case 2:
          LOG4CXX_WARN(logger, ss.str());
          break;
        case 3:
          LOG4CXX_DEBUG(logger, ss.str());
          break;
        case 4:
          LOG4CXX_ERROR(logger, ss.str());
          break;
        default:
          LOG4CXX_INFO(logger, ss.str());
          break;
      }
#else
      *stream << e;
      stream->flush();
#endif
    }
}

// ============================================================
// Function : throwError()
// ------------------------------------------------------------
// 
// ============================================================
void errorHandler::throwError(const std::string &function,
                              const std::string &message,
                              const int type)
{
    nMessages++;
    if (message.length() > 1) {
      error e(function, message, type);
      throwError(e);
    }
}

} // MTKpp namespace
