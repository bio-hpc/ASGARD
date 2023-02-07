/*! 
   \file errorHandler.h
   \brief Error Handling within MTK++
   \author Martin Peters

   $Date: 2010/03/29 20:26:54 $
   $Revision: 1.6 $

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

#ifndef ERRORHANDLER_H
#define ERRORHANDLER_H

#include <iostream>
#include <string>
#include <deque>

#include "Utils/constants.h"

namespace MTKpp
{

class error;

// ============================================================
// Class : errorHandler()
// ------------------------------------------------------------
/*! 
   \class errorHandler
   \brief Error Handling within MTK++
   \author Martin Peters
   \version 0.1
   \date 2007
*/
// ============================================================

class errorHandler
{
public:

    //! errorHandler Constructor
    errorHandler();

    //! errorHandler Destructor
    virtual ~errorHandler();

    /*!
       \brief Throw error message
       \param function Function where message originates
       \param message Error/warning message
       \param type error type
    */
    void throwError(const std::string &function, const std::string &message,
                    int type = 4);

    /*!
       \brief Set the output stream
       \param os Stream
    */
    void setStream(std::ostream *os) { stream = os; }

    /*!
       \brief Get the output stream
    */
    std::ostream* getStream() { return stream; }

    /*!
       \brief Flush output
    */
    void flush() { stream->flush(); }

    /*!
       \brief Turn logging on/off
       \param l On/off
    */
    void setLogging(bool l) { log = l; }

    /*!
       \brief Set level of output
       \param l Level
    */
    void setLevel(int l) { level = l; }

    /*!
       \brief Set maximum number of messages
       \param m maximum number of messages
    */
    void setmaxMessages(int m) { this->maxMessages = m; }

    /*!
       \brief Get maximum number of messages
    */
    unsigned int getmaxMessages() { return this->maxMessages; }

    /*!
       \brief Get total number of messages
    */
    unsigned int getnMessages() { return this->nMessages; }

protected:

    /*!
       \brief Throw error message from errorHandler
       \param err error 
    */
    void throwError(error err);

    //! Level of output
    int level;

    //! Output stream
    std::ostream *stream;

    //! Write output boolean
    bool log;

    //! Stored for a total of maxMessages messages
    std::deque<error> messages;

    //! Total number of messages created
    unsigned int nMessages;

    //! Maximum allowed number of message to store
    unsigned int maxMessages;
};

//! error logger global variable
extern errorHandler errorLogger;

} // MTKpp namespace

#endif // ERRORHANDLER_H
