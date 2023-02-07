/*! 
   \file error.h
   \brief Error object used in errorHandler
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

#ifndef ERROR_H
#define ERROR_H

#include <iostream>
#include <string>
#include "Utils/constants.h"

namespace MTKpp
{

// ============================================================
// Class : error()
// ------------------------------------------------------------
/*! 
   \class error
   \brief Error object used in errorHandler
   \author Martin Peters
   \version 0.1
   \date 2007
*/
// ============================================================

class error
{
public:

    /*!
       \brief error Constructor
       \param function Function where message originates
       \param message Error/warning message
       \param type Error type
    */
    error(const std::string &function = "",
          const std::string &message = "",
          const int type = 4);

    //! error Destructor
    virtual ~error();
    
    /*!
       \brief Get formatted error message
    */
    std::string getFormattedMessage(void) const;

    /*!
       \brief Output a formatted message string
    */
    friend std::ostream& operator<< ( std::ostream &os, const error &e )
      { return os << e.getFormattedMessage(); };

    /*!
       \brief Get associate function for error
    */
    std::string getFunction() { return function; }

    /*!
       \brief Get message for error
    */
    std::string getMessage() { return message; }

    /*!
       \brief Get type of error
    */
    int getType() { return type; }

protected:

    //! function where error was called
    std::string function;

    //! message programmer used
    std::string message;

    /*
       \brief error type
       1 - Error
       2 - Warning
       3 - Debug
       4 - Info
    */
    int type;
};

} // MTKpp namespace

#endif // ERROR_H
