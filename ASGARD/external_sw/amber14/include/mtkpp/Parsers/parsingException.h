/*!
   \file parsingException.h
   \brief Parses mol files
   \author Roger Martin

   $Date: 2010/04/22 21:19:51 $
   $Revision: 1.1 $

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

#ifndef _PARSINGEXCEPTION_H
#define	_PARSINGEXCEPTION_H

#include <exception>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/exception.hpp>
#include <errno.h>
#include <iostream>

using namespace std;

namespace MTKpp
{
//public boost::exception,
class parsingException: public boost::exception {
public:
    string message;

    parsingException();
    parsingException(string message);
    parsingException(const parsingException& orig);
    virtual ~parsingException() throw();
public:
  virtual const char* what() const throw();
protected:

};

} // MTKpp namespace

#endif	/* _PARSINGEXCEPTION_H */

