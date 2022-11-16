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

#include "parsingException.h"

namespace MTKpp
{

parsingException::parsingException() {
}

parsingException::parsingException(string message) {
    this->message=message;
}

parsingException::parsingException(const parsingException& orig) {
}

parsingException::~parsingException() throw() {
}

const char* parsingException::what() const throw()
{
return message.c_str();
}

}
