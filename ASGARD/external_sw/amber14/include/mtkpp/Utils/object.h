/*! 
   \file object.h
   \brief Base object
   \author Martin Peters

   $Date: 2009/04/08 11:02:50 $
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

#ifndef OBJECT_H
#define OBJECT_H

#include <boost/smart_ptr.hpp>

namespace MTKpp {
// ============================================================
// Class : object()
// ------------------------------------------------------------
/*! 
   \class object
   \brief Base object to enforce a virtual destructor
   \author Martin Peters
*/
// ============================================================
class object
{
public:
    object(void) {}
    virtual ~object(void) {}
};

typedef boost::shared_ptr< object > object_ptr;

}

#endif
