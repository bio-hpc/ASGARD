/*!
   \file edge.h
   \brief Container for edge information
   \author Martin Peters

   $Date: 2007/09/14 07:19:22 $
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

#ifndef EDGE_H
#define EDGE_H

#include <iostream>
#include <string>

namespace MTKpp
{

class vertex;

// ============================================================
// Struct : edge
// ------------------------------------------------------------
/*!
   \struct edge
   \brief Container for edge info
   \author Martin Peters
   \date 2006
*/
// ============================================================
struct edge
{
    //! First vertex in edge
    vertex* v1;

    //! Second vertex in edge
    vertex* v2;

    //! ?
    bool visited;

    //! value
    double itsValue;
};

} // MKT++ namespace

#endif // EDGE_H
