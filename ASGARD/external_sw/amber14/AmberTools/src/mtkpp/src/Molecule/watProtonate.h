/*!
   \file watProtonate.h
   \brief Protonates water molecule
   \author Martin Peters

   $Date: 2010/03/29 20:46:06 $
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

#ifndef WATPROTONATE_H
#define WATPROTONATE_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include <boost/array.hpp>

namespace MTKpp
{

class collection;
class atom;
class vector3d;

// ============================================================
// Class : watProtonate()
// ------------------------------------------------------------
/*!
   \class watProtonate

   \brief Class to add hydrogens to water molecules

   \author Martin Peters

   \version 0.1

   \date 2005
*/
// ============================================================
class watProtonate
{
public:

    /*!
       \brief watProtonate Constructor
    */
    watProtonate();

    //! watProtonate Destructor
    virtual ~watProtonate();

    /*!
       \brief Add hydrogen atoms to a water molecules
    */
    void run(collection* pCol);

    /*!
       \brief 
    */
    void setup();

protected: // DATA

    //!
    int nSphereGridPts;

    //!
    int biggestNumberOfHits;

    //!
    int *sphereGridCoords;

    //!
    int *maxPairs;

    //!
    int *Hs;

    //!
    int *lonePairs;
};

} // MTKpp namespace

#endif // WATPROTONATE_H
