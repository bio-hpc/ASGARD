/*!
   \file angle.h
   \brief Container for angle information
   \author Martin Peters

   $Date: 2007/09/14 11:35:23 $
   $Revision: 1.4 $

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

#ifndef ANGLE_h
#define ANGLE_h

#include <iostream>
#include <string>

namespace MTKpp
{

class atom;
struct angleParam;

// ============================================================
// Struct : Angle
// ------------------------------------------------------------
/*! 
   \struct Angle
   \brief Container for Angle info
   \author Martin Peters
   \date 2005
*/
// ============================================================
struct Angle
{
    //! First Atom in Angle
    atom* atom1;

    //! Second Atom in Angle
    atom* atom2;

    //! Third Atom in Angle
    atom* atom3;

    //! Size of Angle
    double size;

    //! Pointer to Angle Parameters
    angleParam* pAngleParam;
};

} // MTKpp namespace

#endif // ANGLE_h

