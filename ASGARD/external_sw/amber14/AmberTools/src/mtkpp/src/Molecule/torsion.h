/*!
   \file torsion.h
   \brief Container for Torsion info
   \author Martin Peters

   $Date: 2008/02/28 14:06:06 $
   $Revision: 1.5 $

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

#ifndef TORSION_h
#define TORSION_h

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace MTKpp
{

class atom;
struct torsionParam;

// ============================================================
// Struct : Torsion
// ------------------------------------------------------------
/*! 
   \struct Torsion
   \brief Container for Torsion info
   \author Martin Peters
   \date 2005
*/
// ============================================================
struct Torsion
{
    //! First Atom in Torsion
    atom* atom1;

    //! Second Atom in Torsion
    atom* atom2;

    //! Third Atom in Torsion
    atom* atom3;

    //! Fourth Atom in Torsion
    atom* atom4;

    //! Size of Torsion
    double size;

    //! Parameters assigned
    bool parametersAssigned;

    //! Vector of Torsion Parameters
    std::vector<torsionParam*> pTorsionParamList;
};

} // MTKpp namespace

#endif // TORSION_h
