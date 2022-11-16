/*!
   \file utility.h
   \brief Contains string functions
   \author Martin Peters

   $Date: 2010/03/29 20:46:06 $
   $Revision: 1.9 $

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

#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <vector>
#include <algorithm>

namespace MTKpp
{

// ============================================================
// Function : i2s()
// ------------------------------------------------------------
//
// ============================================================
inline std::string i2s(int i)
{
    std::stringstream number;
    number << i;
    return number.str();
}

// ============================================================
// Function : d2s()
// ------------------------------------------------------------
//
// ============================================================
inline std::string d2s(double d)
{
    std::stringstream number;
    number << d;
    return number.str();
}

} // MTKpp namespace

#endif // UTILITY_H

