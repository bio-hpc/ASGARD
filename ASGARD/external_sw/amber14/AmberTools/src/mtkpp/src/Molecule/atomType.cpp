/*!
   \file atomType.cpp
   \brief Container of atom types
   \author Martin Peters

   $Date: 2007/09/14 11:35:24 $
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

#include "atomType.h"

namespace MTKpp
{

// ============================================================
// Function : atomTypes()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
atomTypes::atomTypes()
{
    itsName = "";
    pAtomType = 0;
}

// ============================================================
// Function : ~atomTypes()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
atomTypes::~atomTypes() {}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void atomTypes::setName(const std::string name)
{
    this->itsName = name;
}

// ============================================================
// Function : addatomType()
// ------------------------------------------------------------
// 
// ============================================================
atomTypeTMP* atomTypes::addAtomType()
{
    pAtomType = new atomTypeTMP();
    return pAtomType;
}

// ============================================================
// Function : setAtomTypeName()
// ------------------------------------------------------------
// 
// ============================================================
void atomTypes::setAtomTypeName(const std::string name)
{
    pAtomType->name = name;
    itsAtomTypeMap[name] = pAtomType;
}

// ============================================================
// Function : getAtomType()
// ------------------------------------------------------------
// 
// ============================================================
atomTypeTMP* atomTypes::getAtomType(const std::string atomType)
{
    atomTypeMapIterator e = itsAtomTypeMap.find(atomType);

    if (e != itsAtomTypeMap.end()){
         return e->second;
    }
    return NULL;
}

} // MTKpp namespace
