/*!
   \file atomType.h
   \brief Container of atom types
   \author Martin Peters

   $Date: 2010/03/29 20:42:27 $
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

#ifndef ATOMTYPE_H
#define ATOMTYPE_H


#include <iostream>
#include <string>
#include <map>
#include <algorithm>

namespace MTKpp
{

// ============================================================
// Struct : atomType()
// ------------------------------------------------------------
/*! 
   \struct atomType
   \brief 
   \author Martin Peters
   \date 2005
*/
// ============================================================
struct atomTypeTMP
{
    //! 
    std::string name;

    //!
    std::string element;

    //! 
    double mass;

    //!
    std::string hybridization;

    //!
    std::string description;

    //!
    double rvalue;

    //!
    double evalue;

    //! Atomic polarizability in A**3
    double atomPolarizability;
};

// ============================================================
// Class : atomTypes()
// ------------------------------------------------------------
/*! 
   \class atomTypes
   \brief 
   \author Martin Peters
   \date 2005
*/
// ============================================================
class atomTypes
{
public:
     atomTypes();
     virtual ~atomTypes();

     void      setName(const std::string);

     atomTypeTMP* addAtomType();
     void      setAtomTypeName(const std::string);
     atomTypeTMP* getAtomType(const std::string);

protected:
    std::string itsName;
    // - MAP CONTAINER - //
    std::map<std::string, atomTypeTMP*> itsAtomTypeMap;

    // - MAP ITERATORS - //
    typedef std::map<std::string, atomTypeTMP*>::iterator atomTypeMapIterator;

    // - POINTERS - //
    atomTypeTMP* pAtomType;
};

} // MTKpp namespace

#endif // ATOMTYPE_H
