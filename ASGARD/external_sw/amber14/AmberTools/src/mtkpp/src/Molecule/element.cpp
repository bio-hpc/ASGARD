/*!
   \file element.cpp
   \brief Container for element information
   \author Martin Peters

   $Date: 2008/02/28 14:06:06 $
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


#include "element.h"

namespace MTKpp
{

// ============================================================
// Function : element()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
elements::elements()
{
    pElement = 0;
}

// ============================================================
// Function : ~element()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
elements::~elements() {}

// ============================================================
// Function : addElement()
// ------------------------------------------------------------
// 
// ============================================================
element* elements::addElement()
{
    pElement = new element();

    pElement->number = -1;
    pElement->symbol = "";
    pElement->name = "";
    pElement->mass = -1.0;
    pElement->group = -1;
    pElement->period = -1;
    pElement->red = -1.0;
    pElement->green = -1.0;
    pElement->blue = -1.0;
    pElement->valence = -1;
    pElement->filledShell = -1;
    pElement->covalentRadius = -1.0;
    pElement->vdWRadius = -1.0;
    pElement->paulingEN = -1.0;

    return pElement;
}

// ============================================================
// Function : setElementName()
// ------------------------------------------------------------
// 
// ============================================================
void elements::setElementName(const std::string symbol)
{
    pElement->symbol = symbol;
    this->itsElementMap[symbol] = pElement;
}

// ============================================================
// Function : getElement()
// ------------------------------------------------------------
// 
// ============================================================
element* elements::getElement(const std::string element)
{
    ElementMapIterator e = this->itsElementMap.find(element);

    if (e != this->itsElementMap.end()) {
      return e->second;
    }
    return 0;
}

// ============================================================
// Function : getElementNumber()
// ------------------------------------------------------------
// Returns atomic number
// ============================================================
int elements::getElementNumber(const std::string element)
{
    ElementMapIterator e = this->itsElementMap.find(element);

    if (e != this->itsElementMap.end()) {
      return e->second->number;
    }
    return -1;
}

// ============================================================
// Function : getElementName()
// ------------------------------------------------------------
// Returns element name
// ============================================================
std::string elements::getElementName(const std::string element)
{
    ElementMapIterator e = this->itsElementMap.find(element);

    if (e != itsElementMap.end()) {
      return e->second->name;
    }
    return 0;
}

// ============================================================
// Function : getElementMass()
// ------------------------------------------------------------
// Returns atomic mass
// ============================================================
double elements::getElementMass(const std::string element)
{
    ElementMapIterator e = this->itsElementMap.find(element);

    if (e != this->itsElementMap.end()) {
      return e->second->mass;
    }
    return -1;
}

// ============================================================
// Function : hasSEHamiltonian()
// ------------------------------------------------------------
// Returns atomic mass
// ============================================================
bool elements::hasSEHamiltonian(const std::string element, const std::string ham)
{
    ElementMapIterator e = this->itsElementMap.find(element);

    if (e != this->itsElementMap.end()) {
      if (e->second->seHams.size() == 0) return false;
      strIterator a = std::find(e->second->seHams.begin(), e->second->seHams.end(), ham);
      if (a == e->second->seHams.end()) {
        //std::cout << " elements::hasSEHamiltonian " << element << " not in " << ham << std::endl;
        return false;
      }
    }
    return true;
}

} // MTKpp namespace

