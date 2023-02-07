/*!
   \file stdGroup.cpp
   \brief Container for standard fragments
   \author Martin Peters

   $Date: 2009/04/02 15:01:28 $
   $Revision: 1.11 $

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

#include "stdGroup.h"
#include "stdFrag.h"
#include "stdLibrary.h"

#include "molecule.h"
#include "submolecule.h"

namespace MTKpp
{

// ============================================================
// Function : stdGroup()
// ------------------------------------------------------------
// Constructor for the class.  
// ============================================================
stdGroup::stdGroup(stdLibrary *parent):pParent(parent)
{
    this->itsName = "";
    this->itsInfo = "";
    this->pStdFrag = 0;
    this->pMolecule = 0;
}

// ============================================================
// Function : ~stdGroup()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
stdGroup::~stdGroup()
{
    for (stdFragIterator c=itsStdFragList.begin(); c != itsStdFragList.end(); c++){
      this->pStdFrag = *c;
      delete this->pStdFrag;
    }
    this->itsStdFragList.clear();
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
// 
// ============================================================
stdLibrary* stdGroup::getParent()
{ 
    return this->pParent;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void stdGroup::setName(const std::string &name)
{
    this->itsName = name;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string stdGroup::getName()
{
    return this->itsName;
}

// ============================================================
// Function : setInfo()
// ------------------------------------------------------------
// 
// ============================================================
void stdGroup::setInfo(const std::string &i)
{
    this->itsInfo = i;
}

// ============================================================
// Function : getInfo()
// ------------------------------------------------------------
// 
// ============================================================
std::string stdGroup::getInfo()
{
    return this->itsInfo;
}

// ============================================================
// Function : addStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
stdFrag* stdGroup::addStdFrag()
{
    this->pStdFrag = new stdFrag(this);
    this->itsStdFragList.push_back(pStdFrag);
    return this->pStdFrag; 
}

// ============================================================
// Function : addStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
stdFrag* stdGroup::addStdFrag(stdFrag* f)
{
    if (f) {
      this->pStdFrag = new stdFrag(f, this);
      this->itsStdFragList.push_back(pStdFrag);
      return this->pStdFrag; 
    }
    return 0;
}

// ============================================================
// Function : appendStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
void stdGroup::appendStdFrag(stdFrag* f)
{
    if (f) {
      this->pStdFrag = new stdFrag(f, this);
      this->itsStdFragList.push_back(pStdFrag);
    }
}

// ============================================================
// Function : getStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
stdFrag* stdGroup::getStdFrag(const std::string &name)
{
    for (stdFragIterator c = itsStdFragList.begin(); c != itsStdFragList.end(); c++) {
      this->pStdFrag = *c;
      if (pStdFrag->getSymbol() == name) {
        return this->pStdFrag;
      }
    }
    return 0;
}

// ============================================================
// Function : hasStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
bool stdGroup::hasStdFrag(const std::string &name)
{
    for (stdFragIterator c = itsStdFragList.begin(); c != itsStdFragList.end(); c++) {
      this->pStdFrag = *c;
      if (this->pStdFrag->getSymbol() == name) {
        return true;
      }
    }
    return false;
}

// ============================================================
// Function : getStdFragList()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<stdFrag*> stdGroup::getStdFragList()
{
    return this->itsStdFragList;
}

// ============================================================
// Function : hasStdMolecule()
// ------------------------------------------------------------
// 
// ============================================================
bool stdGroup::hasStdMolecule()
{
    if (pMolecule) {
      return true;
    }
    return false;
}

// ============================================================
// Function : setStdMolecule()
// ------------------------------------------------------------
// 
// ============================================================
void stdGroup::setStdMolecule(molecule* m)
{
    pMolecule = m;
}

// ============================================================
// Function : getStdMolecule()
// ------------------------------------------------------------
// 
// ============================================================
molecule* stdGroup::getStdMolecule()
{
    return pMolecule;
}

// ============================================================
// Function : generateSimpleFP()
// ------------------------------------------------------------
// 
// ============================================================
void stdGroup::generateSimpleFP()
{
    for (stdFragIterator c=itsStdFragList.begin(); c != itsStdFragList.end(); c++){
      this->pStdFrag = *c;
      this->pStdFrag->generateSimpleFP();
    }
}

// ============================================================
// Function : generateAdjMatrices()
// ------------------------------------------------------------
// 
// ============================================================
void stdGroup::generateAdjMatrices()
{
    for (stdFragIterator c = itsStdFragList.begin(); c != itsStdFragList.end(); c++) {
      this->pStdFrag = *c;
      this->pStdFrag->generateAdjMatrix();
    }
}

// ============================================================
// Function : getCharge()
// ------------------------------------------------------------
// 
// ============================================================
double stdGroup::getCharge()
{
    double d = 0.0;
    if (pMolecule) {
      std::vector<submolecule*> sList = this->pMolecule->getSubMoleculeList();
      for (unsigned int p = 0; p < sList.size(); p++) {
        this->pStdFrag = this->getStdFrag(sList[p]->getName());
        d += this->pStdFrag->getCharge();
        //std::cout << this->pStdFrag->getSymbol() << " " << this->pStdFrag->getCharge() << std::endl;
      }
    }
    return d;
}

} // MTKpp namespace

