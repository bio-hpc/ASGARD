/*!
   \file stdLibrary.cpp
   \brief Container for standard groups
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
   $Revision: 1.12 $

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

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"

#include "Log/errorHandler.h"

namespace MTKpp
{

stdLibrary* StdLibrary=NULL;

// ============================================================
// Function : getInstance()
// ------------------------------------------------------------
//
// ============================================================
stdLibrary* stdLibrary::getInstance()
{
    if (StdLibrary == NULL)
    {
        StdLibrary = new stdLibrary();
    }
    return StdLibrary;
}

// ============================================================
// Function : stdLibrary()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
stdLibrary::stdLibrary()
{
    itsName = "";
    pStdGroup = 0;
    pStdFrag = 0;

    bSimpleFPGenerated = 0;
    bAdjMatricesGenerated = 0;
    bAtomKindsGenerated = 0;

    LLL2L["ALA"] = "A";
    LLL2L["GLY"] = "G";
    LLL2L["VAL"] = "V";
    LLL2L["LEU"] = "L";
    LLL2L["ILE"] = "I";
    LLL2L["PHE"] = "F";
    LLL2L["TYR"] = "Y";
    LLL2L["TRP"] = "W";
    LLL2L["SER"] = "S";
    LLL2L["CYS"] = "C";
    LLL2L["CYM"] = "C";
    LLL2L["CYX"] = "C";
    LLL2L["THR"] = "T";
    LLL2L["MET"] = "M";
    LLL2L["PRO"] = "P";
    LLL2L["HIS"] = "H";
    LLL2L["HID"] = "H"; // HIS with Hydrogen on delta nitrogen
    LLL2L["HIE"] = "H"; // HIS with Hydrogen on epsilon nitrogen
    LLL2L["HIP"] = "H"; // HIS with Hydrogen on both nitrogens, +ve charge
    LLL2L["LYS"] = "K";
    LLL2L["ARG"] = "R";
    LLL2L["ASP"] = "D";
    LLL2L["GLU"] = "E";
    LLL2L["GLN"] = "Q";
    LLL2L["ASN"] = "N";
    LLL2L["HOH"] = "O";
    LLL2L["WAT"] = "O";
}

// ============================================================
// Function : ~stdLibrary()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
stdLibrary::~stdLibrary()
{
    for (stdGroupIterator c = itsStdGroupList.begin(); c != itsStdGroupList.end(); c++) {
      pStdGroup = *c;
      delete pStdGroup;
    }
    itsStdGroupList.clear();
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void stdLibrary::setName(const std::string &name)
{
    itsName = name;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string stdLibrary::getName()
{
    return itsName;
}

// ============================================================
// Function : addStdGroup()
// ------------------------------------------------------------
// 
// ============================================================
stdGroup* stdLibrary::addStdGroup()
{
    pStdGroup = new stdGroup(this);
    this->itsStdGroupList.push_back(pStdGroup);
    return pStdGroup; 
}

// ============================================================
// Function : getStdGroup()
// ------------------------------------------------------------
// 
// ============================================================
stdGroup* stdLibrary::getStdGroup(const std::string &name)
{
    for (stdGroupIterator c = itsStdGroupList.begin();
         c != itsStdGroupList.end(); c++) {
      pStdGroup = *c;
      if (pStdGroup->getName() == name) {
        return pStdGroup;
      }
    }
    return 0;
}

// ============================================================
// Function : getStdGroupList()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<stdGroup*> stdLibrary::getStdGroupList()
{
    return this->itsStdGroupList;
}

// ============================================================
// Function : getStdFrag()
// ------------------------------------------------------------
// Return first stdFrag with name equals fragName from the 
// entire standard library
// ============================================================
stdFrag* stdLibrary::getStdFrag(const std::string &fragName)
{
    for (unsigned int i = 0; i < itsStdGroupList.size(); i++) {
      pStdGroup = itsStdGroupList[i];
      if (pStdGroup) {
        std::vector<stdFrag*> fragList = pStdGroup->getStdFragList();
        for (unsigned int j = 0; j < fragList.size(); j++) {
          pStdFrag = fragList[j];
          if (pStdFrag->getSymbol() == fragName) {
            return pStdFrag;
          }
        }
      }
      else {
        return 0;
      }
    }
    return 0;
}

// ============================================================
// Function : getStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
stdFrag* stdLibrary::getStdFrag(const std::string &fragName, const std::string &n)
{
    // First check if n is a group
    pStdGroup = this->getStdGroup(n);
    if (pStdGroup) {
      pStdFrag = pStdGroup->getStdFrag(fragName);
      if (pStdFrag) {
        return pStdFrag;
      }
   }
   else {
     for (unsigned int i = 0; i < this->itsStdGroupList.size(); i++) {
       pStdGroup = this->itsStdGroupList[i];
       if (pStdGroup) {
         std::vector<stdFrag*> fragList = pStdGroup->getStdFragList();
         for (unsigned int j = 0; j < fragList.size(); j++) {
           pStdFrag = fragList[j];
           if ((pStdFrag->getSymbol() == fragName) and (pStdFrag->getType() == n)) {
             return pStdFrag;
           }
         }
       }
       else {
         return 0;
       }
     }
   }
   return 0;
}

// ============================================================
// Function : getNumberStdFrag()
// ------------------------------------------------------------
// 
// ============================================================
int stdLibrary::getNumberStdFrag()
{
    int n = 0;
    for (unsigned int i = 0; i < itsStdGroupList.size(); i++) {
      pStdGroup = itsStdGroupList[i];
      if (pStdGroup) {
        std::vector<stdFrag*> fragList = pStdGroup->getStdFragList();
        n += fragList.size();
      }
      else {
      }
    }
    return n;
}

// ============================================================
// Function : setL()
// ------------------------------------------------------------
// 
// ============================================================
void stdLibrary::setL(std::string lll, std::string l)
{
    this->LLL2L[lll] = l;
}


// ============================================================
// Function : getL()
// ------------------------------------------------------------
// 
// ============================================================
std::string stdLibrary::getL(std::string s)
{
    typedef std::map<std::string, std::string>::iterator strStrMapIterator;
    strStrMapIterator b = LLL2L.find(s);
    if (b != LLL2L.end()) {
      return LLL2L[s];
    }
    std::string errMessage = "Can't find 1 Letter code for " + s;
    errorLogger.throwError("stdLibrary::getL", errMessage, WARNING);
    //throw MTKException(errMessage);
    return "X";
}

// ============================================================
// Function : generateSimpleFP()
// ------------------------------------------------------------
// 
// ============================================================
void stdLibrary::generateSimpleFP()
{
    for (stdGroupIterator c = itsStdGroupList.begin(); c != itsStdGroupList.end(); c++) {
      pStdGroup = *c;
      pStdGroup->generateSimpleFP();
    }
    bSimpleFPGenerated = 1;
}

// ============================================================
// Function : generateAdjMatrices()
// ------------------------------------------------------------
// 
// ============================================================
void stdLibrary::generateAdjMatrices()
{
    for (stdGroupIterator c = itsStdGroupList.begin(); c != itsStdGroupList.end(); c++) {
      pStdGroup = *c;
      pStdGroup->generateAdjMatrices();
    }
    bAdjMatricesGenerated = 1;
}

// ============================================================
// Function : generateAtomKinds()
// ------------------------------------------------------------
// 
// ============================================================
void stdLibrary::generateAtomKinds()
{
    for (stdGroupIterator c = itsStdGroupList.begin(); c != itsStdGroupList.end(); c++) {
      pStdGroup = *c;
      if (pStdGroup) {
        std::vector<stdFrag*> fragList = pStdGroup->getStdFragList();
        for (unsigned int j = 0; j < fragList.size(); j++) {
          pStdFrag = fragList[j];
          pStdFrag->generateAtomKinds();
        }
      }
    }
    bAtomKindsGenerated = 1;
}

} // MTKpp namespace

