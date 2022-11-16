/*!
   \file atom.cpp
   \brief Container for atom information
   \author Martin Peters

   $Date: 2010/03/29 20:42:27 $
   $Revision: 1.22 $

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

#include "atom.h"
#include "element.h"
#include "Utils/vector3d.h"
#include "stdFrag.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : atom(atom*, submolecule*)
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
atom::atom(atom* pAt, submolecule *parent):pParent(parent)
{
    this->pElement = pAt->getElement();
    if (!this->pElement) {
      std::cout << " Adding atom with unknown element " << pAt->getName() << " ... exiting " <<  std::endl;
      //exit(0);
      throw MTKException();
    }
    this->pCoords = new vector3d(0.0);
    (*this->pCoords)[0] = (*pAt->getCoords())[0];
    (*this->pCoords)[1] = (*pAt->getCoords())[1];
    (*this->pCoords)[2] = (*pAt->getCoords())[2];

    //pCoords = pAt->getCoords();
    itsElement = pAt->itsElement;
    itsName = pAt->getName();
    itsFileId = pAt->getFileID();
    itsValence = 0;
    itsFormalCharge = 0;
    itsHybridization = pAt->getHybridization();
    itsHydrophilicity = pAt->getHydrophilicity();
    itsType = 0;
    pStdAtom = 0;
    pStdAtom = pAt->getStdAtom();
    //itsIndex = pAt->getIndex();
    //itsCollectionIndex = 0;
}

// ============================================================
// Function : atom(submolecule*)
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
atom::atom(submolecule *parent):pParent(parent)
{
    this->pElement = 0;
    this->pCoords = new vector3d(0.0);
    this->itsElement = "";
    this->itsName = "";
    this->itsIndex = 0;
    this->itsCollectionIndex = 0;
    this->itsFileId = 0;
    this->itsAtomicNum = 0;
    this->itsZcharge = 0.0;
    this->itsValence = 0;
    this->itsFormalCharge = 0;
    this->itsHybridization = 0;
    this->itsHydrophilicity = -1;
    this->itsType = 0;
    this->pStdAtom = 0;
}

// ============================================================
// Function : ~atom()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
atom::~atom()
{
    delete pCoords;
}

// ============================================================
// Function : setElement()
// ------------------------------------------------------------
// Sets the Element Symbol for this atom.
// ============================================================
void atom::setElement(element* ele)
{
    this->pElement = ele;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// Sets the name for this atom.
// ============================================================
void atom::setName(std::string name)
{
    this->itsName = name;
}

// ============================================================
// Function : setCoords()
// ------------------------------------------------------------
// Sets the Coordinates for this atom.
// ============================================================
void atom::setCoords(const double &x, const double &y, const double &z)
{
    this->pCoords->set(x,y,z);
}

// ============================================================
// Function : setIndex()
// ------------------------------------------------------------
// Sets the molecule index for this atom
// ============================================================
void atom::setIndex(const int &n)
{
    this->itsIndex = n;
}

// ============================================================
// Function : setColIndex()
// ------------------------------------------------------------
// Sets the collection index for this atom.
// ============================================================
void atom::setColIndex(const int &n)
{
    this->itsCollectionIndex = n;
}

// ============================================================
// Function : setFileID()
// ------------------------------------------------------------
// Sets the Coordinates for this atom.
// ============================================================
void atom::setFileID(const int &n)
{
    this->itsFileId = n;
}

// ============================================================
// Function : setAtomicNum()
// ------------------------------------------------------------
//
// ============================================================
void atom::setAtomicNum(const int &n)
{
    this->itsAtomicNum = n;
}

// ============================================================
// Function : setZcharge()
// ------------------------------------------------------------
//
// ============================================================
void atom::setZcharge(const double &n)
{
    this->itsZcharge = n;
}

// ============================================================
// Function : setFormalCharge()
// ------------------------------------------------------------
//
// ============================================================
void atom::setFormalCharge(const int &n)
{
    this->itsFormalCharge = n;
}

// ============================================================
// Function : setValence()
// ------------------------------------------------------------
//
// ============================================================
void atom::setValence(const int &n)
{
    this->itsValence = n;
}

// ============================================================
// Function : setHybridization()
// ------------------------------------------------------------
//
// ============================================================
void atom::setHybridization(const int &i)
{
    this->itsHybridization = i;
}

// ============================================================
// Function : setType()
// ------------------------------------------------------------
//
// ============================================================
void atom::setType(const int &n)
{
    this->itsType = n;
}

// ============================================================
// Function : setMengType()
// ------------------------------------------------------------
//
// ============================================================
void atom::setMengType(const std::string &n)
{
    this->itsMengType = n;
}

// ============================================================
// Function : setOccupancy()
// ------------------------------------------------------------
//
// ============================================================
void atom::setOccupancy(const double &o)
{
    this->itsOccupancy = o;
}

// ============================================================
// Function : setTempFactor()
// ------------------------------------------------------------
//
// ============================================================
void atom::setTempFactor(const double &b)
{
    this->itsTempFactor = b;
}

// ============================================================
// Function : setHydrophilicity()
// ------------------------------------------------------------
//
// ============================================================
void atom::setHydrophilicity(const int &h)
{
    this->itsHydrophilicity = h;
}

// ============================================================
// Function : getOccupancy()
// ------------------------------------------------------------
//
// ============================================================
double atom::getOccupancy()
{
    return this->itsOccupancy;
}

// ============================================================
// Function : getTempFactor()
// ------------------------------------------------------------
//
// ============================================================
double atom::getTempFactor()
{
    return this->itsTempFactor;
}

// ============================================================
// Function : getHydrophilicity()
// ------------------------------------------------------------
//
// ============================================================
int atom::getHydrophilicity()
{
    return this->itsHydrophilicity;
}

// ============================================================
// Function : getMengType()
// ------------------------------------------------------------
//
// ============================================================
std::string atom::getMengType()
{
    return this->itsMengType;
}

// ============================================================
// Function : addProperty()
// ------------------------------------------------------------
// Add property
// ============================================================
void atom::addProperty(const std::string &name, double value)
{
    this->itsPropertiesMap[name] = value;
}

// ============================================================
// Function : addProperty()
// ------------------------------------------------------------
// Add property
// ============================================================
void atom::addProperty(const std::string &name, int value)
{
    this->itsIntPropertiesMap[name] = value;
}

// ============================================================
// Function : hasProperty()
// ------------------------------------------------------------
// Has atom got a certain property
// ============================================================
bool atom::hasProperty(const std::string &name)
{
    PropertyMapIterator p1 = this->itsPropertiesMap.find(name);
    if (p1 != this->itsPropertiesMap.end()) {
      return true;
    }

    intPropertyMapIterator p2 = this->itsIntPropertiesMap.find(name);
    if (p2 != this->itsIntPropertiesMap.end()) {
      return true;
    }
    return false;
}

// ============================================================
// Function : getProperty()
// ------------------------------------------------------------
// Get atom property by name
// ============================================================
double atom::getProperty(const std::string &name)
{
    PropertyMapIterator p = this->itsPropertiesMap.find(name);
    if (p != this->itsPropertiesMap.end()){
      return p->second;
    }
    return 0.0;
}

// ============================================================
// Function : getPropertyMap()
// ------------------------------------------------------------
// Get atom property map
// ============================================================
std::map<std::string, double> atom::getPropertyMap()
{
    return this->itsPropertiesMap;
}

// ============================================================
// Function : getIntProperty()
// ------------------------------------------------------------
// Get atom integer property
// ============================================================
int atom::getIntProperty(const std::string &name)
{
    intPropertyMapIterator p = this->itsIntPropertiesMap.find(name);
    if (p != this->itsIntPropertiesMap.end()){
      return p->second;
    }
    return 0;
}

// ============================================================
// Function : getIntPropertyMap()
// ------------------------------------------------------------
// Get atom integer property map
// ============================================================
std::map<std::string, int> atom::getIntPropertyMap()
{
    return this->itsIntPropertiesMap;
}

// ============================================================
// Function : getElement()
// ------------------------------------------------------------
// Gets the Element for this atom.
// ============================================================
element* atom::getElement()
{
    if (this->pElement == 0) {
      return 0;
    }
    return this->pElement;
}

// ============================================================
// Function : getElementSymbol()
// ------------------------------------------------------------
// Gets the Element Symbol for this atom.
// ============================================================
std::string atom::getElementSymbol()
{
    if (this->pElement == 0) {
      return "XX";
    }
    return this->pElement->symbol;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// Gets the name for this atom.
// ============================================================
std::string atom::getName()
{
    return this->itsName;
}

// ============================================================
// Function : getCoords
// ------------------------------------------------------------
// returns the atom's cartesian coordinates
// ============================================================
vector3d* atom::getCoords()
{
    return this->pCoords;
}

// ============================================================
// Function : getIndex()
// ------------------------------------------------------------
// Gets the internal molecule index of the atom
// ============================================================
int atom::getIndex()
{
    return this->itsIndex;
}

// ============================================================
// Function : getColIndex()
// ------------------------------------------------------------
// Gets the internal collection index of the atom
// ============================================================
int atom::getColIndex()
{
    return this->itsCollectionIndex;
}

// ============================================================
// Function : getFileID()
// ------------------------------------------------------------
// Gets the input file index of the atom
// ============================================================
int atom::getFileID()
{
    return this->itsFileId;
}

// ============================================================
// Function : getX
// ------------------------------------------------------------
// return X coordinate
// ============================================================
double atom::getX()
{
    return this->pCoords->getX();
}

// ============================================================
// Function : getY
// ------------------------------------------------------------
// return Y
// ============================================================
double atom::getY()
{
    return this->pCoords->getY();
}

// ============================================================
// Function : getZ
// ------------------------------------------------------------
// return Z
// ============================================================
double atom::getZ()
{
    return this->pCoords->getZ();
}

// ============================================================
// Function : getParent
// ------------------------------------------------------------
// return submolecule
// ============================================================
submolecule* atom::getParent()
{
    return this->pParent;
}

// ============================================================
// Function : getAtomicNum()
// ------------------------------------------------------------
//
// ============================================================
int atom::getAtomicNum()
{
    if (this->pElement == 0) {
      return -1;
    }
    return this->pElement->number;
}

// ============================================================
// Function : getAtomicMass()
// ------------------------------------------------------------
//
// ============================================================
double atom::getAtomicMass()
{
    if (this->pElement == 0) {
      return 0.0;
    }
    return this->pElement->mass;
}

// ============================================================
// Function : getZcharge()
// ------------------------------------------------------------
//
// ============================================================
double atom::getZcharge()
{
    return this->itsZcharge;
}

// ============================================================
// Function : getFormalCharge()
// ------------------------------------------------------------
//
// ============================================================
int atom::getFormalCharge()
{
    return this->itsFormalCharge;
}

// ============================================================
// Function : getValence()
// ------------------------------------------------------------
//
// ============================================================
int atom::getValence()
{
    return this->itsValence;
}

// ============================================================
// Function : getHybridization()
// ------------------------------------------------------------
//
// ============================================================
int atom::getHybridization()
{
     return this->itsHybridization;
}

// ============================================================
// Function : getType()
// ------------------------------------------------------------
//
// ============================================================
int atom::getType()
{
    return this->itsType;
}

// ============================================================
// Function : numBondedOxygens()
// ------------------------------------------------------------
//
// ============================================================
int atom::numBondedOxygens()
{
    int nOxygens = 0;
    for (unsigned int i = 0; i < this->bondedAtoms.size(); i++) {
      if (this->bondedAtoms[i]->getElementSymbol() != "O") {
        nOxygens++;
      }
    }
    return nOxygens;
}

// ================================================ //
// ===                                          === //
// ===  M o l e c u l a r  M e c h a n i c s    === //
// ===                                          === //
// ================================================ //

// ============================================================
// Function : setStdAtom()
// ------------------------------------------------------------
//
// ============================================================
void atom::setStdAtom(stdAtom* sa)
{
    this->pStdAtom = sa;
}

// ============================================================
// Function : getStdAtom()
// ------------------------------------------------------------
//
// ============================================================
stdAtom* atom::getStdAtom()
{
    if (this->pStdAtom) {
      return this->pStdAtom;
    }
    return 0;
}

// ============================================================
// Function : addBondedAtom()
// ------------------------------------------------------------
//
// ============================================================
void atom::addBondedAtom(atom* a)
{
    std::vector<atom*>::iterator result;
    result = std::find(this->bondedAtoms.begin(), this->bondedAtoms.end(), a);
    if (result ==  this->bondedAtoms.end()) {
      this->bondedAtoms.push_back(a);
    }
}

// ============================================================
// Function : addBonded13Atom()
// ------------------------------------------------------------
//
// ============================================================
bool atom::addBonded13Atom(atom* a)
{
    std::vector<atom*>::iterator result;
    result = std::find(this->bonded13Atoms.begin(), this->bonded13Atoms.end(), a);
    if (result == this->bonded13Atoms.end()) {
      this->bonded13Atoms.push_back(a);
      return true;
    }
    return false;
}

// ============================================================
// Function : addBonded14Atom()
// ------------------------------------------------------------
//
// ============================================================
void atom::addBonded14Atom(atom* a)
{
    std::vector<atom*>::iterator result;
    result = std::find(this->bonded14Atoms.begin(), this->bonded14Atoms.end(), a);
    if (result ==  this->bonded14Atoms.end()) {
      this->bonded14Atoms.push_back(a);
    }
}

// ============================================================
// Function : getBondedAtoms()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> atom::getBondedAtoms()
{
    return this->bondedAtoms;
}

// ============================================================
// Function : get13BondedAtoms()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> atom::get13BondedAtoms()
{
    return this->bonded13Atoms;
}

// ============================================================
// Function : get14BondedAtoms()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> atom::get14BondedAtoms()
{
    return this->bonded14Atoms;
}

// ============================================================
// Function : getBondedHeavyAtoms()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> atom::getBondedHeavyAtoms()
{
    std::vector<atom*> heavyAtoms;
    for (unsigned int i = 0; i < this->bondedAtoms.size(); i++) {
      if (this->bondedAtoms[i]->getElementSymbol() != "H") {
        heavyAtoms.push_back(this->bondedAtoms[i]);
      }
    }
    return heavyAtoms;
}

// ============================================================
// Function : getNumBonds()
// ------------------------------------------------------------
//
// ============================================================
int atom::getNumBonds()
{
    return this->bondedAtoms.size();
}

// ============================================================
// Function : getNum13Bonds()
// ------------------------------------------------------------
//
// ============================================================
int atom::getNum13Bonds()
{
    return this->bonded13Atoms.size();
}

// ============================================================
// Function : getNum14Bonds()
// ------------------------------------------------------------
//
// ============================================================
int atom::getNum14Bonds()
{
    return this->bonded14Atoms.size();
}

// ============================================================
// Function : hasBondedAtom()
// ------------------------------------------------------------
//
// ============================================================
bool atom::hasBondedAtom(atom* a)
{
    std::vector<atom*>::iterator result;
    result = std::find(this->bondedAtoms.begin(), this->bondedAtoms.end(), a);
    if (result !=  this->bondedAtoms.end()) {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : has13BondedAtom()
// ------------------------------------------------------------
//
// ============================================================
bool atom::has13BondedAtom(atom* a)
{
    std::vector<atom*>::iterator result;
    result = std::find(this->bonded13Atoms.begin(), this->bonded13Atoms.end(), a);
    if (result !=  this->bonded13Atoms.end()) {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : has14BondedAtom()
// ------------------------------------------------------------
//
// ============================================================
bool atom::has14BondedAtom(atom* a)
{
    std::vector<atom*>::iterator result;
    result = std::find(this->bonded14Atoms.begin(), this->bonded14Atoms.end(), a);
    if (result !=  this->bonded14Atoms.end()) {
      return 1;
    }
    return 0;
}

} // MTKpp namespace
