/*!
   \file submolecule.cpp
   \brief Container for atoms and bonds
   \author Martin Peters

   $Date: 2010/04/29 18:59:17 $
   $Revision: 1.18 $

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

#include "submolecule.h"
#include "collection.h"
#include "molecule.h"
#include "atom.h"
#include "Utils/vector3d.h"

// - bond, angle, torsion, improper
#include "bond.h"
#include "angle.h"
#include "torsion.h"
#include "improper.h"

#include "stdFrag.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : submolecule()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
submolecule::submolecule(molecule *parent):pParent(parent)
{
    itsName            = "";
    its1LName          = "";
    itsiCode           = " ";
    itsColIndex        = 0;
    itsIndex           = 0;
    itsSubMolId        = 0;
    itsNumAtoms        = 0;
    itsNumBonds        = 0;
    pAtom              = 0;
    pBond              = 0;
    pStdFrag           = 0;
    pCenterMass        = 0;
}

// ============================================================
// Function : ~submolecule()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
submolecule::~submolecule()
{
    for (AtomIterator c = this->itsAtomList.begin();
                      c != this->itsAtomList.end(); c++) {
      pAtom = *c;
      delete pAtom;
    }
    itsAtomList.clear();
}

// ============================================================
// Function : copy()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::copy(submolecule* pSubMol)
{
    this->setName(pSubMol->getName());
    this->setIndex(pSubMol->getIndex());
    this->setSubMolId(pSubMol->getSubMolId());
    this->setiCode(pSubMol->getiCode());
    this->setStdFrag(pSubMol->getStdFrag());

    atom* pAtom1 = 0;
    atom* pAtom2 = 0;
    std::vector<atom*> atomList = pSubMol->getAtomList();
    for (AtomIterator c = atomList.begin(); c != atomList.end(); c++) {
      pAtom1 = *c;
      pAtom2 = this->addAtom();
      pAtom2->setElement(pAtom1->getElement());
      pAtom2->setName(pAtom1->getName());
      pAtom2->setCoords(pAtom1->getX(), pAtom1->getY(), pAtom1->getZ());
      pAtom2->setFileID(pAtom1->getFileID());
      pAtom2->setHybridization(pAtom1->getHybridization());
      pAtom2->setType(pAtom1->getType());

      if (pAtom1->getStdAtom()) {
        pAtom2->setStdAtom(pAtom1->getStdAtom());
      }
    }
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
//
// ============================================================
molecule* submolecule::getParent()
{
    return pParent;
}

// ============================================================
// Function : addAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* submolecule::addAtom()
{
    pAtom = new atom(this);

    // molecule index
    int iMol = pParent->getAtomIndex();
    if (static_cast<unsigned int>(iMol) > MAXATOMS) {
      std::string errMessage = " Too many atoms, increase MAXATOMS in constant.h and recompile ";
      errorLogger.throwError(" submolecule::addAtom ", errMessage, 1);

      std::stringstream ss;
      ss << " submolecule::addAtom "<< errMessage;
      throw MTKException(ss.str());
    }
    pAtom->setIndex(iMol);
    pParent->setAtomIndex(iMol+1);

    // collection index
    int iCol = pParent->getParent()->getAtomIndex();
    pAtom->setColIndex(iCol);
    pAtom->setFileID(iCol);
    pParent->getParent()->setAtomIndex(iCol+1);

    this->itsAtomList.push_back(pAtom);

    return pAtom;
}

// ============================================================
// Function : addAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* submolecule::addAtom(atom* pAt)
{
    pAtom = new atom(pAt, this);

    // molecule index
    int iMol = pParent->getAtomIndex();
    pAtom->setIndex(iMol);
    pParent->setAtomIndex(iMol+1);

    // collection index
    int iCol = pParent->getParent()->getAtomIndex();
    pAtom->setColIndex(iCol);
    pParent->getParent()->setAtomIndex(iCol+1);

    this->itsAtomList.push_back(pAtom);

    return pAtom;
}

// ============================================================
// Function : addBond()
// ------------------------------------------------------------
//
// ============================================================
Bond* submolecule::addBond(atom* at1, atom* at2)
{
    pBond = new Bond();
    pBond->atom1 = at1;
    pBond->atom2 = at2;
    pBond->type = 0;
    pBond->stereo = 0;
    pBond->topology = 0;
    itsBondList.push_back(pBond);

    return pBond;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setName(const std::string &name)
{
    this->itsName = name;
}

// ============================================================
// Function : set1LName()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::set1LName(const std::string &name)
{
    this->its1LName = name;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
//
// ============================================================
std::string submolecule::getName()
{
    return this->itsName;
}

// ============================================================
// Function : get1LName()
// ------------------------------------------------------------
//
// ============================================================
std::string submolecule::get1LName()
{
    return this->its1LName;
}

// ============================================================
// Function : setIndex()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setIndex(const int &id)
{
    this->itsIndex = id;
}

// ============================================================
// Function : getIndex()
// ------------------------------------------------------------
//
// ============================================================
int submolecule::getIndex()
{
    return this->itsIndex;
}

// ============================================================
// Function : setSubMolId()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setSubMolId(const int &id)
{
    this->itsSubMolId = id;
}

// ============================================================
// Function : getSubMolId()
// ------------------------------------------------------------
//
// ============================================================
int submolecule::getSubMolId()
{
    return this->itsSubMolId;
}

// ============================================================
// Function : setColIndex()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setColIndex(const int &id)
{
    this->itsColIndex = id;
}


// ============================================================
// Function : getColIndex()
// ------------------------------------------------------------
//
// ============================================================
int submolecule::getColIndex()
{
    return this->itsColIndex;
}

// ============================================================
// Function : setiCode()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setiCode(const std::string &i)
{
    this->itsiCode = i;
}

// ============================================================
// Function : getiCode()
// ------------------------------------------------------------
//
// ============================================================
std::string submolecule::getiCode()
{
    return this->itsiCode;
}

// ============================================================
// Function : setNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setNumAtoms(const int &natoms)
{
    this->itsNumAtoms = natoms;
}

// ============================================================
// Function : getNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
int submolecule::getNumAtoms()
{
    return this->itsAtomList.size();
}

// ============================================================
// Function : setNumBonds()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setNumBonds(const int &nbonds)
{
    this->itsNumBonds = nbonds;
}

// ============================================================
// Function : getNumBonds()
// ------------------------------------------------------------
//
// ============================================================
int submolecule::getNumBonds()
{
    return this->itsNumBonds;
}

// ============================================================
// Function : setStdFrag()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::setStdFrag(stdFrag* sf)
{
    this->pStdFrag = sf;
}

// ============================================================
// Function : getStdFrag()
// ------------------------------------------------------------
//
// ============================================================
stdFrag* submolecule::getStdFrag()
{
    return this->pStdFrag;
}

// ============================================================
// Function : hasStdFrag()
// ------------------------------------------------------------
//
// ============================================================
bool submolecule::hasStdFrag()
{
    if (this->pStdFrag) return true;
    return false;
}

// ============================================================
// Function : getAtomList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> submolecule::getAtomList()
{
    return this->itsAtomList;
}

// ============================================================
// Function : getBondList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<Bond*> submolecule::getBondList()
{
    return this->itsBondList;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* submolecule::getAtom(int number)
{
    for (AtomIterator c = itsAtomList.begin(); c != itsAtomList.end(); c++) {
      pAtom = *c;
      if (pAtom->getIndex() == number) {
        return pAtom;
      }
    }
    return 0;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* submolecule::getAtom(const std::string &name)
{
    for (AtomIterator c = itsAtomList.begin(); c != itsAtomList.end(); c++) {
      pAtom = *c;
      if (pAtom->getName() == name) {
        return pAtom;
      }
    }
    return 0;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* submolecule::getAtom(stdAtom* pStdAtom)
{
    stdAtom* pStdAt = 0;
    for (AtomIterator d = itsAtomList.begin(); d != itsAtomList.end(); d++) {
      pAtom = *d;
      pStdAt = pAtom->getStdAtom();
      if (pStdAt == pStdAtom) {
        return pAtom;
      }
    }
    return 0;
}

// ============================================================
// Function : centerOfMass()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::centerOfMass()
{
    vector3d* Coord;
    vector3d  Center;
    std::vector<atom*> atoms = this->getAtomList();
    int nAtoms = atoms.size();

    for ( int j = 0; j < nAtoms; j++) {
      pAtom = atoms[j];
      Coord = pAtom->getCoords();
      Center = Center+(*Coord);
    }

    Center = Center/nAtoms;
    pCenterMass  = new vector3d;
    *pCenterMass =  Center;
}

// ============================================================
// Function : getCenterOfMass()
// ------------------------------------------------------------
//
// ============================================================
vector3d* submolecule::getCenterOfMass()
{
    if (!pCenterMass) {
      this->centerOfMass();
    }
    return pCenterMass;
}

// ============================================================
// Function : numHeavyAtoms()
// ------------------------------------------------------------
//
// ============================================================
int submolecule::numHeavyAtoms()
{
    std::vector<atom*> atoms = this->getAtomList();
    int nAtoms = atoms.size();
    int nHeavy = 0;
    for (int j = 0; j < nAtoms; j++) {
      pAtom = atoms[j];
      if (pAtom->getElementSymbol() != "H") {
        nHeavy++;
      }
    }
    return nHeavy;
}

// ============================================================
// Function : getMolecularWeight()
// ------------------------------------------------------------
//
// ============================================================
double submolecule::getMolecularWeight()
{
    double w = 0;
    std::vector<atom*> atoms = this->getAtomList();
    int nAtoms = atoms.size();
    for (int j = 0; j < nAtoms; j++) {
      pAtom = atoms[j];
      w += pAtom->getAtomicMass();
    }
    return w;
}

// ============================================================
// Function : print()
// ------------------------------------------------------------
//
// ============================================================
void submolecule::print()
{
    std::cout << " Submolecule " << this->getName() << std::endl;
    std::vector<atom*> atoms = this->getAtomList();
    int nAtoms = atoms.size();
    for (int j = 0; j < nAtoms; j++) {
      pAtom = atoms[j];
      std::cout << pAtom->getName() << std::endl;
    }
}

} // MTKpp namespace

