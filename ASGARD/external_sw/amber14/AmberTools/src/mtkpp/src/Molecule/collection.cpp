/*!
   \file collection.cpp
   \brief Container for molecules
   \author Martin Peters

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.25 $

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

#include <sstream>

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "bond.h"
#include "angle.h"
#include "torsion.h"
#include "improper.h"

#include "Utils/vector3d.h"
#include "element.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "atomTyper.h"
#include "parameters.h"

#include "metalCenter.h"
#include "functionalize.h"
#include "superimpose.h"

#include "utility.h"
#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace MTKpp
{

// ============================================================
// Function : collection()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
collection::collection()
{
    this->itsName        = "";
    this->itsAtomIndex   = 1;
    this->itsSubMoleculeIndex = 1;

    // ELEMENTS
    this->pElements = new elements();

    // ============
    //      MM
    // ============
    this->pStdLibrary = 0;
    this->pParameters = 0;
    this->pAtomTyper  = 0;

    // =====================
    //      Metalloproteins
    // =====================
    availableMetals.push_back("Na");
    availableMetals.push_back("Mg");
    availableMetals.push_back("K");
    availableMetals.push_back("Ca");
    availableMetals.push_back("Mn");
    availableMetals.push_back("Fe");
    availableMetals.push_back("Co");
    availableMetals.push_back("Ni");
    availableMetals.push_back("Cu");
    availableMetals.push_back("Zn");

    // Ref. Harding, M. M., Acta Cryst. (2001). D57, 401-411
    // Ref: Harding, M. M., Acta Cryst. (2006). D62, 678-682
    metalDonorDists["Na:HOH"] = 2.41;
    metalDonorDists["Mg:HOH"] = 2.07;
    metalDonorDists["K:HOH"]  = 2.81;
    metalDonorDists["Ca:HOH"] = 2.39;
    metalDonorDists["Mn:HOH"] = 2.19;
    metalDonorDists["Fe:HOH"] = 2.09;
    metalDonorDists["Co:HOH"] = 2.09;
    metalDonorDists["Ni:HOH"] = 2.09;
    metalDonorDists["Cu:HOH"] = 2.13;
    metalDonorDists["Zn:HOH"] = 2.09;

    metalDonorDists["Na:w"] = 2.41;
    metalDonorDists["Mg:w"] = 2.07;
    metalDonorDists["K:w"]  = 2.81;
    metalDonorDists["Ca:w"] = 2.39;
    metalDonorDists["Mn:w"] = 2.19;
    metalDonorDists["Fe:w"] = 2.09;
    metalDonorDists["Co:w"] = 2.09;
    metalDonorDists["Ni:w"] = 2.09;
    metalDonorDists["Cu:w"] = 2.13;
    metalDonorDists["Zn:w"] = 2.09;

    metalDonorDists["Na:ASP"] = 2.41;
    metalDonorDists["Mg:ASP"] = 2.07;
    metalDonorDists["K:ASP"]  = 2.82;
    metalDonorDists["Ca:ASP"] = 2.36;
    metalDonorDists["Mn:ASP"] = 2.15;
    metalDonorDists["Fe:ASP"] = 2.04;
    metalDonorDists["Co:ASP"] = 2.05;
    metalDonorDists["Ni:ASP"] = 2.05;
    metalDonorDists["Cu:ASP"] = 1.99;
    metalDonorDists["Zn:ASP"] = 1.99;

    metalDonorDists["Na:D"] = 2.41;
    metalDonorDists["Mg:D"] = 2.07;
    metalDonorDists["K:D"]  = 2.82;
    metalDonorDists["Ca:D"] = 2.36;
    metalDonorDists["Mn:D"] = 2.15;
    metalDonorDists["Fe:D"] = 2.04;
    metalDonorDists["Co:D"] = 2.05;
    metalDonorDists["Ni:D"] = 2.05;
    metalDonorDists["Cu:D"] = 1.99;
    metalDonorDists["Zn:D"] = 1.99;

    metalDonorDists["Na:GLU"] = 2.41;
    metalDonorDists["Mg:GLU"] = 2.07;
    metalDonorDists["K:GLU"]  = 2.82;
    metalDonorDists["Ca:GLU"] = 2.36;
    metalDonorDists["Mn:GLU"] = 2.15;
    metalDonorDists["Fe:GLU"] = 2.04;
    metalDonorDists["Co:GLU"] = 2.05;
    metalDonorDists["Ni:GLU"] = 2.05;
    metalDonorDists["Cu:GLU"] = 1.99;
    metalDonorDists["Zn:GLU"] = 1.99;

    metalDonorDists["Na:E"] = 2.41;
    metalDonorDists["Mg:E"] = 2.07;
    metalDonorDists["K:E"]  = 2.82;
    metalDonorDists["Ca:E"] = 2.36;
    metalDonorDists["Mn:E"] = 2.15;
    metalDonorDists["Fe:E"] = 2.04;
    metalDonorDists["Co:E"] = 2.05;
    metalDonorDists["Ni:E"] = 2.05;
    metalDonorDists["Cu:E"] = 1.99;
    metalDonorDists["Zn:E"] = 1.99;

    metalDonorDists["Na:CRL"] = 2.38;
    metalDonorDists["Mg:CRL"] = 2.26;
    metalDonorDists["K:CRL"]  = 2.74;
    metalDonorDists["Ca:CRL"] = 2.36;
    metalDonorDists["Mn:CRL"] = 2.19;
    metalDonorDists["Fe:CRL"] = 2.04;
    metalDonorDists["Co:CRL"] = 2.08;
    metalDonorDists["Ni:CRL"] = 2.08;
    metalDonorDists["Cu:CRL"] = 2.04;
    metalDonorDists["Zn:CRL"] = 2.07;

    metalDonorDists["Mn:HIS"] = 2.21;
    metalDonorDists["Fe:HIS"] = 2.16;
    metalDonorDists["Co:HIS"] = 2.14;
    metalDonorDists["Ni:HIS"] = 2.14;
    metalDonorDists["Cu:HIS"] = 2.02;
    metalDonorDists["Zn:HIS"] = 2.03;

    metalDonorDists["Mn:H"] = 2.21;
    metalDonorDists["Fe:H"] = 2.16;
    metalDonorDists["Co:H"] = 2.14;
    metalDonorDists["Ni:H"] = 2.14;
    metalDonorDists["Cu:H"] = 2.02;
    metalDonorDists["Zn:H"] = 2.03;

    metalDonorDists["Mn:CYS"] = 2.35;
    metalDonorDists["Fe:CYS"] = 2.30;
    metalDonorDists["Co:CYS"] = 2.25;
    metalDonorDists["Ni:CYS"] = 2.25;
    metalDonorDists["Cu:CYS"] = 2.15;
    metalDonorDists["Zn:CYS"] = 2.31;

    metalDonorDists["Mn:C"] = 2.35;
    metalDonorDists["Fe:C"] = 2.30;
    metalDonorDists["Co:C"] = 2.25;
    metalDonorDists["Ni:C"] = 2.25;
    metalDonorDists["Cu:C"] = 2.15;
    metalDonorDists["Zn:C"] = 2.31;

    metalDonorDists["Mn:MET"] = 2.35;
    metalDonorDists["Fe:MET"] = 2.30;
    metalDonorDists["Co:MET"] = 2.25;
    metalDonorDists["Ni:MET"] = 2.25;
    metalDonorDists["Cu:MET"] = 2.15;
    metalDonorDists["Zn:MET"] = 2.31;

    metalDonorDists["Mn:M"] = 2.35;
    metalDonorDists["Fe:M"] = 2.30;
    metalDonorDists["Co:M"] = 2.25;
    metalDonorDists["Ni:M"] = 2.25;
    metalDonorDists["Cu:M"] = 2.15;
    metalDonorDists["Zn:M"] = 2.31;

    metalDonorDists["Ca:SER"] = 2.43;
    metalDonorDists["Mg:SER"] = 2.10;
    metalDonorDists["Mn:SER"] = 2.25;
    metalDonorDists["Fe:SER"] = 2.13;
    metalDonorDists["Co:SER"] = 2.10;
    metalDonorDists["Ni:SER"] = 2.10;
    metalDonorDists["Cu:SER"] = 2.00;
    metalDonorDists["Zn:SER"] = 1.14; // should be 2.14

    metalDonorDists["Ca:S"] = 2.43;
    metalDonorDists["Mg:S"] = 2.10;
    metalDonorDists["Mn:S"] = 2.25;
    metalDonorDists["Fe:S"] = 2.13;
    metalDonorDists["Co:S"] = 2.10;
    metalDonorDists["Ni:S"] = 2.10;
    metalDonorDists["Cu:S"] = 2.00;
    metalDonorDists["Zn:S"] = 1.14; // should be 2.14

    metalDonorDists["Ca:THR"] = 2.43;
    metalDonorDists["Mg:THR"] = 2.10;
    metalDonorDists["Mn:THR"] = 2.25;
    metalDonorDists["Fe:THR"] = 2.13;
    metalDonorDists["Co:THR"] = 2.10;
    metalDonorDists["Ni:THR"] = 2.10;
    metalDonorDists["Cu:THR"] = 2.00;
    metalDonorDists["Zn:THR"] = 2.14;

    metalDonorDists["Ca:T"] = 2.43;
    metalDonorDists["Mg:T"] = 2.10;
    metalDonorDists["Mn:T"] = 2.25;
    metalDonorDists["Fe:T"] = 2.13;
    metalDonorDists["Co:T"] = 2.10;
    metalDonorDists["Ni:T"] = 2.10;
    metalDonorDists["Cu:T"] = 2.00;
    metalDonorDists["Zn:T"] = 2.14;

    metalDonorDists["Ca:TYR"] = 2.20;
    metalDonorDists["Mg:TYR"] = 1.87;
    metalDonorDists["Mn:TYR"] = 1.88;
    metalDonorDists["Fe:TYR"] = 1.93;
    metalDonorDists["Co:TYR"] = 1.90;
    metalDonorDists["Ni:TYR"] = 1.90;
    metalDonorDists["Cu:TYR"] = 1.90;
    metalDonorDists["Zn:TYR"] = 1.95;

    metalDonorDists["Ca:Y"] = 2.20;
    metalDonorDists["Mg:Y"] = 1.87;
    metalDonorDists["Mn:Y"] = 1.88;
    metalDonorDists["Fe:Y"] = 1.93;
    metalDonorDists["Co:Y"] = 1.90;
    metalDonorDists["Ni:Y"] = 1.90;
    metalDonorDists["Cu:Y"] = 1.90;
    metalDonorDists["Zn:Y"] = 1.95;
}

// ============================================================
// Function : ~collection()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
collection::~collection()
{
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      delete pMolecule;
    }
    this->itsMoleculeList.clear();

    if (pElements) {
      delete pElements;
    }
/*
    if (pStdLibrary) {
      delete pStdLibrary;
    }

    if (pParameters) {
      delete pParameters;
    }
*/
    if (pAtomTyper) {
      delete pAtomTyper;
    }
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
//
// ============================================================
void collection::setName(std::string name)
{
    this->itsName = name;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
//
// ============================================================
std::string collection::getName()
{
    return this->itsName;
}

// ============================================================
// Function : addMolecule()
// ------------------------------------------------------------
// Add Molecule to Collection
// ============================================================
molecule* collection::addMolecule()
{
    pMolecule = new molecule(this);

    this->itsMoleculeList.push_back(pMolecule);

    return pMolecule;
}

// ============================================================
// Function : delMolecule()
// ------------------------------------------------------------
// Delete Molecule from Collection
// ============================================================
void collection::delMolecule(molecule* pMol)
{
    moleculeIterator c = std::find(this->itsMoleculeList.begin(),
          this->itsMoleculeList.end(), pMol);

    if (c != this->itsMoleculeList.end()) {
      delete pMol;
      this->itsMoleculeList.erase(c);
      pMol = 0;
    }
}

// ============================================================
// Function : delAllMolecules()
// ------------------------------------------------------------
// Delete all Molecules from Collection
// ============================================================
void collection::delAllMolecules()
{
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      delete pMolecule;
    }
    itsMoleculeList.clear();
    pMolecule = 0;
    pSubMolecule = 0;
    itsAtomIndex   = 1;
}

// ============================================================
// Function : clear()
// ------------------------------------------------------------
// Clears collection
// ============================================================
void collection::clear()
{
    this->delAllMolecules();

    this->itsMetalAtoms.erase(this->itsMetalAtoms.begin(), this->itsMetalAtoms.end());
    this->itsMetalAtoms.clear();

    for (MetalGroupIterator c = this->itsMetalGroups.begin();
          c != this->itsMetalGroups.end(); c++) {
      if (*c) {
        metalGroup* pMG = *c;
        delete pMG;
      }
    }
    itsMetalGroups.clear();

    for (MetalCenterIterator c = this->itsMetalCenters.begin();
          c != this->itsMetalCenters.end(); c++) {
      if (*c) {
        metalCenter* pMC = *c;
        delete pMC;
      }
    }
    itsMetalCenters.clear();
}

// ============================================================
// Function : getMolecule()
// ------------------------------------------------------------
//
// ============================================================
molecule* collection::getMolecule()
{
    return this->itsMoleculeList.back();
}

// ============================================================
// Function : getMolecule()
// ------------------------------------------------------------
//
// ============================================================
molecule* collection::getMolecule(int id)
{
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      if (pMolecule->getMolId() == id) {
        return pMolecule;
      }
    }
    return 0;
}

// ============================================================
// Function : getMolecule()
// ------------------------------------------------------------
//
// ============================================================
molecule* collection::getMolecule(std::string name)
{
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      if (pMolecule->getName() == name) {
        return pMolecule;
      }
    }
    return 0;
}

// ============================================================
// Function : getLastAddedMolecule()
// ------------------------------------------------------------
//
// ============================================================
molecule* collection::getLastAddedMolecule()
{
    if (!this->itsMoleculeList.size() < 1) {
      return this->itsMoleculeList.back();
    }
    return 0;
}

// ============================================================
// Function : getMoleculeList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<molecule*> collection::getMoleculeList()
{
    return this->itsMoleculeList;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* collection::getAtom(int number, bool atomIndex, bool fileId, bool atomColIndex)
{
    for (moleculeIterator m = this->itsMoleculeList.begin();
         m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;

      if (pMolecule->getName() == "Reference") {
        continue;  
      }

      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (sMolIterator c = subMolList.begin(); c != subMolList.end(); c++) {
        pSubMolecule = *c;
        std::vector<atom*> atList = pSubMolecule->getAtomList();
        for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
          atom * pAtom = *d;
          if (atomColIndex && (pAtom->getColIndex() == number)) {
            return pAtom;
          }
          if (atomIndex && (pAtom->getIndex() == number)) {
            return pAtom;
          }
          if (fileId && (pAtom->getFileID() == number)) {
            return pAtom;
          }
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getAtomList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> collection::getAtomList()
{
    std::vector<atom*> atList;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<atom*> molAtList = pMolecule->getAtomList();
      for (unsigned int i = 0; i < molAtList.size(); i++) {
        atList.push_back(molAtList[i]);
      }
    }
    return atList;
}

// ============================================================
// Function : getNumberMolecules()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumberMolecules()
{
    int nMols = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      //std::cout << "|" << pMolecule->getName() << "| " << std::endl;
      if (pMolecule->getName() != "Reference") {
        nMols++;
      }
    }

    return nMols;
    //return itsMoleculeList.size();
}

// ============================================================
// Function : getNumberMolecules()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumberSubMolecules()
{
    int nsmol = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      nsmol += pMolecule->getNumSubMolecules();
    }
    return nsmol;
}

// ============================================================
// Function : addStdLibrary()
// ------------------------------------------------------------
//
// ============================================================
void collection::addStdLibrary()
{
    if (pStdLibrary == 0) {
      //pStdLibrary = new stdLibrary();
      pStdLibrary = stdLibrary::getInstance();
    }
}

// ============================================================
// Function : getStdLibrary()
// ------------------------------------------------------------
//
// ============================================================
stdLibrary* collection::getStdLibrary()
{
    return pStdLibrary;
}

// ============================================================
// Function : addAtomTyper()
// ------------------------------------------------------------
//
// ============================================================
void collection::addAtomTyper()
{
    if (pAtomTyper == 0) {
      pAtomTyper = new atomTyper();
    }
}

// ============================================================
// Function : getAtomTyper()
// ------------------------------------------------------------
//
// ============================================================
atomTyper* collection::getAtomTyper()
{
    return pAtomTyper;
}

// ============================================================
// Function : addParameters()
// ------------------------------------------------------------
//
// ============================================================
void collection::addParameters()
{
    if (this->pParameters == 0) {
      //this->pParameters = new parameters(pElements);
      this->pParameters = parameters::getInstance(pElements);
    }
}

// ============================================================
// Function : getParameters()
// ------------------------------------------------------------
//
// ============================================================
parameters* collection::getParameters()
{
    if (!this->pParameters) {
      errorLogger.throwError("collection::getParameters", "Can't find parameters", 2);
      return 0;
    }
    return this->pParameters;
}

// ============================================================
// Function : getNumMMnonBondedPairs()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumMMnonBondedPairs(double cutoff)
{
    errorLogger.throwError("collection::getNumMMnonBondedPairs", " ", 4);

    int nBPs = 0;
    unsigned int nMolecules = this->itsMoleculeList.size();

    double cutoff2 = cutoff*cutoff;
    double subMolCutoff = (cutoff*2)*(cutoff*2);
    std::vector<submolecule*> subMoleculeList;

    // Compute all submolecule centroids
    for (moleculeIterator c = itsMoleculeList.begin(); c != itsMoleculeList.end(); c++) {
      pMolecule = *c;
      subMoleculeList = pMolecule->getSubMoleculeList();
      for (unsigned int d = 0; d < subMoleculeList.size(); d++) {
        subMoleculeList[d]->centerOfMass();
      }
    }

    for (unsigned int i = 0; i < nMolecules; i++) {
      std::vector<submolecule*> subMolList1 = itsMoleculeList[i]->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList1.size(); j++) {
        std::vector<atom*> atList1 = subMolList1[j]->getAtomList();
        for (unsigned int k = 0; k < atList1.size(); k++) {
          for (unsigned int l = i; l < nMolecules; l++) {
            std::vector<submolecule*> subMolList2 = itsMoleculeList[l]->getSubMoleculeList();
            // SAME MOLECULE
            if (itsMoleculeList[i] == itsMoleculeList[l]) {
              for (unsigned int q = j; q < subMolList1.size(); q++) {
                std::vector<atom*> atList2 = subMolList1[q]->getAtomList();
                // INTRA RESIDUE
                if ((subMolList1[j] == subMolList1[q]) and (atList1.size() > 1)) {
                  for (unsigned int w = k+1; w < atList1.size(); w++) {
                    if (atList1[k]->getCoords()->dist2( (*atList1[w]->getCoords())) < cutoff2) {
                      if ( (!atList1[k]->hasBondedAtom(atList1[w])) && (!atList1[k]->has13BondedAtom(atList1[w]))) {
                        nBPs++;
                      }
                    }
                  }
                }
                else {
                  if (subMolList1[j]->getCenterOfMass()->dist2( (*subMolList1[q]->getCenterOfMass())) < subMolCutoff) {
                    for (unsigned int w = 0; w < atList2.size(); w++) {
                      if (atList1[k]->getCoords()->dist2( (*atList2[w]->getCoords())) < cutoff2) {
                        if ( (!atList1[k]->hasBondedAtom(atList2[w])) && (!atList1[k]->has13BondedAtom(atList2[w]))) {
                          nBPs++;
                        }
                      }
                    }
                  }
                }
              }
            }
            // DIFFERENT MOLECULES
            else {
              std::vector<submolecule*> subMolList2 = itsMoleculeList[l]->getSubMoleculeList();
              for (unsigned int q = 0; q < subMolList2.size(); q++) {
                if (subMolList1[j]->getCenterOfMass()->dist2( (*subMolList2[q]->getCenterOfMass())) < subMolCutoff) {
                  std::vector<atom*> atList2 = subMolList2[q]->getAtomList();
                  for (unsigned int w = 0; w < atList2.size(); w++) {
                    if (atList1[k]->getCoords()->dist2( (*atList2[w]->getCoords())) < cutoff2) {
                      if ( (!atList1[k]->hasBondedAtom(atList2[w])) && (!atList1[k]->has13BondedAtom(atList2[w]))) {
                        nBPs++;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return nBPs;
}

// ============================================================
// Function : getMMnonBondedPairs()
// ------------------------------------------------------------
//
// ============================================================
int collection::getMMnonBondedPairs(int nonBonded[], int nonBondedPtrs[],
                                    double nonBondedParams[],
                                    int nonBonded14Ptrs[], double cutoff)
{
    errorLogger.throwError("collection::getMMnonBondedPairs", " ", 4);

    unsigned int nMolecules = this->itsMoleculeList.size();

    double cutoff2 = cutoff*cutoff;
    double subMolCutoff = (cutoff*2)*(cutoff*2);
    std::vector<submolecule*> subMoleculeList;

    int index = 0;
    int indexParam = 0;
    int nNBPtr = 0;
    int n = 0;

    try {
      for (unsigned int i = 0; i < nMolecules; i++) {
        std::vector<submolecule*> subMolList1 = this->itsMoleculeList[i]->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList1.size(); j++) {
          std::vector<atom*> atList1 = subMolList1[j]->getAtomList();
          for (unsigned int k = 0; k < atList1.size(); k++) {
            nNBPtr = 0;
            for (unsigned int l = i; l < nMolecules; l++) {
              std::vector<submolecule*> subMolList2 = itsMoleculeList[l]->getSubMoleculeList();
              // SAME MOLECULE
              if (itsMoleculeList[i] == itsMoleculeList[l]) {
                for (unsigned int q = j; q < subMolList1.size(); q++) {
                  std::vector<atom*> atList2 = subMolList1[q]->getAtomList();
                  // INTRA RESIDUE
                  if ((subMolList1[j] == subMolList1[q]) and (atList1.size() > 1)) {
                    for (unsigned int w = k+1; w < atList1.size(); w++) {
                      if (atList1[k]->getCoords()->dist2(
                           (*atList1[w]->getCoords())) < cutoff2) {
                        if ( (!atList1[k]->hasBondedAtom(atList1[w])) &&
                             (!atList1[k]->has13BondedAtom(atList1[w]))) {
                          pLJ612SE = pParameters->getLJ612SE(atList1[k]->getStdAtom()->type,
                                     atList1[w]->getStdAtom()->type);
                          if (pLJ612SE) {
                            nonBonded[index] = atList1[w]->getColIndex()-1;
                            nonBondedParams[indexParam] = pLJ612SE->sigma;
                            nonBondedParams[indexParam+1] = pLJ612SE->epsilon;
                            if (atList1[k]->has14BondedAtom(atList1[w])) {
                              nonBonded14Ptrs[index] = 1;
                            }
                            else {
                              nonBonded14Ptrs[index] = 0;
                            }
                            nNBPtr++;
                            index++;
                            indexParam+=2;
                          }
                          else {
                            nonBonded[index] = -1;
                            nonBondedParams[indexParam] = 0;
                            nonBondedParams[indexParam+1] = 0;
                            nonBonded14Ptrs[index] = 0;
                            nNBPtr++;
                            index++;
                            indexParam+=2;
                          }
                        }
                      }
                    }
                  }
                  // INTER RESIDUE
                  else {
                    if (subMolList1.size() > 1) {
                      if (subMolList1[j]->getCenterOfMass()->dist2(
                              (*subMolList1[q]->getCenterOfMass()))
                               < subMolCutoff) {
                        for (unsigned int w = 0; w < atList2.size(); w++) {
                          if (atList1[k]->getCoords()->dist2(
                               (*atList2[w]->getCoords())) < cutoff2) {
                            if ( (!atList1[k]->hasBondedAtom(atList2[w])) &&
                                 (!atList1[k]->has13BondedAtom(atList2[w]))) {
                              pLJ612SE = pParameters->getLJ612SE(
                                          atList1[k]->getStdAtom()->type,
                                          atList2[w]->getStdAtom()->type);

                              if (pLJ612SE) {
                                nonBonded[index] = atList2[w]->getColIndex()-1;
                                nonBondedParams[indexParam] = pLJ612SE->sigma;
                                nonBondedParams[indexParam+1] = pLJ612SE->epsilon;
                                if (atList1[k]->has14BondedAtom(atList2[w])) {
                                  nonBonded14Ptrs[index] = 1;
                                }
                                else {
                                  nonBonded14Ptrs[index] = 0;
                                }
                                nNBPtr++;
                                index++;
                                indexParam+=2;
                              }
                              else {
                                nonBonded[index] = -1;
                                nonBondedParams[indexParam] = 0;
                                nonBondedParams[indexParam+1] = 0;
                                nonBonded14Ptrs[index] = 0;
                                nNBPtr++;
                                index++;
                                indexParam+=2;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              // DIFFERENT MOLECULES
              else {
                std::vector<submolecule*> subMolList2 = itsMoleculeList[l]->getSubMoleculeList();
                for (unsigned int q = 0; q < subMolList2.size(); q++) {
                  if (subMolList1[j]->getCenterOfMass()->dist2(
                       (*subMolList2[q]->getCenterOfMass())) < subMolCutoff) {
                    std::vector<atom*> atList2 = subMolList2[q]->getAtomList();
                    for (unsigned int w = 0; w < atList2.size(); w++) {
                      if (atList1[k]->getCoords()->dist2( (*atList2[w]->getCoords())) < cutoff2) {
                        if ( (!atList1[k]->hasBondedAtom(atList2[w])) &&
                             (!atList1[k]->has13BondedAtom(atList2[w]))) {
                          if (atList1[k]->getStdAtom() and atList2[w]->getStdAtom()) {
                            pLJ612SE = pParameters->getLJ612SE(atList1[k]->getStdAtom()->type,
                                       atList2[w]->getStdAtom()->type);
                          }
                          else {
                            pLJ612SE = 0;
                          }
                          if (pLJ612SE) {
                            nonBonded[index] = atList2[w]->getColIndex()-1;

                            nonBondedParams[indexParam] = pLJ612SE->sigma;
                            nonBondedParams[indexParam+1] = pLJ612SE->epsilon;
                            if (atList1[k]->has14BondedAtom(atList2[w])) {
                              nonBonded14Ptrs[index] = 1;
                            }
                            else {
                              nonBonded14Ptrs[index] = 0;
                            }
                            nNBPtr++;
                            index++;
                            indexParam+=2;
                          }
                          else {
                            nonBonded[index] = -1;
                            nonBondedParams[indexParam] = 0;
                            nonBondedParams[indexParam+1] = 0;
                            nonBonded14Ptrs[index] = 0;
                            nNBPtr++;
                            index++;
                            indexParam+=2;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            nonBondedPtrs[n] = nNBPtr;
            n++;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumAtoms()
{
    int nAtoms = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      nAtoms += pMolecule->getNumAtoms();
    }
    return nAtoms;
}

// ============================================================
// Function : getNumBonds()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumBonds()
{
    int nBonds = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      nBonds += pMolecule->numBonds();
    }

    // Add metal cluster bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      nBonds += this->itsMetalCenters[m]->numBonds();
    }
    return nBonds;
}

// ============================================================
// Function : getBonds()
// ------------------------------------------------------------
//
// ============================================================
int collection::getBonds(int bonds[])
{
    int i = 0;
    std::map<int, Bond*> bondMap;

    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
           c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        bondMap = pMolecule->getBondMap();

        if (!bondMap.empty()) {
          for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
            pBond = b->second;
            bonds[i]   = pBond->atom1->getColIndex()-1;
            bonds[i+1] = pBond->atom2->getColIndex()-1;
            i+=2;
          }
        }
      }

      // Add metal cluster bonds
      if (this->itsMetalAtoms.size() != this->itsMetalCenters.size()) {
        std::cout << " Error in getBonds " << std::endl;
        //exit(0);
        throw MTKException(" Error in getBonds ");
      }
      for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
        atom* met = this->itsMetalAtoms[m];
        std::vector<atom*> primShellAtoms;
        this->itsMetalCenters[m]->getPrimaryShellAtoms(primShellAtoms);
        for (unsigned int mm = 0; mm < primShellAtoms.size(); mm++) {
          bonds[i]   = met->getColIndex()-1;
          bonds[i+1] = primShellAtoms[mm]->getColIndex()-1;
          i+=2;
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getBondParams()
// ------------------------------------------------------------
//
// ============================================================
int collection::getBondParams(double bondParams[])
{
    int i = 0;
    std::map<int, Bond*> bondMap;

    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
           c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        bondMap = pMolecule->getBondMap();

        if (!bondMap.empty()) {
          for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
            pBond = b->second;
            pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                   pBond->atom2->getStdAtom()->type);

            if (pBondParam) {
              bondParams[i]   = pBondParam->req;
              bondParams[i+1] = pBondParam->keq;
            }
            else {
              bondParams[i]   = 0.0;
              bondParams[i+1] = 0.0;
            }
            i+=2;
          }
        }
      }

      // Add metal cluster bonds
      for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
        atom* met = this->itsMetalAtoms[m];
        std::vector<atom*> primShellAtoms;
        this->itsMetalCenters[m]->getPrimaryShellAtoms(primShellAtoms);
        for (unsigned int mm = 0; mm < primShellAtoms.size(); mm++) {
          pBondParam = pParameters->getBondParam(met->getStdAtom()->type,
                                                 primShellAtoms[mm]->getStdAtom()->type);
          if (pBondParam) {
            bondParams[i]   = pBondParam->req;
            bondParams[i+1] = pBondParam->keq;
          }
          else {
            bondParams[i]   = 0.0;
            bondParams[i+1] = 0.0;
          }
          i+=2;
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumBondsWithH()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumBondsWithH()
{
    std::map<int, Bond*> bondMap;

    int nBondWithH = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      bondMap = pMolecule->getBondMap();

      if (!bondMap.empty()) {
        for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
          pBond = b->second;
          pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                 pBond->atom2->getStdAtom()->type);
          if (pBondParam) {
            if ((pBond->atom1->getElementSymbol() == "H") or
                (pBond->atom2->getElementSymbol() == "H")) {
              nBondWithH++;
            }
          }
          else {
            std::cout << " Can't find bond parameter " << std::endl;
          }
        }
      }
    }
    return nBondWithH;
}

// ============================================================
// Function : getNumBondsWithOutH()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumBondsWithOutH()
{
    std::map<int, Bond*> bondMap;

    int nBondWithOutH = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      bondMap = pMolecule->getBondMap();

      if (!bondMap.empty()) {
        for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
          pBond = b->second;
          pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                 pBond->atom2->getStdAtom()->type);
          if (pBondParam) {
            if (pBond->atom1->getElementSymbol() != "H" and
                pBond->atom2->getElementSymbol() != "H") {
              nBondWithOutH++;
            }
          }
          else {
            std::cout << " Can't find bond parameter " << std::endl;
          }
        }
      }
    }

    // Add metal cluster bonds, assuming no metal-H bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      bondMap = this->itsMetalCenters[m]->getBondMap();

      if (!bondMap.empty()) {
        for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
          pBond = b->second;
          pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                 pBond->atom2->getStdAtom()->type);
          if (pBondParam) {
            if (pBond->atom1->getElementSymbol() != "H" and
                pBond->atom2->getElementSymbol() != "H") {
              pBond->pBondParam = pBondParam;
              nBondWithOutH++;
            }
          }
          else {
            std::cout << " Can't find bond parameter " << std::endl;
          }
        }
      }
    }
    return nBondWithOutH;
}

// ============================================================
// Function : getNumUniqueBondTypes()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumUniqueBondTypes()
{
    std::map<int, Bond*> bondMap;
    typedef std::vector<bondParam*>::iterator bpVectorIterator;

    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      bondMap = pMolecule->getBondMap();

      if (!bondMap.empty()) {
        for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
          pBond = b->second;
          pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                 pBond->atom2->getStdAtom()->type);
          if (pBondParam) {
            bpVectorIterator r = std::find(uniqueBondParams.begin(), uniqueBondParams.end(), pBondParam);
            if (r == uniqueBondParams.end()) {
              uniqueBondParams.push_back(pBondParam);
            }
          }
          else {
            std::cout << " Can't find bond parameter " << std::endl;
          }
        }
      }
    }

    // Add metal cluster bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      bondMap = this->itsMetalCenters[m]->getBondMap();

      if (!bondMap.empty()) {
        for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
          pBond = b->second;
          pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                 pBond->atom2->getStdAtom()->type);
          if (pBondParam) {
            bpVectorIterator r = std::find(uniqueBondParams.begin(), uniqueBondParams.end(), pBondParam);
            if (r == uniqueBondParams.end()) {
              uniqueBondParams.push_back(pBondParam);
            }
          }
          else {
            std::cout << " Can't find bond parameter " << std::endl;
          }
        }
      }
    }
    return static_cast<int>(uniqueBondParams.size());
}


// ============================================================
// Function : getNumAngles()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumAngles()
{
    int nAngles = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      nAngles += pMolecule->numAngles();
    }

    // Add metal cluster angles
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      nAngles += this->itsMetalCenters[m]->numAngles();
    }

    return nAngles;
}

// ============================================================
// Function : getAngles()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAngles(int angles[])
{
    int i = 0;
    std::map<ULONG_KIND, Angle*> angleMap;

    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        angleMap = pMolecule->getAngleMap();

        if (!angleMap.empty()) {
          for (AngleMapIterator a = angleMap.begin();
               a != angleMap.end(); a++) {
            pAngle = a->second;
            angles[i]   = pAngle->atom1->getColIndex()-1;
            angles[i+1] = pAngle->atom2->getColIndex()-1;
            angles[i+2] = pAngle->atom3->getColIndex()-1;
/*
std::cout << " angle: " << pAngle->atom1->getColIndex()-1 << "-" 
          << pAngle->atom2->getColIndex()-1 << "-"
          << pAngle->atom3->getColIndex()-1 <<std::endl;
*/
            i+=3;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getAngleParams()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAngleParams(double angleParams[])
{
    int i = 0;
    std::map<ULONG_KIND, Angle*> angleMap;

    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        angleMap = pMolecule->getAngleMap();

        if (!angleMap.empty()) {
          for (AngleMapIterator a = angleMap.begin();
               a != angleMap.end(); a++) {
            pAngle = a->second;
            pAngleParam = pParameters->getAngleParam(
               pAngle->atom1->getStdAtom()->type,
               pAngle->atom2->getStdAtom()->type,
               pAngle->atom3->getStdAtom()->type);

            if (pAngleParam) {
              angleParams[i]   = pAngleParam->req;
              angleParams[i+1] = pAngleParam->keq;
            }
            else {
              angleParams[i]   = 0.0;
              angleParams[i+1] = 0.0;
            }
            i+=2;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumAnglesWithH()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumAnglesWithH()
{
    int nAnglesWithH = 0;

    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();

      if (!angleMap.empty()) {
        for (AngleMapIterator a = angleMap.begin();
             a != angleMap.end(); a++) {
          pAngle = a->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            if ((pAngle->atom1->getElementSymbol() == "H") or
                (pAngle->atom3->getElementSymbol() == "H")) {
              nAnglesWithH++;
            }
          }
          else {
            std::cout << pAngle->atom1->getStdAtom()->type << " "
                      << pAngle->atom2->getStdAtom()->type << " "
                      << pAngle->atom3->getStdAtom()->type << std::endl;
            std::cout << " collection::getNumAnglesWithH:1:Can't find angle parameter " << std::endl;
          }
        }
      }
    }

    // Add metal cluster bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      std::map<ULONG_KIND, Angle*> angleMap = this->itsMetalCenters[m]->getAngleMap();

      if (!angleMap.empty()) {
        for (AngleMapIterator a = angleMap.begin();
             a != angleMap.end(); a++) {
          pAngle = a->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            if ((pAngle->atom1->getElementSymbol() == "H") or
                (pAngle->atom3->getElementSymbol() == "H")) {
              pAngle->pAngleParam = pAngleParam;
              nAnglesWithH++;
            }
          }
          else {
            std::cout << pAngle->atom1->getStdAtom()->type << " "
                      << pAngle->atom2->getStdAtom()->type << " "
                      << pAngle->atom3->getStdAtom()->type << std::endl;
            std::cout << " collection::getNumAnglesWithH:2:Can't find angle parameter " << std::endl;
          }
        }
      }
    }
    return nAnglesWithH;
}

// ============================================================
// Function : getNumAnglesWithOutH()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumAnglesWithOutH()
{
    int nAnglesWithOutH = 0;

    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();

      if (!angleMap.empty()) {
        for (AngleMapIterator a = angleMap.begin();
             a != angleMap.end(); a++) {
          pAngle = a->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            if ((pAngle->atom1->getElementSymbol() != "H") and
                (pAngle->atom3->getElementSymbol() != "H")) {
              nAnglesWithOutH++;
            }
          }
          else {
            std::cout << pAngle->atom1->getStdAtom()->type << " "
                      << pAngle->atom2->getStdAtom()->type << " "
                      << pAngle->atom3->getStdAtom()->type << std::endl;
            std::cout << " collection::getNumAnglesWithOutH:1:Can't find angle parameter " << std::endl;
          }
        }
      }
    }

    // Add metal cluster bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      std::map<ULONG_KIND, Angle*> angleMap = this->itsMetalCenters[m]->getAngleMap();

      if (!angleMap.empty()) {
        for (AngleMapIterator a = angleMap.begin();
             a != angleMap.end(); a++) {
          pAngle = a->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            if ((pAngle->atom1->getElementSymbol() != "H") and
                (pAngle->atom3->getElementSymbol() != "H")) {
              std::cout <<" collection::getNumAnglesWithOutH: Adding " << pAngle->atom1->getName() << " "
                      << pAngle->atom2->getName() << " "
                      << pAngle->atom3->getName() << std::endl;
              pAngle->pAngleParam = pAngleParam;
              nAnglesWithOutH++;
            }
          }
          else {
            std::cout << pAngle->atom1->getStdAtom()->type << " "
                      << pAngle->atom2->getStdAtom()->type << " "
                      << pAngle->atom3->getStdAtom()->type << std::endl;
            std::cout << " collection::getNumAnglesWithOutH:2:Can't find angle parameter " << std::endl;
          }
        }
      }
    }

    return nAnglesWithOutH;
}

// ============================================================
// Function : getNumUniqueAngleTypes()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumUniqueAngleTypes()
{
    typedef std::vector<angleParam*>::iterator apVectorIterator;

    for (moleculeIterator c = this->itsMoleculeList.begin();
         c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();

      if (!angleMap.empty()) {
        for (AngleMapIterator a = angleMap.begin();
             a != angleMap.end(); a++) {
          pAngle = a->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            apVectorIterator r = std::find(uniqueAngleParams.begin(), uniqueAngleParams.end(), pAngleParam);
            if (r == uniqueAngleParams.end()) {
              uniqueAngleParams.push_back(pAngleParam);
            }
          }
          else {
            std::cout << pAngle->atom1->getStdAtom()->type << " "
                      << pAngle->atom2->getStdAtom()->type << " "
                      << pAngle->atom3->getStdAtom()->type << std::endl;
            std::cout << " collection::getNumUniqueAngleTypes:1:Can't find angle parameter " << std::endl;
          }
        }
      }
    }

    // Add metal cluster bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      std::map<ULONG_KIND, Angle*> angleMap = this->itsMetalCenters[m]->getAngleMap();

      if (!angleMap.empty()) {
        for (AngleMapIterator a = angleMap.begin(); a != angleMap.end(); a++) {
          pAngle = a->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            apVectorIterator r = std::find(uniqueAngleParams.begin(), uniqueAngleParams.end(), pAngleParam);
            if (r == uniqueAngleParams.end()) {
              uniqueAngleParams.push_back(pAngleParam);
            }
          }
          else {
            std::cout << pAngle->atom1->getStdAtom()->type << " "
                      << pAngle->atom2->getStdAtom()->type << " "
                      << pAngle->atom3->getStdAtom()->type << std::endl;
            std::cout << " collection::getNumUniqueAngleTypes:2:Can't find angle parameter " << std::endl;
          }
        }
      }
    }

    return static_cast<int>(uniqueAngleParams.size());
}

// ============================================================
// Function : getNumTorsions()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumTorsions()
{
    int nTorsions = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      nTorsions += pMolecule->numTorsions();
    }
    return nTorsions;
}

// ============================================================
// Function : getNumMMTorsions()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumMMTorsions()
{
    int nTorsions = 0;
    std::map<ULONG_KIND, Torsion*> torsionMap;

    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      torsionMap = pMolecule->getTorsionMap();
      if (!torsionMap.empty()) {
        for (TorsionMapIterator b = torsionMap.begin(); b != torsionMap.end(); b++) {
          pTorsion = b->second;
          std::vector<torsionParam*> torsionParamList = pParameters->getTorsionParamList(
            pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
            pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);
          if (!torsionParamList.empty()) {
///
            if (pTorsion->atom1->getParent()->getName() == "PRO") {
              pParameters->removeProlineTorsion(torsionParamList);
            }
///
            nTorsions += torsionParamList.size();
          }
        }
      }
    }
    return nTorsions;
}

// ============================================================
// Function : getMMTorsions()
// ------------------------------------------------------------
//
// ============================================================
int collection::getMMTorsions(int torsions[], double torsionParams[])
{
    int i = 0;
    int j = 0;
    std::map<ULONG_KIND, Torsion*> torsionMap;

    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    try {
      for (moleculeIterator m = this->itsMoleculeList.begin();
            m != this->itsMoleculeList.end(); m++) {
        pMolecule = *m;
        torsionMap = pMolecule->getTorsionMap();
        if (!torsionMap.empty()) {
          for (TorsionMapIterator b = torsionMap.begin();
               b != torsionMap.end(); b++) {
            pTorsion = b->second;
            std::vector<torsionParam*> torsionParamList =
              pParameters->getTorsionParamList(
              pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
              pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

            if (!torsionParamList.empty()) {
///
              if (pTorsion->atom1->getParent()->getName() == "PRO") {
                pParameters->removeProlineTorsion(torsionParamList);
              }
///
              for (torsionParamIterator c = torsionParamList.begin();
                    c != torsionParamList.end(); c++) {
                pTorsionParam = *c;
                torsions[i]   = pTorsion->atom1->getColIndex()-1;
                torsions[i+1] = pTorsion->atom2->getColIndex()-1;
                torsions[i+2] = pTorsion->atom3->getColIndex()-1;
                torsions[i+3] = pTorsion->atom4->getColIndex()-1;

                torsionParams[j]   = pTorsionParam->Vn / pTorsionParam->npth;
                torsionParams[j+1] = pTorsionParam->Nt;
                torsionParams[j+2] = pTorsionParam->gamma;
                i+=4;
                j+=3;
              }
            }
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }

    return 0;
}

// ============================================================
// Function : getNumImpropers()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumImpropers()
{
    int nImpropers = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      nImpropers += pMolecule->numImpropers();
    }
    return nImpropers;
}

// ============================================================
// Function : getNumMMImpropers()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumMMImpropers()
{
    int nImpropers = 0;
    std::map<int, Improper*> improperMap;
    std::vector<std::vector<int> > order;

    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      improperMap = pMolecule->getImproperMap();
      if (!improperMap.empty()) {
        for (ImproperMapIterator b = improperMap.begin();
             b != improperMap.end(); b++) {
          pImproper = b->second;
          std::vector<improperParam*> improperParamList =
            pParameters->getImproperParamList(
            pImproper->atom1->getStdAtom()->type,
            pImproper->atom2->getStdAtom()->type,
            pImproper->atom3->getStdAtom()->type,
            pImproper->atom4->getStdAtom()->type, order);
          if (!improperParamList.empty()) {
            nImpropers += improperParamList.size();
          }
        }
      }
    }
    return nImpropers;
}

// ============================================================
// Function : getMMImpropers()
// ------------------------------------------------------------
//
// ============================================================
int collection::getMMImpropers(int impropers[], double improperParams[])
{
    int i = 0;
    int j = 0;
    int index = 0;
    std::vector<std::vector<int> > order;
    std::map<int, Improper*> improperMap;
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    try {
      for (moleculeIterator m = this->itsMoleculeList.begin();
            m != this->itsMoleculeList.end(); m++) {
        pMolecule = *m;
        improperMap = pMolecule->getImproperMap();

        if (!improperMap.empty()) {
          for (ImproperMapIterator b = improperMap.begin(); b != improperMap.end(); b++) {
            pImproper = b->second;
            index = 0;
            std::vector<improperParam*> improperParamList = pParameters->getImproperParamList(
              pImproper->atom1->getStdAtom()->type, pImproper->atom2->getStdAtom()->type,
              pImproper->atom3->getStdAtom()->type, pImproper->atom4->getStdAtom()->type, order);

            if (!improperParamList.empty()) {
              atom* tempImproper[4] = {pImproper->atom1, pImproper->atom2, pImproper->atom3, pImproper->atom4};
              for (improperParamIterator c = improperParamList.begin(); c != improperParamList.end(); c++) {
                pImproperParam = *c;
                impropers[i]   = tempImproper[order[index][0]]->getIndex()-1;
                impropers[i+1] = tempImproper[order[index][1]]->getIndex()-1;
                impropers[i+2] = tempImproper[order[index][2]]->getIndex()-1;
                impropers[i+3] = tempImproper[order[index][3]]->getIndex()-1;
                i+=4;
                improperParams[j]   = pImproperParam->Vn;
                improperParams[j+1] = pImproperParam->Nt;
                improperParams[j+2] = pImproperParam->gamma;
                j+=3;
                index++;
              }
            }
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumDihedralsWithH()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumDihedralsWithH()
{
    int nDihedralsWithH = 0;

    std::map<ULONG_KIND, Torsion*> torsionMap;

    typedef std::vector<torsionParam*>::iterator torsionParamIterator;
    for (moleculeIterator m = this->itsMoleculeList.begin();
          m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;
      torsionMap = pMolecule->getTorsionMap();
      if (!torsionMap.empty()) {
        for (TorsionMapIterator b = torsionMap.begin();
             b != torsionMap.end(); b++) {
          pTorsion = b->second;
          std::vector<torsionParam*> torsionParamList =
          pParameters->getTorsionParamList(
          pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
          pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

          if (!torsionParamList.empty()) {
            for (torsionParamIterator c = torsionParamList.begin();
                 c != torsionParamList.end(); c++) {
              pTorsionParam = *c;

              if ((pTorsion->atom1->getElementSymbol() == "H") or
                  (pTorsion->atom4->getElementSymbol() == "H")) {
                nDihedralsWithH++;
              }
            }
          }
        }
      }
    }

    std::vector<std::vector<int> > order;
    std::map<int, Improper*> improperMap;
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    for (moleculeIterator m = this->itsMoleculeList.begin();
         m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;
      improperMap = pMolecule->getImproperMap();

      if (!improperMap.empty()) {
        for (ImproperMapIterator b = improperMap.begin(); b != improperMap.end(); b++) {
          pImproper = b->second;
          std::vector<improperParam*> improperParamList = pParameters->getImproperParamList(
          pImproper->atom1->getStdAtom()->type, pImproper->atom2->getStdAtom()->type,
          pImproper->atom3->getStdAtom()->type, pImproper->atom4->getStdAtom()->type, order);

          if (!improperParamList.empty()) {
            atom* tempImproper[4] = {pImproper->atom1, pImproper->atom2, pImproper->atom3, pImproper->atom4};
            int index = 0;
            for (improperParamIterator c = improperParamList.begin(); c != improperParamList.end(); c++) {
              pImproperParam = *c;

              if ((tempImproper[order[index][0]]->getElementSymbol() == "H") or
                  (tempImproper[order[index][1]]->getElementSymbol() == "H") or
                  (tempImproper[order[index][2]]->getElementSymbol() == "H") or
                  (tempImproper[order[index][3]]->getElementSymbol() == "H")) {
                nDihedralsWithH++;
              }
              index++;
            }
          }
        }
      }
    }
    return nDihedralsWithH;
}

// ============================================================
// Function : getNumDihedralsWithOutH()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumDihedralsWithOutH()
{
    int nDihedralsWithOutH = 0;

    std::map<ULONG_KIND, Torsion*> torsionMap;
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    std::vector<torsionParam*> toBeDeleted;
    std::vector<torsionParam*>::iterator vecIterator;

    for (moleculeIterator m = this->itsMoleculeList.begin();
          m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;
      torsionMap = pMolecule->getTorsionMap();
      if (!torsionMap.empty()) {
        for (TorsionMapIterator b = torsionMap.begin();
             b != torsionMap.end(); b++) {
          pTorsion = b->second;

          std::vector<torsionParam*> torsionParamList =
          pParameters->getTorsionParamList(
          pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
          pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

          if (!torsionParamList.empty()) {
///
            if (pTorsion->atom1->getParent()->getName() == "PRO") {
              pParameters->removeProlineTorsion(torsionParamList);
            }
///
            for (torsionParamIterator c = torsionParamList.begin();
                 c != torsionParamList.end(); c++) {
              pTorsionParam = *c;

              if ((pTorsion->atom1->getElementSymbol() != "H") and
                  (pTorsion->atom4->getElementSymbol() != "H")) {
                nDihedralsWithOutH++;
              }
            }
          }
        }
      }
    }

    int index = 0;
    std::vector<std::vector<int> > order;
    std::map<int, Improper*> improperMap;
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    for (moleculeIterator m = this->itsMoleculeList.begin();
         m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;
      improperMap = pMolecule->getImproperMap();

      if (!improperMap.empty()) {
        for (ImproperMapIterator b = improperMap.begin(); b != improperMap.end(); b++) {
          pImproper = b->second;
          index = 0;
          std::vector<improperParam*> improperParamList = pParameters->getImproperParamList(
          pImproper->atom1->getStdAtom()->type, pImproper->atom2->getStdAtom()->type,
          pImproper->atom3->getStdAtom()->type, pImproper->atom4->getStdAtom()->type, order);

          if (!improperParamList.empty()) {
            atom* tempImproper[4] = {pImproper->atom1, pImproper->atom2, pImproper->atom3, pImproper->atom4};
            for (improperParamIterator c = improperParamList.begin(); c != improperParamList.end(); c++) {
              pImproperParam = *c;

              if ((tempImproper[order[index][0]]->getElementSymbol() != "H") and
                  (tempImproper[order[index][1]]->getElementSymbol() != "H") and
                  (tempImproper[order[index][2]]->getElementSymbol() != "H") and
                  (tempImproper[order[index][3]]->getElementSymbol() != "H")) {
                nDihedralsWithOutH++;
              }
            }
          }
        }
      }
    }
    return nDihedralsWithOutH;
}

// ============================================================
// Function : getNumUniqueDihedralTypes()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumUniqueDihedralTypes()
{
    std::map<ULONG_KIND, Torsion*> torsionMap;
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    for (moleculeIterator m = this->itsMoleculeList.begin();
          m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;
      torsionMap = pMolecule->getTorsionMap();
      if (!torsionMap.empty()) {
        for (TorsionMapIterator b = torsionMap.begin();
             b != torsionMap.end(); b++) {
          pTorsion = b->second;
/*
          std::cout << pTorsion->atom1->getName() << "-"
                    << pTorsion->atom2->getName() << "-"
                    << pTorsion->atom3->getName() << "-"
                    << pTorsion->atom4->getName() << std::endl;
*/
          std::vector<torsionParam*> torsionParamList =
          pParameters->getTorsionParamList(
          pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
          pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

          if (!torsionParamList.empty()) {
///
            if (pTorsion->atom1->getParent()->getName() == "PRO") {
              pParameters->removeProlineTorsion(torsionParamList);
            }
///
            for (torsionParamIterator c = torsionParamList.begin();
                 c != torsionParamList.end(); c++) {
              pTorsionParam = *c;
              torsionParamIterator r = std::find(uniqueTorsionParams.begin(), uniqueTorsionParams.end(), pTorsionParam);
              if (r == uniqueTorsionParams.end()) {
                uniqueTorsionParams.push_back(pTorsionParam);
/*
                std::cout << "   UNIQUE TORSION: " << pTorsionParam->atomType1 << "-"
                          << pTorsionParam->atomType2 << "-"
                          << pTorsionParam->atomType3 << "-"
                          << pTorsionParam->atomType4 << " Nt="
                          << pTorsionParam->Nt << " Vn=" << pTorsionParam->Vn << " Gamma="
                          << pTorsionParam->gamma << " npth=" << pTorsionParam->npth << std::endl;
*/
              }
            }
          }
        }
      }
    }

    int index = 0;
    std::vector<std::vector<int> > order;
    std::map<int, Improper*> improperMap;
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    for (moleculeIterator m = this->itsMoleculeList.begin();
         m != this->itsMoleculeList.end(); m++) {
      pMolecule = *m;
      improperMap = pMolecule->getImproperMap();

      if (!improperMap.empty()) {
        for (ImproperMapIterator b = improperMap.begin(); b != improperMap.end(); b++) {
          pImproper = b->second;
          index = 0;
/*
          std::cout << pImproper->atom1->getName() << "-"
                          << pImproper->atom2->getName() << "-"
                          << pImproper->atom3->getName() << "-"
                          << pImproper->atom4->getName() << std::endl;
*/
          std::vector<improperParam*> improperParamList = pParameters->getImproperParamList(
          pImproper->atom1->getStdAtom()->type, pImproper->atom2->getStdAtom()->type,
          pImproper->atom3->getStdAtom()->type, pImproper->atom4->getStdAtom()->type, order);

          if (!improperParamList.empty()) {
            for (improperParamIterator c = improperParamList.begin(); c != improperParamList.end(); c++) {
              pImproperParam = *c;

              pImproperParam = *c;
              improperParamIterator r = std::find(uniqueImproperParams.begin(), uniqueImproperParams.end(), pImproperParam);
              if (r == uniqueImproperParams.end()) {
                uniqueImproperParams.push_back(pImproperParam);
/*
                std::cout << " UNIQUE IMPROPER: " << pImproperParam->atomType1 << "-"
                          << pImproperParam->atomType2 << "-"
                          << pImproperParam->atomType3 << "-"
                          << pImproperParam->atomType4 << " "
                          << pImproperParam->Nt << " " << pImproperParam->Vn << " "
                          << pImproperParam->gamma << std::endl;
*/
              }
            }
          }
        }
      }
    }
    return static_cast<int>(uniqueTorsionParams.size()) + static_cast<int>(uniqueImproperParams.size());
}

// ============================================================
// Function : getCoordinates()
// ------------------------------------------------------------
//
// ============================================================
int collection::getCoordinates(double coords[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            coords[i] = atList[k]->getX();
            coords[i+1] = atList[k]->getY();
            coords[i+2] = atList[k]->getZ();
            i+=3;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : setCoordinates()
// ------------------------------------------------------------
//
// ============================================================
int collection::setCoordinates(double coords[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            atList[k]->setCoords(coords[i], coords[i+1], coords[i+2]);
            i+=3;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getMMCharges()
// ------------------------------------------------------------
//
// ============================================================
int collection::getMMCharges(double charges[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            if (atList[k]->getStdAtom()) {
              charges[i] = atList[k]->getStdAtom()->atmCharge * E2KCAL;
            }
            else {
              charges[i] = 0.0;
            }
            i++;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getAtomSymbols()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAtomSymbols(char symbols[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            std::string e = atList[k]->getElementSymbol();
            if (e.size() == 2) {
              symbols[i] = e[0];
              symbols[i+1] = e[1];
            }
            else {
              symbols[i] = e[0];
              symbols[i+1] = ' ';
            }
            i+=2;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getAtomNames()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAtomNames(char names[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            std::string e = atList[k]->getName();
            if (e.size() == 4) {
              names[i]   = e[0];
              names[i+1] = e[1];
              names[i+2] = e[2];
              names[i+3] = e[3];
            }
            i+=4;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getAtomMasses()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAtomMasses(double masses[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            element* pEl = atList[k]->getElement();
            if (pEl) {
              masses[i] = pEl->mass;
            }
            else {
              masses[i] = 0.0;
            }
            i++;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumUniqueAtomTypes()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumUniqueAtomTypes()
{
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          stdAtom* pStdA = atList[k]->getStdAtom();
          if (pStdA) {
            std::string atType = pStdA->type;
            //std::cout << pStdA->identity << " " << pStdA->type << std::endl;
            stringVectorIterator r = std::find(atomTypesUsed.begin(), atomTypesUsed.end(), atType);
            if (r == atomTypesUsed.end()) {
              atomTypesUsed.push_back(atType);
              //std::cout << "           Unique Atom Type " << std::endl;
            }
          }
        }
      }
    }
    return static_cast<int>(atomTypesUsed.size());
}

// ============================================================
// Function : getUniqueAtomTypes()
// ------------------------------------------------------------
//
// ============================================================
std::vector<std::string> collection::getUniqueAtomTypes()
{
    return this->atomTypesUsed;
}

// ============================================================
// Function : getUniqueBondTypes()
// ------------------------------------------------------------
//
// ============================================================
std::vector<bondParam*> collection::getUniqueBondTypes()
{
    return this->uniqueBondParams;
}

// ============================================================
// Function : getUniqueAngleTypes()
// ------------------------------------------------------------
//
// ============================================================
std::vector<angleParam*> collection::getUniqueAngleTypes()
{
    return this->uniqueAngleParams;
}

// ============================================================
// Function : getUniqueTorsionTypes()
// ------------------------------------------------------------
//
// ============================================================
std::vector<torsionParam*> collection::getUniqueTorsionTypes()
{
    return this->uniqueTorsionParams;
}

// ============================================================
// Function : getUniqueImproperTypes()
// ------------------------------------------------------------
//
// ============================================================
std::vector<improperParam*> collection::getUniqueImproperTypes()
{
    return this->uniqueImproperParams;
}

// ============================================================
// Function : getAtomTypes()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAtomTypes(int atomTypes[])
{
    int atomIndex = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          stdAtom* pStdA = atList[k]->getStdAtom();
          if (pStdA) {
            std::string atType = pStdA->type;
            bool bFound = false;
            for (unsigned int i = 0; i < atomTypesUsed.size(); i++) {
              if (atomTypesUsed[i] == atType) {
                atomTypes[atomIndex] = i;
                atomIndex++;
                bFound = true;
              }
            }
            if (!bFound) {
              std::cout << " No Type Found ... Error " << std::endl;
              //exit(1);
              throw MTKException(" No Type Found ... Error ");
            }
          }
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getAtomTypes()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAtomTypes(char atomTypes[])
{
    int atomIndex = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          stdAtom* pStdA = atList[k]->getStdAtom();
          if (pStdA) {
            std::string atType = pStdA->type;

            bool bFound = false;
            for (unsigned int i = 0; i < atomTypesUsed.size(); i++) {
              if (atomTypesUsed[i] == atType) {
                if (pStdA->type.size() == 2) {
                  atomTypes[i] = pStdA->type[0];
                  atomTypes[i+1] = pStdA->type[1];
                }
                else {
                  atomTypes[i] = pStdA->type[0];
                  atomTypes[i+1] = ' ';
                }
                atomIndex+=2;
                bFound = true;
              }
            }
            if (!bFound) {
              std::cout << " No Type Found ... Error " << std::endl;
              //exit(1);
              throw MTKException(" No Type Found ... Error ");
            }
          }
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getLJParams()
// ------------------------------------------------------------
//
// ============================================================
int collection::getLJParams(double r6[], double r12[])
{
    errorLogger.throwError("collection::getLJParams", " ", 4);

    unsigned int n = atomTypesUsed.size();

    int indices[n*n];
    int l = 0;
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < n; j++) {
        if (i >= j) {
          indices[j * n + i] = l;
          indices[i * n + j] = l;
          l++;
        }
      }
    }

    int counter = 0;
    for (unsigned int i = 0; i < n; i++) {
      std::string atType1 = atomTypesUsed[i];
      for (unsigned int j = i; j < n; j++) {
        std::string atType2 = atomTypesUsed[j];
        pLJ612SE = pParameters->getLJ612SE(atType1, atType2);
        if (pLJ612SE) {
          double lR6 = 2.0 * pLJ612SE->epsilon * pLJ612SE->sigma;
          double lR12 = pLJ612SE->epsilon * pLJ612SE->sigma * pLJ612SE->sigma;
          int paramIndex = indices[i * n + j];
          r6[paramIndex] = lR6;
          r12[paramIndex] = lR12;
          counter++;
        }
        else {
          std::string errMessage = " Can't find Parameters for " + atType1 + " " + atType2;
          errorLogger.throwError("collection::getLJParams", errMessage, 1);
          //exit(1);
          std::stringstream ss;
          ss << "collection::getLJParams"<< errMessage;
          throw MTKException(ss.str());
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getNumExcludedAtoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumExcludedAtoms()
{
    std::vector<atom*> atomList = this->getAtomList();
    unsigned int nAtoms = atomList.size();
    int excludedAtoms = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if (pAt1->hasBondedAtom(pAt2) or
            pAt1->has13BondedAtom(pAt2) or
            pAt1->has14BondedAtom(pAt2)) {
          excludedAtoms++;
        }
      }
    }

/*
    int nE = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          nE = 0;
          nE = atList[k]->getNumBonds();
          nE += atList[k]->getNum13Bonds();
          nE += atList[k]->getNum14Bonds();
          excludedAtoms += nE;
        }
      }
    }
*/
    std::stringstream number;
    number << excludedAtoms;
    std::string errMessage = " Number of Excluded Atoms = " + number.str();
    MTKpp::errorLogger.throwError("collection::getNumExcludedAtoms", errMessage, 4);

    return excludedAtoms;
}

// ============================================================
// Function : getNumExcluded14Atoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumExcluded14Atoms()
{
    std::vector<atom*> atomList = this->getAtomList();
    unsigned int nAtoms = atomList.size();
    int excluded14Atoms = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if ((pAt1->has14BondedAtom(pAt2)) and
            (!pAt1->has13BondedAtom(pAt2))) {
          excluded14Atoms++;
        }
      }
    }

/*
    int nE = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          nE = atList[k]->getNum14Bonds();
          excluded14Atoms += nE;
        }
      }
    }
*/
    std::stringstream number;
    number << excluded14Atoms;
    std::string errMessage = " Number of Excluded 1-4 Atoms = " + number.str();
    MTKpp::errorLogger.throwError("collection::getNumExcluded14Atoms", errMessage, 4);

    return excluded14Atoms;
}

// ============================================================
// Function : getNumExcludedAtoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumExcludedAtoms(int excludedAtoms[])
{
    std::vector<atom*> atomList = this->getAtomList();
    unsigned int nAtoms = atomList.size();
    int totalExcludedAtoms = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];

      int nE = 0;
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if (pAt1->hasBondedAtom(pAt2) or
            pAt1->has13BondedAtom(pAt2) or
            pAt1->has14BondedAtom(pAt2)) {
          nE++;
        }
      }
      excludedAtoms[i] = nE;
      totalExcludedAtoms+=nE;
    }

/*
    int eAtoms = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          eAtoms = 0;
          eAtoms = atList[k]->getNumBonds();
          eAtoms += atList[k]->getNum13Bonds();
          eAtoms += atList[k]->getNum14Bonds();
          excludedAtoms[atIndex] = eAtoms;
          totalExcludedAtoms += eAtoms;
          atIndex++;
        }
      }
    }
*/
    std::stringstream number;
    number << totalExcludedAtoms;
    std::string errMessage = " Number of Excluded Atoms = " + number.str();
    MTKpp::errorLogger.throwError("collection::getNumExcludedAtoms", errMessage, 4);

    return 0;
}

// ============================================================
// Function : getNumExcluded14Atoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumExcluded14Atoms(int excluded14Atoms[])
{
    std::vector<atom*> atomList = this->getAtomList();
    unsigned int nAtoms = atomList.size();
    int totalExcluded14Atoms = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];

      int nE = 0;
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if ((pAt1->has14BondedAtom(pAt2)) and
            (!pAt1->has13BondedAtom(pAt2))) {
          nE++;
        }
      }
      excluded14Atoms[i] = nE;
      totalExcluded14Atoms+=nE;
    }

/*
    int atIndex = 0;
    int eAtoms = 0;
    int totalExcluded14Atoms = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          eAtoms = atList[k]->getNum14Bonds();
          excluded14Atoms[atIndex] = eAtoms;
          totalExcluded14Atoms += eAtoms;
          atIndex++;
        }
      }
    }
*/
    std::stringstream number;
    number << totalExcluded14Atoms;
    std::string errMessage = " Number of Excluded 1-4 Atoms = " + number.str();
    MTKpp::errorLogger.throwError("collection::getNumExcluded14Atoms", errMessage, 4);

    return 0;
}

// ============================================================
// Function : getExcludedAtoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getExcludedAtoms(int excludedAtoms[])
{
    std::vector<atom*> atomList = this->getAtomList();
    unsigned int nAtoms = atomList.size();
    int atIndex = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];

      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if (pAt1->hasBondedAtom(pAt2) or
            pAt1->has13BondedAtom(pAt2) or
            pAt1->has14BondedAtom(pAt2)) {
          excludedAtoms[atIndex] = pAt2->getColIndex()-1;
          atIndex++;
        }
      }
    }

/*
    int eAtoms = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        eAtoms = 0;
        for (unsigned int k = 0; k < atList.size(); k++) {
          std::vector<int> exAtoms;
          std::vector<atom*> bondedAtoms = atList[k]->getBondedAtoms();
          for (unsigned int l = 0; l < bondedAtoms.size(); l++) {
            exAtoms.push_back(bondedAtoms[l]->getColIndex()-1);
          }
          std::vector<atom*> bonded13Atoms = atList[k]->get13BondedAtoms();
          for (unsigned int l = 0; l < bonded13Atoms.size(); l++) {
            exAtoms.push_back(bonded13Atoms[l]->getColIndex()-1);
          }
          std::vector<atom*> bonded14Atoms = atList[k]->get14BondedAtoms();
          for (unsigned int l = 0; l < bonded14Atoms.size(); l++) {
            exAtoms.push_back(bonded14Atoms[l]->getColIndex()-1);
          }
          std::sort(exAtoms.begin(), exAtoms.end());
          for (unsigned int l = 0; l < exAtoms.size(); l++) {
            excludedAtoms[atIndex] = exAtoms[l];
            atIndex++;
          }
        }
      }
    }
*/
    return 0;
}

// ============================================================
// Function : getExcluded14Atoms()
// ------------------------------------------------------------
//
// ============================================================
int collection::getExcluded14Atoms(int excluded14Atoms[])
{
    std::vector<atom*> atomList = this->getAtomList();
    unsigned int nAtoms = atomList.size();
    int atIndex = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];

      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if ((pAt1->has14BondedAtom(pAt2)) and
            (!pAt1->has13BondedAtom(pAt2))) {
          excluded14Atoms[atIndex] = pAt2->getColIndex()-1;
          atIndex++;
        }
      }
    }

/*
    int atIndex = 0;
    int eAtoms = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        eAtoms = 0;
        for (unsigned int k = 0; k < atList.size(); k++) {
          std::vector<int> exAtoms;

          std::vector<atom*> bonded14Atoms = atList[k]->get14BondedAtoms();
          for (unsigned int l = 0; l < bonded14Atoms.size(); l++) {
            exAtoms.push_back(bonded14Atoms[l]->getColIndex()-1);
          }
          std::sort(exAtoms.begin(), exAtoms.end());
          for (unsigned int l = 0; l < exAtoms.size(); l++) {
            excluded14Atoms[atIndex] = exAtoms[l];
            atIndex++;
          }
        }
      }
    }
*/
    return 0;
}

// ============================================================
// Function : getResidueNames()
// ------------------------------------------------------------
//
// ============================================================
int collection::getResidueNames(char resNames[])
{
    int i = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          std::string e = subMolList[j]->getName();
          if (e.size() == 3) {
            resNames[i]   = e[0];
            resNames[i+1] = e[1];
            resNames[i+2] = e[2];
          }
          else {
            std::cout << " ERROR ... exiting " << std::endl;
            //exit(0);
            throw MTKException(" ERROR ... exiting ");
          }
          i+=3;
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getResiduePointers()
// ------------------------------------------------------------
//
// ============================================================
int collection::getResiduePointers(int resPointers[])
{
    int resIndex = 0;
    int atomIndex = 0;
    try {
      for (moleculeIterator c = this->itsMoleculeList.begin();
            c != this->itsMoleculeList.end(); c++) {
        pMolecule = *c;
        std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
        for (unsigned int j = 0; j < subMolList.size(); j++) {
          resPointers[resIndex] = atomIndex;
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            atomIndex++;
          }
          resIndex++;
        }
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Out of bounds Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : setAtomIndex()
// ------------------------------------------------------------
// Used in submolecule when an atom is added
// ============================================================
void collection::setAtomIndex(const int &n)
{
    this->itsAtomIndex = n;
}

// ============================================================
// Function : getAtomIndex()
// ------------------------------------------------------------
//
// ============================================================
int collection::getAtomIndex()
{
    return this->itsAtomIndex;
}

// ============================================================
// Function : setSubMoleculeIndex()
// ------------------------------------------------------------
// Used in molecule when a submolecule is added
// ============================================================
void collection::setSubMoleculeIndex(const int &n)
{
    this->itsSubMoleculeIndex = n;
}

// ============================================================
// Function : getSubMoleculeIndex()
// ------------------------------------------------------------
//
// ============================================================
int collection::getSubMoleculeIndex()
{
    return this->itsSubMoleculeIndex;
}

// ============================================================
// Function : getFormalCharge()
// ------------------------------------------------------------
//
// ============================================================
int collection::getFormalCharge()
{
    typedef boost::numeric::converter<int, double> Double2Int;
    double formalCharge = 0.0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        stdFrag* pStdFrag = subMolList[j]->getStdFrag();
        if (pStdFrag) {
          std::vector<atom*> atomList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atomList.size(); k++) {
            if (atomList[k]->getStdAtom()) {
              formalCharge += atomList[k]->getStdAtom()->atmCharge;
            }
            else {
              std::cout << " Error " << std::endl;
            }
          }
          //formalCharge += pStdFrag->getCharge();
        }
        else {
          std::vector<atom*> atList = subMolList[j]->getAtomList();
          for (unsigned int k = 0; k < atList.size(); k++) {
            formalCharge += atList[k]->getFormalCharge();
          }
        }
      }
    }
    int iCharge2 = static_cast<int>((formalCharge < double(0)) ? (formalCharge - double(.5)) : (formalCharge + double(.5)));

    int iCharge = Double2Int()(formalCharge);
    if (formalCharge > 0) {
      if (iCharge + 0.5 < formalCharge) iCharge++;
    }
    else {
      if (iCharge - 0.5 > formalCharge) iCharge--;
    }
    std::cout << " formal Charge = " << formalCharge << " iCharge = " << iCharge << "   iCharge2 = " << iCharge2 << std::endl;
    return iCharge;
}

// ============================================================
// Function : getNumNeighbors()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumNeighbors(atom* pAt, double distance)
{
    submolecule* pAtSubMol = pAt->getParent();
    int nNeighbors = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        if (pAtSubMol == subMolList[j]) continue;
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          if (!pAt->hasBondedAtom(atList[k]) or
              !pAt->has13BondedAtom(atList[k]) or
              !pAt->has14BondedAtom(atList[k])) {
            if (atList[k]->getCoords()->dist((*pAt->getCoords())) < distance) {
              nNeighbors++;
            }
          }
        }
      }
    }
    return nNeighbors;
}

// ============================================================
// Function : getNumHeavyNeighbors()
// ------------------------------------------------------------
//
// ============================================================
int collection::getNumHeavyNeighbors(atom* pAt, double distance)
{
    submolecule* pAtSubMol = pAt->getParent();
    int nNeighbors = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        if (pAtSubMol == subMolList[j]) continue;
        std::vector<atom*> atList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atList.size(); k++) {
          if (atList[k]->getElementSymbol() == "H") continue;
          if (!pAt->hasBondedAtom(atList[k]) or
              !pAt->has13BondedAtom(atList[k]) or
              !pAt->has14BondedAtom(atList[k])) {
            if (atList[k]->getCoords()->dist((*pAt->getCoords())) < distance) {
              nNeighbors++;
            }
          }
        }
      }
    }
    return nNeighbors;
}

// ============================================================
// Function : largestResidueSize()
// ------------------------------------------------------------
//
// ============================================================
int collection::largestResidueSize()
{
    int maxRes = 0;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        int nAtoms = subMolList[j]->getNumAtoms();
        if (nAtoms > maxRes) {
          maxRes = nAtoms;
        }
      }
    }
    return maxRes;
}

// ============================================================
// Function : renumber()
// ------------------------------------------------------------
//
// ============================================================
void collection::renumber()
{
    int nRes = 1;
    int nAtom = 1;
    for (moleculeIterator c = this->itsMoleculeList.begin();
          c != this->itsMoleculeList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> subMolList = pMolecule->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        subMolList[j]->setSubMolId(nRes);
        nRes++;
        std::vector<atom*> atomList = subMolList[j]->getAtomList();
        for (unsigned int k = 0; k < atomList.size(); k++) {
          atomList[k]->setFileID(nAtom);
          nAtom++;
        }
      }
    }
}

    /////////////////////////////
    // -   Metalloproteins   - //
    /////////////////////////////

// ============================================================
// Function : hasMetal()
// ------------------------------------------------------------
//
// ============================================================
bool collection::hasMetal()
{
    for (unsigned int m = 0; m < this->itsMoleculeList.size(); m++) {
      std::vector<atom*> atomList = this->itsMoleculeList[m]->getAtomList();
      for (unsigned int a = 0; a < atomList.size(); a++) {
        std::string currentElement = atomList[a]->getElementSymbol();
        stringVectorIterator r = std::find(availableMetals.begin(), availableMetals.end(), currentElement);
        if (r != availableMetals.end()) {
          MTKpp::errorLogger.throwError("collection::hasMetal", "\n   Found Metal: "+currentElement, INFO);
          return true;
        }
      }
    }
    return false;
}

// ============================================================
// Function : findMetals()
// ------------------------------------------------------------
//
// ============================================================
void collection::findMetals()
{
    std::string errMessage = "\n Metal Found: ";
    for (unsigned int m = 0; m < this->itsMoleculeList.size(); m++) {
      std::vector<atom*> atomList = this->itsMoleculeList[m]->getAtomList();
      if (this->itsMoleculeList[m]->getName() == "Reference") continue;
      for (unsigned int a = 0; a < atomList.size(); a++) {
        std::string currentElement = atomList[a]->getElementSymbol();
        stringVectorIterator r = std::find(availableMetals.begin(), availableMetals.end(), currentElement);
        if (r != availableMetals.end()) {
          this->itsMetalAtoms.push_back(atomList[a]);
          errMessage += "\n    " + atomList[a]->getElementSymbol() +
                        ":" + i2s(atomList[a]->getFileID());
        }
      }
    }
    //errMessage += "\n";
    MTKpp::errorLogger.throwError("collection::findMetals", errMessage, INFO);
}

// ============================================================
// Function : determineMetalEnvironments()
// ------------------------------------------------------------
//
// ============================================================
void collection::determineMetalEnvironments()
{
    std::string errMessage = "\n Metal Environments: \n";
    typedef std::map<stdAtom*, atom*>::iterator funcGroupAtomIterator;
    for (unsigned int m = 0; m < this->itsMetalAtoms.size(); m++) {

      metalCenter* metCen = new metalCenter(this->itsMetalAtoms[m], metalDonorDists);
      vector3d* metalCoords = this->itsMetalAtoms[m]->getCoords();
      std::string cMetal = this->itsMetalAtoms[m]->getElementSymbol();

      if (this->itsMetalAtoms[m]->getParent()->getParent()->getName() == "Reference") continue;

      //errMessage += "\n" + cMetal + " " + i2s(this->itsMetalAtoms[m]->getFileID());

      for (unsigned int n = 0; n < this->itsMoleculeList.size(); n++) {
        if (this->itsMoleculeList[n]->getName() == "Reference") continue;
        std::vector<atom*> atomList = this->itsMoleculeList[n]->getAtomList();
        for (unsigned int a = 0; a < atomList.size(); a++) {
          std::string atSymbol = atomList[a]->getElementSymbol();
          if (atomList[a] == this->itsMetalAtoms[m]) continue;
          if (atSymbol == "C" or atSymbol == "H" or atSymbol == "P") continue;
          double myDistance = metalCoords->dist( (*atomList[a]->getCoords()) );
          if (myDistance < 0.5) {
            metCen->setError(5);
          }
          else if (myDistance < 3.6) {
            std::string resName = atomList[a]->getParent()->getName();

            if (atomList[a]->getParent()->getStdFrag()) {
              resName = atomList[a]->getParent()->getStdFrag()->getCharacter();
            }

            std::string atName = atomList[a]->getName();
            std::string label = cMetal + ":" + resName;

            if (atName == " O  ") { // CARBONYL
              label = cMetal + ":CRL";
            }

            //std::cout << "   " << label << " " << myDistance << "\n";

            strDbMapIterator b = metalDonorDists.find(label);
            if (b != metalDonorDists.end()) {
              if (myDistance < metalDonorDists[label] + 1.0) {
                metCen->addAtom(atomList[a]);
                metCen->addLabel(label);
                //errMessage += "\n    M-L: " + atomList[a]->getParent()->getName() + "@|" + 
                //              atomList[a]->getName() + "| dist = " + d2s(myDistance) + " [" + label + "]";
              }
            }
            else {
              if (atSymbol == "O") {
                label = cMetal + ":SER";
              }
              else if (atSymbol == "N") {
                label = cMetal + ":HIS";
              }
              else {
                label = cMetal + ":CYS";
              }
              if (atomList[a]->getNumBonds() > 2) {
/*
                std::cout << atomList[a]->getParent()->getName() << " "
                          << atomList[a]->getName() << " has "
                          << atomList[a]->getNumBonds() << " bonds" << std::endl;
*/
                continue;
              }
/*
              std::cout << "|" << atomList[a]->getParent()->getName() << "| |"
                        << atomList[a]->getName() << "| |"
                        << atomList[a]->getParent()->getStdFrag()->getCharacter() << "| |"
                        << label << "| "
                        << myDistance << std::endl;
*/
              strDbMapIterator b2 = metalDonorDists.find(label);
              if (b2 != metalDonorDists.end()) {
                if (myDistance < metalDonorDists[label] + 0.5) {
                  metCen->addAtom(atomList[a]);
                  metCen->addLabel(label);
                  //errMessage += "\n    M-L: " + atomList[a]->getParent()->getName() + "@|" + 
                  //            atomList[a]->getName() + "| dist = " + d2s(myDistance) + " [" + label + "]";
                }
              }
            }
          }
        }
      }

      std::string stdGroupName = cMetal + "-";
      if (!metCen->getError()) {
        metCen->assignBondType();
        metCen->assignGeometry();

        std::string priShell = metCen->getPrimaryShell();
        std::sort(priShell.begin(), priShell.end());
        std::string secShell = metCen->getSecondaryShell();
        std::sort(secShell.begin(), secShell.end());

        std::string g1 = metCen->getGeometry();
        double gr = metCen->getGeometryRMS();
/*
        std::cout << "   collection::determineMetalEnvironments: " << cMetal + "-"
                  << priShell << " " << secShell << std::endl;
*/
        errMessage += "  GROUP:" + cMetal + "-[" + priShell + "][" + secShell + "] geom1 = " + g1
                    + " geom1RMSD = " + d2s(gr) +
                      "\n";
        this->itsMetalCenters.push_back(metCen);
        stdGroupName += priShell;
      }
      else {
        std::cout << " metCen error " << metCen->getError() << "\n";
      }
    }

    // std::cout << "collection::determineMetalEnvironments " << errMessage << std::endl;

    MTKpp::errorLogger.throwError("collection::determineMetalEnvironments", errMessage, INFO);

    this->assignMetalPartners();

    this->assignMetalParameters();
    // std::cout << "collection::determineMetalEnvironments END" << std::endl;
}

// ============================================================
// Function : assignMetalPartners()
// ------------------------------------------------------------
// Determine if metal centers share ligating residues
// ============================================================
void collection::assignMetalPartners()
{
    // std::cout << "collection::assignMetalPartners start " << std::endl;
    std::string errMessage = "\n";

    //
    std::vector<std::set<int> > partners;

    if (this->itsMetalCenters.size() > 1) {
      // std::cout << "nMetalCenters " << this->itsMetalCenters.size() << std::endl;

      for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {

        std::set<int> mmPartner;
        mmPartner.insert(m);
        mmPartner.insert(m);
        partners.push_back(mmPartner);

        metalCenter* metCen1 = this->itsMetalCenters[m];
        // std::cout << "metal center " << this->itsMetalCenters[m] << std::endl;

        // Primary Shell
        std::vector<atom*> primShellAtoms1;
        metCen1->getPrimaryShellAtoms(primShellAtoms1);

        std::vector<submolecule*> primShellResidues1;
        std::vector<molecule*> primShellMolecules1;

        for (unsigned int i = 0; i < primShellAtoms1.size(); i++) {
          primShellResidues1.push_back(primShellAtoms1[i]->getParent());
          primShellMolecules1.push_back(primShellResidues1[i]->getParent());
        }
        // std::cout << "nPrimaryResidues " << primShellResidues1.size() << std::endl;

        // Secondary Shell
        std::vector<atom*> secShellAtoms1;
        metCen1->getSecondaryShellAtoms(secShellAtoms1);

        std::vector<submolecule*> secShellResidues1;
        std::vector<molecule*> secShellMolecules1;

        for (unsigned int i = 0; i < secShellAtoms1.size(); i++) {
          secShellResidues1.push_back(secShellAtoms1[i]->getParent());
          secShellMolecules1.push_back(secShellResidues1[i]->getParent());
        }
        // std::cout << "nSecondaryResidues " << secShellResidues1.size() << std::endl;

        for (unsigned int n = m; n < this->itsMetalCenters.size(); n++) {
          if (n == m) continue;
          metalCenter* metCen2 = this->itsMetalCenters[n];

          // Primary Shell
          std::vector<atom*> primShellAtoms2;
          metCen2->getPrimaryShellAtoms(primShellAtoms2);

          std::vector<submolecule*> primShellResidues2;
          std::vector<molecule*> primShellMolecules2;

          for (unsigned int i = 0; i < primShellAtoms2.size(); i++) {
            primShellResidues2.push_back(primShellAtoms2[i]->getParent());
            primShellMolecules2.push_back(primShellResidues2[i]->getParent());
          }

          // Secondary Shell
          std::vector<atom*> secShellAtoms2;
          metCen2->getSecondaryShellAtoms(secShellAtoms2);

          std::vector<submolecule*> secShellResidues2;
          std::vector<molecule*> secShellMolecules2;

          for (unsigned int i = 0; i < secShellAtoms2.size(); i++) {
            secShellResidues2.push_back(secShellAtoms2[i]->getParent());
            secShellMolecules2.push_back(secShellResidues2[i]->getParent());
          }
          // std::cout << "000" <<std::endl;

          bool bCommonResidue = false;
          for (unsigned int i = 0; i < primShellResidues1.size(); i++) {
            for (unsigned int j = 0; j < primShellResidues2.size(); j++) {
              if (primShellResidues1[i] == primShellResidues2[j]) {
                bCommonResidue = true;
//std::cout << " MATCH1 " << primShellResidues1[i]->getName() << primShellResidues1[i]->getSubMolId() << "\n";
              }
            }
            for (unsigned int j = 0; j < secShellResidues2.size(); j++) {
              if (primShellResidues1[i] == secShellResidues2[j]) {
                bCommonResidue = true;
//std::cout << " MATCH2 " << primShellResidues1[i]->getName() << primShellResidues1[i]->getSubMolId() << "\n";
              }
            }
          }

          if (bCommonResidue) {
            std::set<int> partner;
            partner.insert(m);
            partner.insert(n);
            partners.push_back(partner);
          }
        }
      }
    }
    else {
      std::set<int> mmPartner;
      mmPartner.insert(0);
      mmPartner.insert(0);
      partners.push_back(mmPartner);
    }

/*

      Metal Ligating residues
      1     5-6
      2       6-7
      3         7-8
      4     5-    8

      Group
      Metals: 1,2,3,4
      Residues: 5,6,7,8

      Partners:
       1-2
       2-3
       3-4
       1-4

         1  2  3  4
      1  x  x  -  x
      2     x  x  -
      3        x  x
      4           x

       and/or
         5  6  7  8     5  6  7  8      5  6  7  8
      1  1  1  0  0
      2  0  1  1  0 --> 1  1  1  0 
      3  0  0  1  1                 --> 1  1  1  1
      4  1  0  0  1

      set_difference             computes the difference between two sets
      set_intersection           computes the intersection of two sets
      set_symmetric_difference   computes the symmetric difference between two sets
      set_union                  computes the union of two sets

*/

    metalGroup* pMetalGroup;

    if (partners.size() > 0) {

/*
      for (unsigned int i = 0; i < partners.size(); i++) {
        std::cout << " Partner: " << i << " (";
        for (std::set<int>::const_iterator iter = partners[i].begin(); iter != partners[i].end(); ++iter ) {
          std::cout << *iter << " ";
        }
        std::cout << ")\n";
      }
*/

      for (unsigned int i = 0; i < partners.size(); i++) {
        if (itsMetalGroups.size() == 0) {
          pMetalGroup = new metalGroup();
          itsMetalGroups.push_back(pMetalGroup);
          for (std::set<int>::const_iterator iter = partners[i].begin(); iter != partners[i].end(); ++iter) {
            pMetalGroup->addMetalCenter(this->itsMetalCenters[*iter]);
          }
        }
        else {
          bool bGroup = false;

          for (unsigned int j = 0; j < itsMetalGroups.size(); j++) {
            std::vector<metalCenter*> metCenters = itsMetalGroups[j]->getMetalCenters();

            for (std::set<int>::const_iterator iter = partners[i].begin(); iter != partners[i].end(); ++iter) {
              for (unsigned int k = 0; k < metCenters.size(); k++) {
                if (this->itsMetalCenters[*iter] == metCenters[k]) {
                  bGroup = true;
                }
              }
            }

            if (bGroup) {
              pMetalGroup = itsMetalGroups[j];
              break;
            }
          }

          if (!bGroup) {
            pMetalGroup = new metalGroup();
            itsMetalGroups.push_back(pMetalGroup);
          }

          for (std::set<int>::const_iterator iter = partners[i].begin(); iter != partners[i].end(); ++iter) {
            pMetalGroup->addMetalCenter(this->itsMetalCenters[*iter]);
          }
        }
      }
    }

    errMessage += " Number of Metal Groups = " + i2s(itsMetalGroups.size());
    MTKpp::errorLogger.throwError("collection::assignMetalPartners", errMessage, INFO);
}

// ============================================================
// Function : assignMetalParameters()
// ------------------------------------------------------------
//
// ============================================================
void collection::assignMetalParameters()
{
    std::string errMessage = "";
    // Add primary and secondary bonds
    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {

      metalCenter* metCen = this->itsMetalCenters[m];
      atom* metalAtom = metCen->getMetalAtom();

      errMessage += "\n Metal Center: " + i2s(m+1) + "\n";
      // Primary shell
      std::vector<atom*> primShellAtoms;
      metCen->getPrimaryShellAtoms(primShellAtoms);
      int nBondedAtoms = static_cast<int>(primShellAtoms.size());

      if (nBondedAtoms > 0) errMessage += "    Adding Primary Bonds: \n";

      // Add bonds to metal center
      std::string angleErrMessage = "";
      for (int e = 0; e < nBondedAtoms; e++) {
        Bond* localBond = metCen->addBond(metalAtom, primShellAtoms[e]);
        metalAtom->addBondedAtom(primShellAtoms[e]);
        primShellAtoms[e]->addBondedAtom(metalAtom);
        if (!localBond) {
          std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
          //exit(1);
          std::stringstream ss;
          ss << "collection::assignMetalParameters"<< " Error ";
          throw MTKException(ss.str());
        }
        errMessage += metalAtom->getParent()->getName() + "@|"
                    + metalAtom->getName() + "|-"
                    + primShellAtoms[e]->getParent()->getName()
                    + i2s(primShellAtoms[e]->getParent()->getSubMolId())
                    + primShellAtoms[e]->getParent()->getiCode()
                    + "@|"
                    + primShellAtoms[e]->getName() + "| dist = " + d2s(localBond->size) +"\n";

        // Angles
        for (int d = 0; d < nBondedAtoms; d++) {
          if (e == d) continue;
          Angle* localAngle = metCen->addAngle(primShellAtoms[e], this->itsMetalAtoms[m], primShellAtoms[d]);
          if (!localAngle) {
            std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
          }
          bool added = primShellAtoms[e]->addBonded13Atom(primShellAtoms[d]);
          bool added2 = primShellAtoms[d]->addBonded13Atom(primShellAtoms[e]);
          if (added or added2) {
            angleErrMessage += primShellAtoms[e]->getParent()->getName() + "@|" + primShellAtoms[e]->getName() + "|-"
                             + this->itsMetalAtoms[m]->getParent()->getName() + "@|" + this->itsMetalAtoms[m]->getName() + "|-"
                             + primShellAtoms[d]->getParent()->getName() + "@|" + primShellAtoms[d]->getName() 
                             + "| angle = " + d2s(localAngle->size * RAD2DEG) + "\n";
          }
          //metCen->add13Bond(primShellAtoms[e], primShellAtoms[d]);
        }

        std::vector<atom*> lBondedAtoms = primShellAtoms[e]->getBondedAtoms();
        for (unsigned int d = 0; d < lBondedAtoms.size(); d++) {
          if (this->itsMetalAtoms[m] == lBondedAtoms[d]) continue;
          Angle* localAngle = metCen->addAngle(this->itsMetalAtoms[m], primShellAtoms[e], lBondedAtoms[d]);
          if (!localAngle) {
            std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
          }
          bool added = this->itsMetalAtoms[m]->addBonded13Atom(lBondedAtoms[d]);
          bool added2 = lBondedAtoms[d]->addBonded13Atom(this->itsMetalAtoms[m]);
          if (added or added2) {
            angleErrMessage += this->itsMetalAtoms[m]->getParent()->getName() + "@|" + primShellAtoms[e]->getName() + "|-"
                             + primShellAtoms[e]->getParent()->getName() + "@|" + primShellAtoms[e]->getName() + "|-"
                             + lBondedAtoms[d]->getParent()->getName() + "@|" + lBondedAtoms[d]->getName() 
                             + "| angle = " + d2s(localAngle->size * RAD2DEG) + "\n";
          }
        }
      }

      errMessage += angleErrMessage;

      // Secondary shell
      std::vector<atom*> secShellAtoms;
      metCen->getSecondaryShellAtoms(secShellAtoms);
/*
      int nSecBondedAtoms = static_cast<int>(secShellAtoms.size());

      if (nSecBondedAtoms > 0) errMessage += "  Adding Secondary Bond: \n";

      // Add bonds to metal center
      for (int e = 0; e < nSecBondedAtoms; e++) {
        Bond* localBond = metCen->addBond(metalAtom, secShellAtoms[e]);
        metalAtom->addBondedAtom(secShellAtoms[e]);
        secShellAtoms[e]->addBondedAtom(metalAtom);
        if (!localBond) {
          std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
          exit(1);
        }
        errMessage += metalAtom->getParent()->getName() + "@|"
                    + metalAtom->getName() + "|-"
                    + secShellAtoms[e]->getParent()->getName() + "@|"
                    + secShellAtoms[e]->getName() + "| dist = " + d2s(localBond->size) +"\n";
      }
*/
    }

    errMessage += "\n number of metal atoms = ";
    errMessage += i2s(this->itsMetalAtoms.size());
    errMessage += "\n number of metal centers = ";
    errMessage += i2s(this->itsMetalCenters.size());
    errMessage += "\n number of metal groups = ";
    errMessage += i2s(itsMetalGroups.size());

    MTKpp::errorLogger.throwError("collection::assignMetalParameters", errMessage, INFO);

    if (itsMetalGroups.size() > 0) {
      if (this->pStdLibrary) {
        std::vector<stdGroup*> stdGrpList = this->pStdLibrary->getStdGroupList();

        for (unsigned int i = 0; i < stdGrpList.size(); i++) {
//std::cout << "\n\ncollection::assignMetalParameters " << i << " " << stdGrpList.size() << " " << stdGrpList[i]->getName() << " "
//          << stdGrpList[i]->hasStdMolecule() << std::endl;
          if (stdGrpList[i]->hasStdMolecule()) {
            //std::cout << "stdgroup name: " << stdGrpList[i]->getName() << std::endl;
            //stdGroup* pStdGroup = stdGrpList[i];
            molecule* pStdMolecule = stdGrpList[i]->getStdMolecule();

/*
            std::map<int, Bond*> stdBondMap = pStdMolecule->getBondMap();
            std::cout << " \n\n STANDARD MOLECULE \n";
            if (!stdBondMap.empty()) {
              for (BondMapIterator b = stdBondMap.begin(); b != stdBondMap.end(); b++) {
                Bond* pBond = b->second;
                atom* pBdAtom1 = pBond->atom1;
                atom* pBdAtom2 = pBond->atom2;
                std::cout << pBdAtom1->getParent()->getName() << ":|" << pBdAtom1->getName() << "|-"
                          << pBdAtom2->getParent()->getName() << ":|" << pBdAtom2->getName() << "|\n";
              }
            }
*/

            // loop over metal groups
            for (unsigned int j = 0; j < itsMetalGroups.size(); j++) {

              //std::cout << "\n Metal Group " << j << std::endl;

              // create new molecule
              molecule* pNewMolecule = this->addMolecule();
              std::vector<atom*> atomList;

              // copy residues
              std::vector<metalCenter*> metCens = itsMetalGroups[j]->getMetalCenters();

              std::vector<submolecule*> addedResidues;
              std::vector<submolecule*>::iterator result;

              for (unsigned int k = 0; k < metCens.size(); k++) {
                metalCenter* metCen = metCens[k];

                // Primary Shell
                std::vector<atom*> primShellAtoms;
                metCen->getPrimaryShellAtoms(primShellAtoms);
 
                for (unsigned int l = 0; l < primShellAtoms.size(); l++) {
                  result = std::find(addedResidues.begin(),
                                addedResidues.end(), primShellAtoms[l]->getParent());
                  if (result == addedResidues.end()) {
                    addedResidues.push_back(primShellAtoms[l]->getParent());
                    submolecule* pNewSubmolecule = pNewMolecule->addSubMolecule();
                    pNewSubmolecule->copy(primShellAtoms[l]->getParent());
                    std::vector<atom*> atomListTmp = primShellAtoms[l]->getParent()->getAtomList();
                    for (unsigned int i = 0 ; i < atomListTmp.size(); i++) {
                        atomList.push_back(atomListTmp[i]);
                    }
/*
std::cout << "Primary Shell " << primShellAtoms[l]->getParent()->getName() << " "
          << primShellAtoms[l]->getParent()->getSubMolId() << " "
          << primShellAtoms[l]->getParent()->getNumAtoms() << " "
          << std::endl;
*/
                  }
                }

                // Secondary Shell
                std::vector<atom*> secShellAtoms;
                metCen->getSecondaryShellAtoms(secShellAtoms);

                for (unsigned int l = 0; l < secShellAtoms.size(); l++) {
                  result = std::find(addedResidues.begin(),
                                addedResidues.end(), secShellAtoms[l]->getParent());
                  if (result == addedResidues.end()) {
                    addedResidues.push_back(secShellAtoms[l]->getParent());
                    submolecule* pNewSubmolecule = pNewMolecule->addSubMolecule();
                    pNewSubmolecule->copy(secShellAtoms[l]->getParent());

                    std::vector<atom*> atomListTmp = secShellAtoms[l]->getParent()->getAtomList();
                    for (unsigned int i = 0 ; i < atomListTmp.size(); i++) {
                        atomList.push_back(atomListTmp[i]);
                    }
/*
std::cout << "Secondary Shell " << secShellAtoms[l]->getParent()->getName() << " "
          << secShellAtoms[l]->getParent()->getSubMolId()
          << secShellAtoms[l]->getParent()->getNumAtoms() << " "
          << std::endl;
*/
                  }
                }

                // Metal Atom
                atom* metalAtom = metCen->getMetalAtom();
                result = std::find(addedResidues.begin(),
                                addedResidues.end(), metalAtom->getParent());
                if (result == addedResidues.end()) {
                  addedResidues.push_back(metalAtom->getParent());
                  submolecule* pNewSubmolecule = pNewMolecule->addSubMolecule();
                  pNewSubmolecule->copy(metalAtom->getParent());

                    std::vector<atom*> atomListTmp = metalAtom->getParent()->getAtomList();
                    for (unsigned int i = 0 ; i < atomListTmp.size(); i++) {
                        atomList.push_back(atomListTmp[i]);
                    }

/*
std::cout << "Metal Atom " << metalAtom->getParent()->getName() << " "
          << metalAtom->getParent()->getSubMolId()
          << metalAtom->getParent()->getNumAtoms() << " "
          << std::endl;
*/
                }
              }

//std::cout << "NEW MOLECULE \n";

              // Copy bonds
              if (pNewMolecule) {
                for (unsigned int m = 0; m < this->itsMoleculeList.size(); m++) {
                  if (this->itsMoleculeList[m]->getName() == "Reference") {
                    continue;
                  }
                  std::map<int, Bond*> bondMap = this->itsMoleculeList[m]->getBondMap();

                  if (!bondMap.empty()) {
                    for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
                      Bond* pBond = b->second;
                      atom* pBdAtom1 = pBond->atom1;
                      atom* pBdAtom2 = pBond->atom2;

                      atom* pNewAtom1 = pNewMolecule->getAtom(pBdAtom1->getFileID(), 0, 1, 0);
                      atom* pNewAtom2 = pNewMolecule->getAtom(pBdAtom2->getFileID(), 0, 1, 0);
                      if (pNewAtom1 and pNewAtom2) {
                        pNewMolecule->addBond(pNewAtom1, pNewAtom2, pBond->type,
                                      pBond->stereo, pBond->topology, pBond->size);

/*
std::cout << "\n\n" << pBdAtom1->getParent()->getParent()->getName() << " "
          << pBdAtom1->getParent()->getName() << ":|" << pBdAtom1->getName() << " " << pBdAtom1->getFileID() << "|-" 
          << pBdAtom2->getParent()->getName() << ":|" << pBdAtom2->getName() << " " << pBdAtom2->getFileID() << "|" << std::endl;

std::cout << pNewAtom1->getParent()->getName() << ":|" << pNewAtom1->getName() << "|-" 
          << pNewAtom2->getParent()->getName() << ":|" << pNewAtom2->getName() << "|" << std::endl;
*/
                      }
                    }
                  }
                }

                for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
                  std::map<int, Bond*> bondMap = this->itsMetalCenters[m]->getBondMap();
                  if (!bondMap.empty()) {
                    for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
                      Bond* pBond = b->second;
                      atom* pBdAtom1 = pBond->atom1;
                      atom* pBdAtom2 = pBond->atom2;

                      atom* pNewAtom1 = pNewMolecule->getAtom(pBdAtom1->getFileID(), 0, 1, 0);
                      atom* pNewAtom2 = pNewMolecule->getAtom(pBdAtom2->getFileID(), 0, 1, 0);

//std::cout << "asdf |" << pBdAtom1->getName() << "|-|" << pBdAtom2->getName() << "| " 
//          << pBdAtom1->getFileID() << " " << pBdAtom2->getFileID() << " "
//          << pNewAtom1 << " " << pNewAtom2
//          << std::endl;

                      if (pNewAtom1 and pNewAtom2) {
                        pNewMolecule->addBond(pNewAtom1, pNewAtom2, pBond->type,
                                      pBond->stereo, pBond->topology, pBond->size);
/*
std::cout << pNewAtom1->getParent()->getName() << ":|" << pNewAtom1->getName() << "|-" 
          << pNewAtom2->getParent()->getName() << ":|" << pNewAtom2->getName() << "|" << std::endl;
std::cout << "|" << pNewAtom1->getName() << "|-|" << pNewAtom2->getName() << "|" << std::endl;
*/
                      }
                    }
                  }
                }
              }

//std::cout << " After copy " << std::endl;

              int nAtomsA = pStdMolecule->getNumAtoms();
              int nAtomsB = pNewMolecule->getNumAtoms();

              int nBondsA = pStdMolecule->getNumBonds();
              int nBondsB = pNewMolecule->getNumBonds();

              int nHeavyAtomsA = pStdMolecule->getNumHeavyAtoms();
              int nHeavyAtomsB = pNewMolecule->getNumHeavyAtoms();

//std::cout << "CHECK " << nAtomsA << " " << nAtomsB << " " << nBondsA << " " << nBondsB << "\n";
              int w = 0;
              if (nAtomsA != nAtomsB) w++;
              int h = 0;
              if (nHeavyAtomsA != nHeavyAtomsB) h++;

//std::cout << " w = " << w << " h = " << h << std::endl;

              if (w or h) {
                std::string eMes = " " + stdGrpList[i]->getName() + " is not a match \n";
                eMes +=   "   std mol nAtoms: " + i2s(nAtomsA) +
                        "\n   new mol nAtoms: " + i2s(nAtomsB) +
                        "\n   std mol nHeavyAtoms: " + i2s(nHeavyAtomsA) +
                        "\n   new mol nHeavyAtoms: " + i2s(nHeavyAtomsB) +
                        "\n   std mol nBonds: " + i2s(pStdMolecule->numBonds()) +
                        "\n   new mol nBonds: " + i2s(pNewMolecule->numBonds()) + "\n";
                MTKpp::errorLogger.throwError("collection::assignMetalParameters", eMes, INFO);
              }
              else {
                if (w == 0) {
                  // Determine corresponences between two molecules
                  std::vector<std::vector<int> > correspondenceMatrices;
                  //int cor = 0;
//std::cout << " superimpose " << std::endl;
                  superimpose* pSuperImpose = new superimpose();
                  int f = 0;
                  try {
                    f = pSuperImpose->initializeCorrespondences(pStdMolecule, pNewMolecule, 3, correspondenceMatrices);
                  }
                  catch (...) {
                    f = 1;
                  }

//std::cout << " after superimpose " << std::endl;

                  if (f != 0) {
                    std::string eMes = " Error initializing correspondences\n " + stdGrpList[i]->getName() + " is not a match \n";
                    eMes +=   "   std mol nAtoms: " + i2s(nAtomsA) +
                            "\n   new mol nAtoms: " + i2s(nAtomsB) +
                            "\n   std mol nHeavyAtoms: " + i2s(nHeavyAtomsA) +
                            "\n   new mol nHeavyAtoms: " + i2s(nHeavyAtomsB) +
                            "\n   std mol nBonds: " + i2s(pStdMolecule->numBonds()) +
                            "\n   new mol nBonds: " + i2s(pNewMolecule->numBonds()) + "\n";
                    MTKpp::errorLogger.throwError("collection::assignMetalParameters", eMes, INFO);
                  }
                  else if (f == 0) {
                    std::vector< vector3d > CoordsB;
                    pNewMolecule->getCoordinates(CoordsB);
                    int cor = 0;

//std::cout << " f = 0 " << std::endl;

                    double dRMSD = pSuperImpose->rmsd(pStdMolecule, CoordsB, correspondenceMatrices, cor);
                    std::string eMes = " " + stdGrpList[i]->getName() + " is a match with rmsd of " + d2s(dRMSD) + "\n";
                    MTKpp::errorLogger.throwError("collection::assignMetalParameters", eMes, INFO);
                    delete pSuperImpose;

                    if (h == 0) {
                      int nHeavyAtoms = pNewMolecule->getNumHeavyAtoms();

                      std::vector<atom*> molHeavyAtoms1 = pStdMolecule->getHeavyAtomList();
                      std::vector<atom*> molHeavyAtoms2 = pNewMolecule->getHeavyAtomList();

                      for (int t = 0; t < nHeavyAtoms; t++) {
                        for (int t2 = 0; t2 < nHeavyAtoms; t2++) {
                          if (correspondenceMatrices[cor][t*nHeavyAtoms+t2]) {
                            atom* pA1 = molHeavyAtoms1[t];
                            atom* pA2 = molHeavyAtoms2[t2];

                            if (pA1 and pA2) {
                              // getAtom(int number, bool atomIndex, bool fileId, bool atomColIndex)
                              atom* pOrgAtom = this->getAtom(pA2->getFileID(), 0, 1, 0);
                              // atom* pOrgAtom = this->getAtom(pA2->getColIndex(), 0, 0, 1);


//std::cout << pA1->getName() << "|" << pA2->getName() << "|" << pOrgAtom->getName() << " "
//          << pA1->getParent()->getName() << "|" << pA2->getParent()->getName() << "|" << pOrgAtom->getParent()->getName() << " |"
//          << pA1->getFileID() << "|" << pA2->getFileID() << "|" << pOrgAtom->getFileID() << "|"
//          << pA1->getColIndex() << "|" << pA2->getColIndex() << "|" << pOrgAtom->getColIndex() 
//                        << " |\n";

                              //pOrgAtom->getParent()->setName(pA1->getParent()->getName());
                              
                              for (unsigned int i = 0 ; i < atomList.size(); i++) {
                                  if (atomList[i]->getFileID() == pA2->getFileID()) {
                                    //std::cout <<  "setting " << atomList[i]->getParent()->getName() << " to " 
                                    //          << pA1->getParent()->getName() << "\n";
                                    atomList[i]->getParent()->setName(pA1->getParent()->getName());
                                  }
                              }
//std::cout << pOrgAtom->getParent() << " "  << pOrgAtom->getParent()->getName() << " " 
//<< pA2->getFileID()
//          << " "  << pOrgAtom->getParent()->getColIndex() << "\n";
                            }
                          }
                        }
                      }

                    }
                  }
                }
              }

              // delete new molecule
              if (pNewMolecule) {
                this->delMolecule(pNewMolecule);
              }
            }
          }
        }
      }
    }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
    typedef std::map<stdAtom*, atom*>::iterator funcGroupAtomIterator;

    for (unsigned int m = 0; m < this->itsMetalCenters.size(); m++) {
      std::string errMessage = "";
      metalCenter* metCen = this->itsMetalCenters[m];
      vector3d* metalCoords = metCen->getMetalAtom()->getCoords();

      std::string cMetal = metCen->getMetalAtom()->getElementSymbol();
      std::string stdGroupName = cMetal + "-";
      std::string priShell = metCen->getPrimaryShell();
      std::sort(priShell.begin(), priShell.end());
      stdGroupName += priShell;

      std::string secShell = metCen->getSecondaryShell();
      std::sort(secShell.begin(), secShell.end());

      if (this->pStdLibrary) {
        std::vector<stdGroup*> stdGrpList = this->pStdLibrary->getStdGroupList();
        stdGroup* pStdGroup = 0;

        std::vector<atom*> primShellAtoms;
        std::vector<submolecule*> primShellResidues;
        std::vector<molecule*> primShellMolecules;
        metCen->getPrimaryShellAtoms(primShellAtoms);
        for (unsigned int i = 0; i < primShellAtoms.size(); i++) {
          primShellResidues.push_back(primShellAtoms[i]->getParent());
          primShellMolecules.push_back(primShellResidues[i]->getParent());
        }

        int nBondedAtoms = static_cast<int>(primShellAtoms.size());
        std::vector<int> stdResiduesIndex;
        atom* pAtomX = 0;
        stdAtom* pStdAtomX = 0;
        funcGroup* pAtomXFuncGrp = 0;
        std::string pAtomXFuncGrpName = "";

        // If primary shell has an X (... currently only works for one X residue ...)
        std::string::size_type loc = stdGroupName.find("X", 0);
        if (loc != std::string::npos) {

          for (unsigned int i = 0; i < primShellAtoms.size(); i++) {
            if (pStdLibrary->getL(primShellResidues[i]->getName()) == "X") {

              pAtomX = primShellAtoms[i];
              if (!pAtomX) {
                exit(1);
              }

              pAtomXFuncGrp = primShellMolecules[i]->getFunctionalGroup(pAtomX);
              if (!pAtomXFuncGrp) {
                std::cout << " Error in collection::assignMetalParameters |"
                          << pAtomX->getName() << "|@" << pAtomX->getParent()->getName()
                          << " ... exiting " << std::endl;
                exit(1);
              }

              pAtomXFuncGrpName = pAtomXFuncGrp->pStdFrag->getSymbol();
              if (pAtomXFuncGrpName != "") {
                exit(1);
              }
    
              //std::cout << "    " << pAtomX->getName() << ":" << pAtomX->getFileID()
              //          << ":" << pAtomXFuncGrpName << ":" << pAtomXFuncGrp->pStdFrag->getName() << std::endl;

              stdResiduesIndex.push_back(0);

              if (!pAtomXFuncGrp->atomMap.empty()) {
                for (funcGroupAtomIterator b = pAtomXFuncGrp->atomMap.begin();
                     b != pAtomXFuncGrp->atomMap.end(); b++) {
                  std::cout << b->first->identity << " : " << b->second->getName() << std::endl;
                  if (pAtomX == b->second) {
                    pStdAtomX = b->first;
                  }
                }
              }
              if (!pStdAtomX) {
                std::cout << " collection::assignMetalParameters:Error ... exiting " << std::endl;
                //exit(0);
                throw MTKException(" collection::determineMetalEnvironments:Error ... exiting ");
              }
            }
            else {
              stdResiduesIndex.push_back(1);
            }
          }

          for (unsigned int g = 0; g < stdGrpList.size(); g++) {
            stdGroup* pStdGrp = stdGrpList[g];
            //std::cout << pStdGrp->getName() << "  " <<stdGroupName << std::endl;
            if (pStdGrp->getName() == stdGroupName) {
              std::string stdGrpInfo = pStdGrp->getInfo();
              std::string::size_type loc2 = stdGrpInfo.find(pAtomXFuncGrpName, 0);
              if (loc2 != std::string::npos) {
                //std::cout << " MATCH " << std::endl;
                pStdGroup = stdGrpList[g];
              }
            }
          }
        }
        else {
          stdResiduesIndex.assign(nBondedAtoms, 1);
          pStdGroup = this->pStdLibrary->getStdGroup(stdGroupName);
        }

        if (pStdGroup) {
          metCen->setStdGroup(pStdGroup);
          std::vector<stdFrag*> stdFrags = pStdGroup->getStdFragList();
          stdFrag* pStdMetalFrag = stdFrags[stdFrags.size()-1];
          submolecule* pMetalSubmolecule = this->itsMetalAtoms[m]->getParent();
          pMetalSubmolecule->setStdFrag(pStdMetalFrag);
          pMetalSubmolecule->setName(pStdMetalFrag->getSymbol());

          std::vector<stdAtom*> metalAtoms =  pStdMetalFrag->getStdAtomList();
          if (metalAtoms.size() > 1) {
            errMessage = " Currently only one metal ion is supported";
            MTKpp::errorLogger.throwError("collection::assignMetalParameters", errMessage, ERROR);
            exit(0);
          }
          stdAtom* pStdMetal = metalAtoms[0];
          this->itsMetalAtoms[m]->setStdAtom(pStdMetal);

          // Determine which stdGroup correspond to which attached residue
          int* matchMatrix;
          try {
            matchMatrix   = new int [nBondedAtoms*nBondedAtoms];
          }
          catch (std::bad_alloc) {
            std::cout << " Memory Allocation Failure " << std::endl;
            //exit(0);
            throw MTKException(" Memory Allocation Failure ");
          }

          for (int y = 0; y < nBondedAtoms; y++) {
            for (int x = 0; x < nBondedAtoms; x++) {
              if (matchMatrix[y*nBondedAtoms+x]) {
                matchMatrix[y*nBondedAtoms+x] = 0;
              }
            }
          }
#ifdef DEBUG

          for (unsigned int i = 0; i < stdFrags.size()-1; i++) {
            for (unsigned int j = 0; j < primShellResidues.size(); j++) {
              std::cout << " " << stdFrags[i]->getSymbol() << "/" << primShellResidues[j]->getName();
            }
            std::cout << " \n";
          }
#endif

// i loop stdFrags
// j loop primShellResidues

          for (unsigned int i = 0; i < stdFrags.size()-1; i++) {
            int nStdFragHeavyAtoms = stdFrags[i]->numStdHeavyAtoms();
            for (unsigned int j = 0; j < primShellResidues.size(); j++) {
              if (stdResiduesIndex[j]) {
                int nHeavyAtoms = primShellResidues[j]->numHeavyAtoms();
                std::string atName = primShellAtoms[j]->getName();
                if (nStdFragHeavyAtoms == nHeavyAtoms) {
                  if (stdFrags[i]->hasStdAtom(atName)) {
                    if (primShellResidues[j]->getName() == "HIS" or
                        primShellResidues[j]->getName() == "HIE" or
                        primShellResidues[j]->getName() == "HID") {
                      if (stdFrags[i]->hasStdAtom(" HE2") and atName == " ND1") { // HIE
                        matchMatrix[i*nBondedAtoms+j] = 1;
                      }
                      else if (stdFrags[i]->hasStdAtom(" HD1") and atName == " NE2") { // HID
                        matchMatrix[i*nBondedAtoms+j] = 1;
                      }
                      else {
//
//                        std::cout << " NO MATCH: " << atName << " " << primShellResidues[j]->getName()
//                        << " " << stdFrags[i]->getSymbol() << std::endl;
//
                      }
                    }
                    else {
                      matchMatrix[i*nBondedAtoms+j] = 1;
                    }
                  }
                }
              }
              else {
                int nHeavyAtoms = pAtomXFuncGrp->pStdFrag->numStdHeavyAtoms();
                if (nStdFragHeavyAtoms == nHeavyAtoms) {
                  matchMatrix[i*nBondedAtoms+j] = 1;
                }
              }
            }
          }

          // Check for mismatches
          for (int t = 0; t < nBondedAtoms; t++) {
            int check = 0;
            for (int tt = 0; tt < nBondedAtoms; tt++) {
              check += matchMatrix[t*nBondedAtoms+tt];
            }
            if (check == 0) {
//
//              std::cout << " collection::determineMetalEnvironments:";
//              std::cout << "Error in match matrix ... exiting " << std::endl;
//
            }
          }

#ifdef DEBUG
          std::cout << " collection::assignMetalParameters: match matrix " << std::endl;
          for (int t = 0; t < nBondedAtoms; t++) {
            std::cout <<  primShellAtoms[t]->getName() << " ";
          }
          std::cout << " " << std::endl;
          for (int t = 0; t < nBondedAtoms; t++) {
            for (int tt = 0; tt < nBondedAtoms; tt++) {
              std::cout << matchMatrix[t*nBondedAtoms+tt] << " ";
            }
            std::cout << " " << std::endl;
          }
          std::cout << " " << std::endl;
#endif

          std::vector<int> match;
          std::vector<std::vector<int> > matches;
          this->findMatchings(0, nBondedAtoms, matchMatrix, match, matches);

#ifdef DEBUG
          std::cout << " collection::assignMetalParameters: match matrix " << std::endl;
          for (int t = 0; t < nBondedAtoms; t++) {
            std::cout <<  primShellAtoms[t]->getName() << " ";
          }
          for (unsigned int t = 0; t < matches.size(); t++) {
            std::cout << " \n Match " << t << std::endl;
            int indexT = 0;
            for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
              for (int tttt = 0; tttt < nBondedAtoms; tttt++) {
                std::cout << matches[t][indexT] << " ";
                indexT++;
              }
              std::cout << " " << std::endl;
            }
          }
          std::cout << " " << std::endl;
#endif

          vector3d*      metalCoord = this->itsMetalAtoms[m]->getCoords();
          vector3d*      coord1;
          vector3d*      coord2;
          std::vector<double> bondDistancesPdb;
          for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
            coord1 = primShellAtoms[ttt]->getCoords();
            double d = coord1->dist(*metalCoord);
            bondDistancesPdb.push_back(d);
          }

          std::vector<double> angleSizesPdb;
          for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
            coord1 = primShellAtoms[ttt]->getCoords();
            for (int tttt = 0; tttt < nBondedAtoms; tttt++) {
              if (ttt == tttt) continue;
              coord2 = primShellAtoms[tttt]->getCoords();
              double a = RAD2DEG * angle(*(coord1), *(metalCoord), *(coord2));
              angleSizesPdb.push_back(a);
            }
          }
//
//std::cout << "\n primShell: ";
//for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
//  std::cout << primShellResidues[ttt]->getName()  << " ";
//}

//std::cout << "\n stdFrags: ";
//for (int tttt = 0; tttt < nBondedAtoms; tttt++) {
//  std::cout << stdFrags[tttt]->getSymbol() << " ";
//}
//std::cout << "\n\n ttt/tttt ps/stdFrag\n";
//

          double bestScore = BIGNUM;
          int bestMatch = 0;
          for (unsigned int t = 0; t < matches.size(); t++) {
            double matchScore = 0.0;
            std::vector<stdFrag*> lStdFrags;
            std::vector<stdAtom*> lStdAtoms;
            int nJKasdf = 0;
            int indexT = 0;

// ttt loop stdFrags
// tttt loop primShellResidues

            for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
              for (int tttt = 0; tttt < nBondedAtoms; tttt++) {
                if (matches[t][indexT]) {
                  lStdFrags.push_back(stdFrags[ttt]);

//std::cout << " " << ttt << "/" << tttt << " " 
//<< primShellResidues[tttt]->getName() << "/" << stdFrags[ttt]->getSymbol() << std::endl;

                  if (stdFrags[ttt]->getStdAtom(primShellAtoms[tttt]->getName())) {
                    stdAtom* pStdA = stdFrags[ttt]->getStdAtom(primShellAtoms[tttt]->getName());
                    lStdAtoms.push_back(pStdA);
                  }
                  else {
//std::cout << "      --------------- \n";
                    lStdAtoms.push_back(pStdAtomX);
                    nJKasdf++;
                  }
                }
                indexT++;
              }
            }
//std::cout << " \n";

//
//            for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
//std::cout << lStdFrags[ttt]->getName() << " " << primShellAtoms[ttt]->getName() <<std::endl;
//              if (lStdFrags[ttt]->getStdAtom(primShellAtoms[ttt]->getName())) {
//                stdAtom* pStdA = lStdFrags[ttt]->getStdAtom(primShellAtoms[ttt]->getName());
//                lStdAtoms.push_back(pStdA);
//              }
//              else {
//std::cout << "      --------------- \n";
//                lStdAtoms.push_back(pStdAtomX);
//                nJKasdf++;
//              }
//            }
//
            if (nJKasdf > 1) {
              std::cout << " collection::assignMetalParameters:Error ... EXITING " << std::endl;
              exit(0);
            }
//
//            for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
//                std::cout << lStdAtoms[ttt]->identity << " ";
//            }
//            std::cout << " " << std::endl;
//

            double eBond = 0.0;
            for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
              bondParam* pBondP = pParameters->getBondParam(pStdMetal->type, lStdAtoms[ttt]->type);
              eBond += pow(bondDistancesPdb[ttt] - pBondP->req, 2);
            }

            double eAngle = 0.0;
            int angleIndex = 0;
            for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
              for (int tttt = 0; tttt < nBondedAtoms; tttt++) {
                if (ttt == tttt) continue;
                angleParam* pAngleP = pParameters->getAngleParam(lStdAtoms[ttt]->type,
                            pStdMetal->type, lStdAtoms[tttt]->type);
                if (pAngleP) {
                  eAngle += pow((angleSizesPdb[angleIndex] - pAngleP->req * RAD2DEG), 2);
                  angleIndex++;
                }
                else {
                  std::cout << " collection::assignMetalParameters:Error ... exiting " << std::endl;
                  exit(0);
                }
              }
            }

            //std::cout << " eBOND = " << eBond << " eAngle = " << eAngle << std::endl;
            matchScore = 100.0*eBond + 0.1*eAngle;
            if (matchScore < bestScore) {
              bestScore = matchScore;
              bestMatch = t;
            }
          }

//
//          std::cout << " BEST SCORE = " << bestScore << std::endl;
//          std::cout << " BEST MATCH = " << bestMatch << std::endl;
//

// i loop stdFrags
// j loop primShellResidues

          int indexT = 0;
          for (int ttt = 0; ttt < nBondedAtoms; ttt++) {
            for (int tttt = 0; tttt < nBondedAtoms; tttt++) {
              if (matches[bestMatch][indexT]) {
                if (stdResiduesIndex[ttt]) {

                  // Set standard residue and name
                  primShellResidues[tttt]->setStdFrag(stdFrags[ttt]);
                  primShellResidues[tttt]->setName(stdFrags[ttt]->getSymbol());
//
//                  std::cout << primShellResidues[tttt]->getName() << " "
//                            << primShellResidues[tttt]->getSubMolId() << " "
//                            << stdFrags[ttt]->getSymbol() << std::endl;
//
                  std::vector<atom*> resAtomList = primShellResidues[tttt]->getAtomList();
                  for (unsigned int u = 0; u < resAtomList.size(); u++) {
                    stdAtom* plocalStdAtom = stdFrags[ttt]->getStdAtom(resAtomList[u]->getName());
                    //if (stdFrags[ttt]->getSymbol() == "CY1") std::cout <<  << " " << <<"\n";
                    resAtomList[u]->setName(plocalStdAtom->identity); // ensures name in pdb file is the same as the prep file
                    if (plocalStdAtom) {
//
//                      stdAtom* pOldStdAtom = resAtomList[u]->getStdAtom();
//                      std::cout << " ATOM NAME:: " << resAtomList[u]->getName() << std::endl;
//                      std::cout << " ATOM TYPE:: " << pOldStdAtom->type << ":" << plocalStdAtom->type << std::endl;
//
                      resAtomList[u]->setStdAtom(plocalStdAtom);
                    }
                    else {
                      //std::cout << " Can't find stdatom (1) for " << resAtomList[u]->getParent()->getName() << " " 
                      //          << resAtomList[u]->getName() << std::endl;
                      //std::cout << "   collection::assignMetalParameters:Error ... exiting " << std::endl;
                      //exit(0);

                      std::stringstream ss;
                      ss << " Can't find stdatom (1) for " << resAtomList[u]->getName() << std::endl;
                      ss << "   collection::determineMetalEnvironments:Error ... exiting " << std::endl;
                      std::cout << ss.str();
                      throw MTKException(ss.str());
                    }
                  }
                }
                else {

                  std::map<atom*, int> atomsAssigned;
                  std::vector<atom*> atomXsmolAtoms = primShellResidues[tttt]->getAtomList();
                  for (unsigned int f = 0; f < atomXsmolAtoms.size(); f++) {
                    atomsAssigned[atomXsmolAtoms[f]] = 0;
                  }

                  if (!pAtomXFuncGrp->atomMap.empty()) {
                    for (funcGroupAtomIterator b = pAtomXFuncGrp->atomMap.begin();
                         b != pAtomXFuncGrp->atomMap.end(); b++) {
                      stdAtom* plocalStdAtom = stdFrags[ttt]->getStdAtom(b->first->identity);
                      if (plocalStdAtom) {
                        stdAtom* pOldStdAtom = b->second->getStdAtom();
                        plocalStdAtom->type = pOldStdAtom->type;
                        //std::cout << " ATOM TYPE:: " << pOldStdAtom->type << ":" << plocalStdAtom->type << std::endl;
                        b->second->setStdAtom(plocalStdAtom);
                        atomsAssigned[b->second] = 1;
                      }
                      else {
                        std::cout << " Can't find stdatom 2v"<< std::endl;
                        std::cout << "   collection::assignMetalParameters:Error ... exiting " << std::endl;
                        //exit(0);
                        throw MTKException(" collection::assignMetalParameters:Error ... EXITING ");
                      }
                    }
                  }

                  unsigned int nStdFuncGrpAtoms = pAtomXFuncGrp->atomMap.size();
                  unsigned int nMolAtoms = primShellResidues[tttt]->getStdFrag()->numStdAtoms();
                  int nOtherStdFuncGrpAtoms = nMolAtoms - nStdFuncGrpAtoms;
                  std::vector<funcGroup*> allFuncGroups = primShellMolecules[tttt]->getFunctionalGroups();
                  for (int unsigned p = 0; p < allFuncGroups.size(); p++) {
                    if (allFuncGroups[p] != pAtomXFuncGrp) {
                      if (allFuncGroups[p]->pStdFrag->numStdAtoms() == nOtherStdFuncGrpAtoms) {
                        bool bOK = true;
                        for (funcGroupAtomIterator b = allFuncGroups[p]->atomMap.begin();
                             b != allFuncGroups[p]->atomMap.end(); b++) {
                          if (atomsAssigned[b->second]) {
                            bOK = false;
                          }
                        }
                        if (bOK) {
                          for (funcGroupAtomIterator b = allFuncGroups[p]->atomMap.begin();
                               b != allFuncGroups[p]->atomMap.end(); b++) {
//
//                            std::cout << " ATOM TYPE:: " << b->second->getStdAtom()->type
//                                      << ":" << b->first->type << std::endl;
//
                            b->first->type = b->second->getStdAtom()->type;
                            b->second->setStdAtom(b->first);
                          }
                          break;
                        }
                      }
                    }
                  }
                }
              }
              indexT++;
            }
          }
///////////////////////////////////////////////////////////////////////////////
          std::vector<atom*> torAtoms;
          errMessage += "\n";
          // Add bonds and angles to metal center
          for (int e = 0; e < nBondedAtoms; e++) {
            torAtoms.push_back(primShellAtoms[e]);

            // Bonds
            Bond* localBond = metCen->addBond(this->itsMetalAtoms[m], primShellAtoms[e]);
            this->itsMetalAtoms[m]->addBondedAtom(primShellAtoms[e]);
            primShellAtoms[e]->addBondedAtom(this->itsMetalAtoms[m]);
            if (!localBond) {
              std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
              //exit(1);
              throw MTKException(" collection::assignMetalParameters:Error ... exiting ");
            }
            errMessage += "    BOND: "
                        + this->itsMetalAtoms[m]->getParent()->getName() + "@|" +
                                this->itsMetalAtoms[m]->getName() + "|-"
                        + primShellAtoms[e]->getParent()->getName() + "@|" +
                          primShellAtoms[e]->getName() + "|             dist = " + d2s(localBond->size) +"\n";

            // Angles
            for (int d = 0; d < nBondedAtoms; d++) {
              if (e == d) continue;
              Angle* localAngle = metCen->addAngle(primShellAtoms[e], this->itsMetalAtoms[m], primShellAtoms[d]);
              if (!localAngle) {
                std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
              }
              bool added = primShellAtoms[e]->addBonded13Atom(primShellAtoms[d]);
              bool added2 = primShellAtoms[d]->addBonded13Atom(primShellAtoms[e]);
              if (added or added2) {
                errMessage += "   ANGLE: "
                            + primShellAtoms[e]->getParent()->getName() + "@|" + primShellAtoms[e]->getName() + "|-"
                            + this->itsMetalAtoms[m]->getParent()->getName() + "@|" + this->itsMetalAtoms[m]->getName() + "|-"
                            + primShellAtoms[d]->getParent()->getName() + "@|" + primShellAtoms[d]->getName() 
                            + "| angle = " + d2s(localAngle->size * RAD2DEG) + "\n";
              }
              //metCen->add13Bond(primShellAtoms[e], primShellAtoms[d]);
            }
            std::vector<atom*> lBondedAtoms = primShellAtoms[e]->getBondedAtoms();
            for (unsigned int d = 0; d < lBondedAtoms.size(); d++) {
              if (this->itsMetalAtoms[m] == lBondedAtoms[d]) continue;
              Angle* localAngle = metCen->addAngle(this->itsMetalAtoms[m], primShellAtoms[e], lBondedAtoms[d]);
              if (!localAngle) {
                std::cout << "  collection::assignMetalParameters:Error ... exiting " << std::endl;
              }
              bool added = this->itsMetalAtoms[m]->addBonded13Atom(lBondedAtoms[d]);
              bool added2 = lBondedAtoms[d]->addBonded13Atom(this->itsMetalAtoms[m]);
              if (added or added2) {
                errMessage += "   ANGLE: "
                            + this->itsMetalAtoms[m]->getParent()->getName() + "@|" + primShellAtoms[e]->getName() + "|-"
                            + primShellAtoms[e]->getParent()->getName() + "@|" + primShellAtoms[e]->getName() + "|-"
                            + lBondedAtoms[d]->getParent()->getName() + "@|" + lBondedAtoms[d]->getName() 
                            + "| angle = " + d2s(localAngle->size * RAD2DEG) + "\n";
              }
            }
          }
          torAtoms.push_back(this->itsMetalAtoms[m]);

          // Add torsions to metal center
          for (unsigned int e = 0; e < torAtoms.size(); e++) {
            atom* pAtom1 = torAtoms[e];
            for (unsigned int l = 0; l < pAtom1->bondedAtoms.size(); l++) {
              atom* pAtom2 = pAtom1->bondedAtoms[l];
              for (unsigned int o = 0; o < pAtom2->bondedAtoms.size(); o++) {
                atom* pAtom3 = pAtom2->bondedAtoms[o];
                if (pAtom1 != pAtom3) {
                  for (unsigned int p = 0; p < pAtom3->bondedAtoms.size(); p++) {
                    atom* pAtom4 = pAtom3->bondedAtoms[p];
                    if (pAtom4 == pAtom2) continue;
                    bool gotTorsion = metCen->hasTorsion(pAtom1, pAtom2, pAtom3, pAtom4);
                    if (!gotTorsion) {
                      double dTorsion = torsion(*(pAtom1->getCoords()), *(pAtom2->getCoords()),
                                         *(pAtom3->getCoords()), *(pAtom4->getCoords()));
                      metCen->addTorsion(pAtom1, pAtom2, pAtom3, pAtom4, dTorsion);
                      pAtom1->addBonded14Atom(pAtom4);
                      pAtom4->addBonded14Atom(pAtom1);
//
//                      std::cout << "   TORSION: "
//                          << pAtom1->getParent()->getName() << "@|" << pAtom1->getName() << "|-"
//                          << pAtom2->getParent()->getName() << "@|" << pAtom2->getName() << "|-"
//                          << pAtom3->getParent()->getName() << "@|" << pAtom3->getName() << "|-"
//                          << pAtom4->getParent()->getName() << "@|" << pAtom4->getName() << "|"
//                          << std::endl;
//
                    }
                  }
                }
              }
            }
          }
//
//          std::string info = "";
//          info = metCen->getInfo();
//          MTKpp::errorLogger.throwError("collection::determineMetalEnvironments", info, INFO);
//
          delete matchMatrix;
        }
        else {
          std::vector<atom*> primShellAtoms;
          metCen->getPrimaryShellAtoms(primShellAtoms);

          errMessage += "\n  stdGroup: " + stdGroupName + " not found. Bonded Atoms = ";
          for (unsigned int x = 0; x < primShellAtoms.size(); x++) {
            errMessage += primShellAtoms[x]->getElementSymbol();
          }
          //errMessage += " \n";

          for (unsigned int x = 0; x < primShellAtoms.size(); x++) {
            errMessage += " \n     M-" + primShellAtoms[x]->getElementSymbol() + " Bond Distance = "
                       + d2s(metalCoords->dist( (*primShellAtoms[x]->getCoords()) ));
          }
        }
      }
      MTKpp::errorLogger.throwError("collection::assignMetalParameters", errMessage, INFO);
    }
*/

/*
    double localTotalCharge = 0.0;
    for (unsigned int m = 0; m < this->itsMoleculeList.size(); m++) {

      std::vector<submolecule*> smolList = this->itsMoleculeList[m]->getSubMoleculeList();
      for (unsigned int s = 0; s < smolList.size(); s++) {
        std::vector<atom*> atomListA = smolList[s]->getAtomList();
        double localMolCharge = 0;
        for (unsigned int a = 0; a < atomListA.size(); a++) {
          localMolCharge += atomListA[a]->getStdAtom()->atmCharge;
        }
        std::cout << "  " << smolList[s]->getName() << localMolCharge;
      }

      std::vector<atom*> atomList = this->itsMoleculeList[m]->getAtomList();
      double localMolCharge = 0;
      for (unsigned int a = 0; a < atomList.size(); a++) {
        if (!atomList[a]->getStdAtom()) {

          //std::cout << " CAN'T FIND STDATOM FOR " << atomList[a]->getParent()->getName()
          //          << " |" << atomList[a]->getName() << "|" << std::endl;
          //exit(0);

          std::stringstream ss;
          ss << " CAN'T FIND STDATOM FOR " << atomList[a]->getParent()->getName()
                    << " |" << atomList[a]->getName() << "|" << std::endl;
          std::cout << ss.str();
          throw MTKException(ss.str());

        }
        std::cout << "   " << atomList[a]->getParent()->getName() << " |" << atomList[a]->getName() << "| "
                  << atomList[a]->getStdAtom()->atmCharge << std::endl;
        localMolCharge += atomList[a]->getStdAtom()->atmCharge;
      }
      localTotalCharge += localMolCharge;
      std::cout << "      Charge = " << localMolCharge << std::endl;
    }
    std::cout << "       TOTAL CHARGE = " << localTotalCharge << std::endl;
*/

}

// ============================================================
// Function : getMetalCenters()
// ------------------------------------------------------------
//
// ============================================================
std::vector<metalCenter*> collection::getMetalCenters()
{
    return this->itsMetalCenters;
}

// ============================================================
// Function : findMatchings()
// ------------------------------------------------------------
//
// ============================================================
void collection::findMatchings(int pos, int nBondedAtoms,
                 int matchMatrix[], std::vector<int>& match,
                 std::vector<std::vector<int> >& matches)
{
    // back up match matrix
    int matchMatrixBackUp[nBondedAtoms*nBondedAtoms];
    for (int i = 0; i < nBondedAtoms; i++) {
      for (int j = 0; j < nBondedAtoms; j++) {
        matchMatrixBackUp[i*nBondedAtoms+j] = matchMatrix[i*nBondedAtoms+j];
      }
    }

    bool mismatch = false;
    for (int i = 0; i < nBondedAtoms; i++) {
      if (matchMatrix[pos*nBondedAtoms+i]) {
        this->updateMatchMatrix(i, pos, nBondedAtoms, matchMatrix);
        this->refineMatchMatrix(nBondedAtoms, matchMatrix, mismatch);

        if (!mismatch) {
          if (pos == nBondedAtoms-1) {
            for (int f = 0; f < nBondedAtoms; f++) {
              for (int g = 0; g < nBondedAtoms; g++) {
                match.push_back(matchMatrix[f*nBondedAtoms+g]);
                matchMatrix[f*nBondedAtoms+g] = matchMatrixBackUp[f*nBondedAtoms+g];
              }
            }
            matches.push_back(match);
            match.clear();
          }
          else {
            findMatchings(pos+1, nBondedAtoms, matchMatrix, match, matches);
          }
        }
        else {
          // reset match matrix
          for (int f = 0; f < nBondedAtoms; f++) {
            for (int j = 0; j < nBondedAtoms; j++) {
              matchMatrix[f*nBondedAtoms+j] = matchMatrixBackUp[f*nBondedAtoms+j];
            }
          }
        }
      }
      // reset match matrix
      for (int f = 0; f < nBondedAtoms; f++) {
        for (int j = 0; j < nBondedAtoms; j++) {
          matchMatrix[f*nBondedAtoms+j] = matchMatrixBackUp[f*nBondedAtoms+j];
        }
      }
    }
}

// ============================================================
// Function : updateMatchMatrix()
// ------------------------------------------------------------
//
// ============================================================
void collection::updateMatchMatrix(int iPos, int jPos, int nBondedAtoms, int matchMatrix[])
{
    // zero out the row
    for (int i = 0; i < nBondedAtoms; i++) {
      if (i == iPos) {
        continue;
      }
      else {
        matchMatrix[jPos*nBondedAtoms+i] = 0;
      }
    }

    // zero out the column
    for (int i = 0; i < nBondedAtoms; i++) {
      if (i == jPos) {
        continue;
      }
      else {
        matchMatrix[i*nBondedAtoms+iPos] = 0;
      }
    }
}

// ============================================================
// Function : refineMatchMatrix()
// ------------------------------------------------------------
//
// ============================================================
void collection::refineMatchMatrix(int nBondedAtoms, int matchMatrix[], bool &mismatch)
{
    for (int i = 0; i < nBondedAtoms; i++) {
      bool rowOk = false;
      for (int j = 0; j < nBondedAtoms; j++) {
        if (matchMatrix[i*nBondedAtoms+j] > 0) rowOk = true;
      }
      if (!rowOk) {
        mismatch = true;
      }
    }
}

} // MTKpp namespace
