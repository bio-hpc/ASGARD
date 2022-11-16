/*!
   \file parameters.cpp
   \brief Container for parameter information
   \author Martin Peters

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.14 $

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

#include "parameters.h"
#include "element.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

parameters* Parameters=NULL;

// ============================================================
// Function : getInstance()
// ------------------------------------------------------------
//
// ============================================================
parameters* parameters::getInstance(elements* p)
{
    if (Parameters == NULL) {
      Parameters=new parameters(p);
    }
    return Parameters;
}

// ============================================================
// Function : paramParser()
// ------------------------------------------------------------
// Constructor for the class
// ============================================================
parameters::parameters(elements* p)
{
    pElements = p;
    bLJ612SE = false;
}

// ============================================================
// Function : parameters()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
parameters::~parameters()
{
    for (atomTypeIterator a = itsTypeList.begin(); a != itsTypeList.end(); a++) {
      pAtomType = *a;
      delete pAtomType;
    }
    itsTypeList.clear();

    for (bondParamIterator b = itsBondList.begin(); b != itsBondList.end(); b++) {
      pBondParam = *b;
      delete pBondParam;
    }
    itsBondList.clear();

    for (angleParamIterator b = itsAngleList.begin(); b != itsAngleList.end(); b++) {
      pAngleParam = *b;
      delete pAngleParam;
    }
    itsAngleList.clear();

    for (torsionParamIterator b = itsTorsionList.begin(); b != itsTorsionList.end(); b++) {
      pTorsionParam = *b;
      delete pTorsionParam;
    }
    itsTorsionList.clear();

    for (improperParamIterator b = itsImproperList.begin(); b != itsImproperList.end(); b++) {
      pImproperParam = *b;
      delete pImproperParam;
    }
    itsImproperList.clear();

    for (hBondParamIterator b = itsHBondList.begin(); b != itsHBondList.end(); b++) {
      pHBondParam = *b;
      delete pHBondParam;
    }
    itsHBondList.clear();
}

// ============================================================
// Function : printAtomTypes()
// ------------------------------------------------------------
// 
// ============================================================
void parameters::printAtomTypes()
{
    for (atomTypeIterator a = itsTypeList.begin(); a != itsTypeList.end(); a++) {
      pAtomType = *a;
      std::cout << pAtomType->name << "\n"
                << "         element = " << pAtomType->element << "\n"
                << "   hybridization = " << pAtomType->hybridization << "\n"
                << "          rvalue = " << pAtomType->rvalue << "\n"
                << "          evalue = " << pAtomType->evalue << "\n"
                << "           group = " << pAtomType->groupName << "\n"
                << "     description = " << pAtomType->description << std::endl;
    }
}

// ============================================================
// Function : addAtomType()
// ------------------------------------------------------------
// 
// ============================================================
atomType* parameters::addAtomType()
{
    bLJ612SE = false;
    pAtomType = new atomType();
    pAtomType->name = "";
    pAtomType->element = "";
    pAtomType->atNum = 0;
    pAtomType->mass = 0.0;
    pAtomType->hybridization = "";
    pAtomType->description = "";
    pAtomType->rvalue = 0.0;
    pAtomType->evalue = 0.0;
    pAtomType->atomPolarizability = 0.0;
    pAtomType->groupName = "";
    pAtomType->optimize = false;
    this->itsTypeList.push_back(pAtomType);
    return pAtomType;
}

// ============================================================
// Function : addAtomType(atomType*)
// ------------------------------------------------------------
// 
// ============================================================
atomType* parameters::addAtomType(atomType* a, std::string n, std::string g)
{
    if (!a) {
      std::cout << " Error in parameters::addAtomType ... exiting " << std::endl;
      //exit(0);
      throw MTKException(" Error in parameters::addAtomType ... exiting ");
    }

    if (this->hasAtomType(n, g)) {
      //std::cout << " not adding: " << n << std::endl;
      return 0;
    }
    //std::cout << " adding: " << n << std::endl;

    pAtomType = this->addAtomType();
    if (!pAtomType) {
      std::cout << " Error in parameters::addAtomType ... exiting " << std::endl;
      //exit(0);
      throw MTKException(" Error in parameters::addAtomType ... exiting ");
    }

    pAtomType->name = n;
    pAtomType->element = a->element;
    pAtomType->mass = a->mass;
    pAtomType->hybridization = a->hybridization;
    pAtomType->description = a->description;
    pAtomType->rvalue = a->rvalue;
    pAtomType->evalue = a->evalue;
    pAtomType->atomPolarizability = a->atomPolarizability;
    pAtomType->groupName = g;
    pAtomType->optimize = a->optimize;

    // Copy bond parameters
    bondParam* pBdPrm = 0;
    bondParam* pBP = 0;
    std::vector<bondParam*> bondsToBeAdded;
    for (bondParamIterator c = this->itsBondList.begin(); c != this->itsBondList.end(); c++) {
      pBdPrm = *c;
      if ((pBdPrm->atomType1 == a->name) or
          (pBdPrm->atomType2 == a->name)) {
        //if (pBdPrm->groupName != g) {
        if (!this->hasBondParam(pBdPrm->atomType1, pBdPrm->atomType2, g)) {
          bondsToBeAdded.push_back(pBdPrm);
        }
      }
    }
    for (bondParamIterator c = bondsToBeAdded.begin(); c != bondsToBeAdded.end(); c++) {
      pBdPrm = *c;

      //std::cout << " Creating new bond type: " << pBdPrm->atomType1 << "-" << pBdPrm->atomType2 << " " << g
      //          << " from " <<  pBdPrm->groupName << std::endl;

      pBP = this->addBondParam();
      if (pBdPrm->atomType1 == a->name) {
        pBP->atomType1 = n;
      }
      else {
        pBP->atomType1 = pBdPrm->atomType1;
      }

      if (pBdPrm->atomType2 == a->name) {
        pBP->atomType2 = n;
      }
      else {
        pBP->atomType2 = pBdPrm->atomType2;
      }
      pBP->keq = pBdPrm->keq;
      pBP->req = pBdPrm->req;
      pBP->groupName = g;
      pBP->optimize = pBdPrm->optimize;
    }

    // Copy angle parameters
    angleParam* pAngPrm = 0;
    angleParam* pAP = 0;
    std::vector<angleParam*> anglesToBeAdded;
    for (angleParamIterator c = this->itsAngleList.begin(); c != this->itsAngleList.end(); c++) {
      pAngPrm = *c;
      if ((pAngPrm->atomType1 == a->name) or
          (pAngPrm->atomType2 == a->name) or
          (pAngPrm->atomType3 == a->name)) {
        if (!this->hasAngleParam(pAngPrm->atomType1, pAngPrm->atomType2, pAngPrm->atomType3, g)) {
          anglesToBeAdded.push_back(pAngPrm);
        }
      }
    }

    for (angleParamIterator c = anglesToBeAdded.begin(); c != anglesToBeAdded.end(); c++) {
      pAngPrm = *c;
      pAP = this->addAngleParam();
      if (pAngPrm->atomType1 == a->name) {
        pAP->atomType1 = n;
      }
      else {
        pAP->atomType1 = pAngPrm->atomType1;
      }

      if (pAngPrm->atomType2 == a->name) {
        pAP->atomType2 = n;
      }
      else {
        pAP->atomType2 = pAngPrm->atomType2;
      }

      if (pAngPrm->atomType3 == a->name) {
        pAP->atomType3 = n;
      }
      else {
        pAP->atomType3 = pAngPrm->atomType3;
      }

      pAP->keq = pAngPrm->keq;
      pAP->req = pAngPrm->req;
      pAP->groupName = g;
      pAP->optimize = pAngPrm->optimize;
    }

    // Copy Torsion parameter
    torsionParam* pTP = 0;
    torsionParam* pTorPrm = 0;
    std::vector<torsionParam*> torsionsToBeAdded;
    for (torsionParamIterator c = this->itsTorsionList.begin(); c != this->itsTorsionList.end(); c++) {
      pTorPrm = *c;
      if ((pTorPrm->atomType1 == a->name) or
          (pTorPrm->atomType2 == a->name) or
          (pTorPrm->atomType3 == a->name) or
          (pTorPrm->atomType4 == a->name)) {
        torsionsToBeAdded.push_back(pTorPrm);
      }
    }

    for (torsionParamIterator c = torsionsToBeAdded.begin(); c != torsionsToBeAdded.end(); c++) {
      pTorPrm = *c;
      pTP = this->addTorsionParam();
      if (pTorPrm->atomType1 == a->name) {
        pTP->atomType1 = n;
      }
      else {
        pTP->atomType1 = pTorPrm->atomType1;
      }

      if (pTorPrm->atomType2 == a->name) {
        pTP->atomType2 = n;
      }
      else {
        pTP->atomType2 = pTorPrm->atomType2;
      }

      if (pTorPrm->atomType3 == a->name) {
        pTP->atomType3 = n;
      }
      else {
        pTP->atomType3 = pTorPrm->atomType3;
      }

      if (pTorPrm->atomType4 == a->name) {
        pTP->atomType4 = n;
      }
      else {
        pTP->atomType4 = pTorPrm->atomType4;
      }
      pTP->Nt = pTorPrm->Nt;
      pTP->Vn = pTorPrm->Vn;
      pTP->gamma = pTorPrm->gamma;
      pTP->npth = pTorPrm->npth;
      pTP->groupName = g;
    }

    // Copy Improper parameter
    improperParam* pIP = 0;
    improperParam* pImpPrm = 0;
    std::vector<improperParam*> impropersToBeAdded;
    for (improperParamIterator c = this->itsImproperList.begin(); c != this->itsImproperList.end(); c++) {
      pImpPrm = *c;
      if ((pImpPrm->atomType1 == a->name) or
          (pImpPrm->atomType2 == a->name) or
          (pImpPrm->atomType3 == a->name) or
          (pImpPrm->atomType4 == a->name)) {
        impropersToBeAdded.push_back(pImpPrm);
      }
    }

    for (improperParamIterator c = impropersToBeAdded.begin(); c != impropersToBeAdded.end(); c++) {
      pImpPrm = *c;
      pIP = this->addImproperParam();
      if (pImpPrm->atomType1 == a->name) {
        pIP->atomType1 = n;
      }
      else {
        pIP->atomType1 = pImpPrm->atomType1;
      }

      if (pImpPrm->atomType2 == a->name) {
        pIP->atomType2 = n;
      }
      else {
        pIP->atomType2 = pImpPrm->atomType2;
      }

      if (pImpPrm->atomType3 == a->name) {
        pIP->atomType3 = n;
      }
      else {
        pIP->atomType3 = pImpPrm->atomType3;
      }

      if (pImpPrm->atomType4 == a->name) {
        pIP->atomType4 = n;
      }
      else {
        pIP->atomType4 = pImpPrm->atomType4;
      }
      pIP->Nt = pImpPrm->Nt;
      pIP->Vn = pImpPrm->Vn;
      pIP->gamma = pImpPrm->gamma;
      pIP->groupName = g;
    }

    return pAtomType;
}

// ============================================================
// Function : setAtomNumber()
// ------------------------------------------------------------
// 
// ============================================================
void parameters::setAtomNumber(atomType* a)
{
    if (a and pElements) {
      int atNum = pElements->getElementNumber(a->element);
      a->atNum = atNum;
    }
}

// ============================================================
// Function : addBondParam()
// ------------------------------------------------------------
// 
// ============================================================
bondParam* parameters::addBondParam()
{
    pBondParam = new bondParam();
    pBondParam->atomType1 = "";
    pBondParam->atomType2 = "";
    pBondParam->keq = 0.0;
    pBondParam->req = 0.0;
    pBondParam->groupName = "";
    pBondParam->optimize = false;

    this->itsBondList.push_back(pBondParam);
    return pBondParam;
}

// ============================================================
// Function : addAngleParam()
// ------------------------------------------------------------
// 
// ============================================================
angleParam* parameters::addAngleParam()
{
    pAngleParam = new angleParam();
    pAngleParam->atomType1 = "";
    pAngleParam->atomType2 = "";
    pAngleParam->atomType3 = "";
    pAngleParam->keq = 0.0;
    pAngleParam->req = 0.0;
    pAngleParam->groupName = "";
    pAngleParam->optimize = false;

    this->itsAngleList.push_back(pAngleParam);
    return pAngleParam;
}

// ============================================================
// Function : addTorsionParam()
// ------------------------------------------------------------
// 
// ============================================================
torsionParam* parameters::addTorsionParam()
{
    pTorsionParam = new torsionParam();
    pTorsionParam->atomType1 = "";
    pTorsionParam->atomType2 = "";
    pTorsionParam->atomType3 = "";
    pTorsionParam->atomType4 = "";
    pTorsionParam->Nt = 0.0;
    pTorsionParam->Vn = 0.0;
    pTorsionParam->gamma = 0.0;
    pTorsionParam->npth = 0;
    pTorsionParam->groupName = "";

    this->itsTorsionList.push_back(pTorsionParam);
    return pTorsionParam;
}

// ============================================================
// Function : addImproperParam()
// ------------------------------------------------------------
// 
// ============================================================
improperParam* parameters::addImproperParam()
{
    pImproperParam = new improperParam();
    pImproperParam->atomType1 = "";
    pImproperParam->atomType2 = "";
    pImproperParam->atomType3 = "";
    pImproperParam->atomType4 = "";
    pImproperParam->Nt = 0.0;
    pImproperParam->Vn = 0.0;
    pImproperParam->gamma = 0.0;
    pImproperParam->groupName = "";

    this->itsImproperList.push_back(pImproperParam);
    return pImproperParam;
}

// ============================================================
// Function : addHBondParam()
// ------------------------------------------------------------
// 
// ============================================================
hBondParam* parameters::addHBondParam()
{
    pHBondParam = new hBondParam();
    pHBondParam->atomType1 = "";
    pHBondParam->atomType2 = "";
    pHBondParam->p10 = 0.0;
    pHBondParam->p12 = 0.0;
    pHBondParam->groupName = "";

    this->itsHBondList.push_back(pHBondParam);
    return pHBondParam;
}

// ============================================================
// Function : addEquivalentAtomsParam()
// ------------------------------------------------------------
// 
// ============================================================
equivalentAtomsParam* parameters::addEquivalentAtomsParam()
{
    pEquivalentAtomsParam = new equivalentAtomsParam();
    pEquivalentAtomsParam->original = "";
    pEquivalentAtomsParam->groupName = "";

    this->itsEquivalentAtomsList.push_back(pEquivalentAtomsParam);
    return pEquivalentAtomsParam;
}

// ============================================================
// Function : calculateSigmaEpsilon()
// ------------------------------------------------------------
// 
// ============================================================
void parameters::calculateSigmaEpsilon()
{
    double sigmaIJ = 0.0;
    double epsilonIJ = 0.0;
    for (unsigned int i = 0; i < this->itsTypeList.size(); i++) {
      for (unsigned int j = i; j < this->itsTypeList.size(); j++) {
        sigmaIJ   = pow((this->itsTypeList[i]->rvalue + 
                         this->itsTypeList[j]->rvalue),6);

        epsilonIJ = sqrt(this->itsTypeList[i]->evalue * 
                         this->itsTypeList[j]->evalue);

        pLJ612SE = new LJ612SE();
        this->itsLJ612SEList.push_back(pLJ612SE);
        pLJ612SE->atomType1 = this->itsTypeList[i]->name;
        pLJ612SE->atomType2 = this->itsTypeList[j]->name;
        pLJ612SE->sigma = sigmaIJ;
        pLJ612SE->epsilon = epsilonIJ;
      }
    }
    bLJ612SE = true;
}

// ============================================================
// Function : getAtomicNum()
// ------------------------------------------------------------
// 
// ============================================================
int parameters::getAtomicNum(const std::string& at)
{
    for (atomTypeIterator a = this->itsTypeList.begin(); a != this->itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->name == at) {
        return pAtomType->atNum;
      }
    }
    return -1;
}

// ============================================================
// Function : hasAtomType()
// ------------------------------------------------------------
// 
// ============================================================
bool parameters::hasAtomType(const std::string& at, const std::string& g)
{
    for (atomTypeIterator a = this->itsTypeList.begin(); a != this->itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->name == at and pAtomType->groupName == g) {
std::cout << " found " << at << " " << pAtomType->groupName << std::endl;
        return true;
      }
    }
    return false;
}

// ============================================================
// Function : getAtomType()
// ------------------------------------------------------------
// 
// ============================================================
atomType* parameters::getAtomType(const std::string& at)
{
    for (atomTypeIterator a = this->itsTypeList.begin(); a != this->itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->name == at) {
        return pAtomType;
      }
    }
    return 0;
}

// ============================================================
// Function : getAtomTypes()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<atomType*> parameters::getAtomTypes()
{
    return this->itsTypeList;
}

// ============================================================
// Function : getBondParams()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<bondParam*> parameters::getBondParams()
{
    return this->itsBondList;
}

// ============================================================
// Function : getAngleParams()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<angleParam*> parameters::getAngleParams()
{
    return this->itsAngleList;
}

// ============================================================
// Function : getTorsionParams()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<torsionParam*> parameters::getTorsionParams()
{
    return this->itsTorsionList;
}

// ============================================================
// Function : getImproperParams()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<improperParam*> parameters::getImproperParams()
{
    return this->itsImproperList;
}

// ============================================================
// Function : getAtomTypeSymbol()
// ------------------------------------------------------------
// Get element symbol using atom type
// ============================================================
std::string parameters::getAtomTypeSymbol(const std::string& at)
{
    for (atomTypeIterator a = this->itsTypeList.begin(); a != this->itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->name == at) {
        return pAtomType->element;
      }
    }
    return "";
}

// ============================================================
// Function : getAtomTypeHybridization()
// ------------------------------------------------------------
// Get atom hybridization using atom type
// ============================================================
std::string parameters::getAtomTypeHybridization(const std::string& at)
{
    for (atomTypeIterator a = this->itsTypeList.begin(); a != this->itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->name == at) {
        return pAtomType->hybridization;
      }
    }
    return "";
}

// ============================================================
// Function : getBondParam()
// ------------------------------------------------------------
// 
// ============================================================
bondParam* parameters::getBondParam(const std::string& at1, const std::string& at2)
{
    for (bondParamIterator c = this->itsBondList.begin(); c != this->itsBondList.end(); c++) {
      pBondParam = *c;

      if (pBondParam->atomType1 == at1 && pBondParam->atomType2 == at2) {
        return pBondParam;
      }
      if (pBondParam->atomType1 == at2 && pBondParam->atomType2 == at1) {
        return pBondParam;
      }
    }
    return 0;
}

// ============================================================
// Function : hasBondParam()
// ------------------------------------------------------------
// 
// ============================================================
bool parameters::hasBondParam(const std::string& at1, const std::string& at2, const std::string& g)
{
    for (bondParamIterator c = this->itsBondList.begin(); c != this->itsBondList.end(); c++) {
      pBondParam = *c;

      if (pBondParam->atomType1 == at1 && pBondParam->atomType2 == at2) {
        if (pBondParam->groupName == g) {
          return true;
        }
      }
      if (pBondParam->atomType1 == at2 && pBondParam->atomType2 == at1) {
        if (pBondParam->groupName == g) {
          return true;
        }      }
    }
    return false;
}

// ============================================================
// Function : getAngleParam()
// ------------------------------------------------------------
// 
// ============================================================
angleParam* parameters::getAngleParam(const std::string& at1, const std::string& at2, const std::string& at3)
{
    for (angleParamIterator c = this->itsAngleList.begin(); c != this->itsAngleList.end(); c++) {
      pAngleParam = *c;
      if (pAngleParam->atomType1 == at1 &&
          pAngleParam->atomType2 == at2 &&
          pAngleParam->atomType3 == at3) {
        return pAngleParam;
      }
      if (pAngleParam->atomType1 == at3 &&
          pAngleParam->atomType2 == at2 &&
          pAngleParam->atomType3 == at1) {
        return pAngleParam;
      }
    }
    return 0;
}

// ============================================================
// Function : hasAngleParam()
// ------------------------------------------------------------
// 
// ============================================================
bool parameters::hasAngleParam(const std::string& at1, const std::string& at2, 
	                           const std::string& at3, const std::string& g)
{
    for (angleParamIterator c = this->itsAngleList.begin(); c != this->itsAngleList.end(); c++) {
      pAngleParam = *c;
      if (pAngleParam->atomType1 == at1 &&
          pAngleParam->atomType2 == at2 &&
          pAngleParam->atomType3 == at3) {
        if (pAngleParam->groupName == g) {
          return true;
        }
      }
      if (pAngleParam->atomType1 == at3 &&
          pAngleParam->atomType2 == at2 &&
          pAngleParam->atomType3 == at1) {
        if (pAngleParam->groupName == g) {
          return true;
        }      }
    }
    return false;
}

// ============================================================
// Function : getTorsionParamList()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<torsionParam*> parameters::getTorsionParamList(
                    const std::string& at1, const std::string& at2,
                    const std::string& at3, const std::string& at4)
{
    std::vector<torsionParam*> torsionParamList;
    std::vector<torsionParam*>::iterator vecIterator;
    for (torsionParamIterator c = this->itsTorsionList.begin(); c != this->itsTorsionList.end(); c++) {
      pTorsionParam = *c;

      if ( (pTorsionParam->atomType2 == at2 && pTorsionParam->atomType3 == at3) ||
           (pTorsionParam->atomType3 == at2 && pTorsionParam->atomType2 == at3)) {
        if ((pTorsionParam->atomType1 == at1 && pTorsionParam->atomType4 == at4)) {
          vecIterator = std::find(torsionParamList.begin(), torsionParamList.end(), pTorsionParam);
          if (vecIterator == torsionParamList.end()) {
            torsionParamList.push_back(pTorsionParam);
          }
        }
        else if ((pTorsionParam->atomType1 == at4 && pTorsionParam->atomType4 == at1)) {
          vecIterator = std::find(torsionParamList.begin(), torsionParamList.end(), pTorsionParam);
          if (vecIterator == torsionParamList.end()) {
            torsionParamList.push_back(pTorsionParam);
          }
        }
      }
    }

    // If no SPECIFIC (A-B-C-D) torsion parameter exists then search for 
    // a GENERAL (X-A-B-X ) torsion parameter (JACS 106, 765-784)
    if (torsionParamList.size() == 0) {
      for (torsionParamIterator c = this->itsTorsionList.begin(); c != this->itsTorsionList.end(); c++) {
        pTorsionParam = *c;
        if ( (pTorsionParam->atomType2 == at2 && pTorsionParam->atomType3 == at3) ||
             (pTorsionParam->atomType3 == at2 && pTorsionParam->atomType2 == at3)) {
          if ((pTorsionParam->atomType1 == "X" && pTorsionParam->atomType4 == "X")) {
            vecIterator = std::find(torsionParamList.begin(), torsionParamList.end(), pTorsionParam);
            if (vecIterator == torsionParamList.end()) {
              torsionParamList.push_back(pTorsionParam);
            }
          }
        }
      }
    }
    return torsionParamList;
}

// ============================================================
// Function : removeProlineTorsion()
// ------------------------------------------------------------
// 
// ============================================================
void parameters::removeProlineTorsion(std::vector<torsionParam*>& t)
{
    std::vector<torsionParam*> toBeDeleted;
    std::vector<torsionParam*>::iterator vecIterator;

    for (unsigned int k = 0; k < t.size(); k++ ) {
      pTorsionParam = t[k];
      if ((pTorsionParam->atomType1 == "CT") and (pTorsionParam->atomType2 == "CT") and
          (pTorsionParam->atomType3 == "N" ) and (pTorsionParam->atomType4 == "C" ) ) {
        if ((pTorsionParam->Nt == -4) or (pTorsionParam->Nt == 1)) {
          toBeDeleted.push_back(pTorsionParam);
        }
      }
    }

    if (!toBeDeleted.empty()) {
      std::sort( toBeDeleted.begin(), toBeDeleted.end() );
      for (unsigned int n = 0; n < toBeDeleted.size(); n++) {
        vecIterator = std::find(t.begin(), t.end(),toBeDeleted[n]);
        t.erase(vecIterator);
      }
    }
    toBeDeleted.clear();
}

// ============================================================
// Function : getImproperParamList()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<improperParam*> parameters::getImproperParamList(
                            const std::string& at1,const std::string& at2,
                            const std::string& at3,const std::string& at4,
                            std::vector<std::vector<int> >& order)
{
    int possibleImpropers[11][4] = {{0,1,2,3}, {0,1,3,2}, {0,3,2,1},
                                   {0,3,1,2}, {1,0,2,3}, {1,0,3,2},
                                   {1,3,2,0}, {3,0,2,1}, {3,1,2,0},
                                   {2,3,1,0}, {0,2,1,3}} ;

    bool gotIt = false;
    std::string improperTypes[4] = {at1, at2, at3, at4};
    std::vector<improperParam*> improperParamList;

    for (improperParamIterator c = this->itsImproperList.begin(); c != this->itsImproperList.end(); c++) {
      pImproperParam = *c;
      if (gotIt) {
        return improperParamList;
      }
      // X-A-B-C
      if (pImproperParam->atomType1 == "X" && pImproperParam->atomType2 != "X") {
        for (unsigned int i = 0; i < 11; i++) {
          if ((improperTypes[possibleImpropers[i][1]] == pImproperParam->atomType2) &&
              (improperTypes[possibleImpropers[i][2]] == pImproperParam->atomType3) &&
              (improperTypes[possibleImpropers[i][3]] == pImproperParam->atomType4)) {
            improperParamList.push_back(pImproperParam);
            std::vector <int> temp;
            for (unsigned int x = 0; x < 4; x++) {
              temp.push_back(possibleImpropers[i][x]);
            }
            order.push_back(temp);
            gotIt = true;
            break;
          }
        }
      }

      // X-X-A-B
      else if (pImproperParam->atomType1 == "X" && pImproperParam->atomType2 == "X" ) {
        for (unsigned int i = 0; i < 11; i++) {
          if ((improperTypes[possibleImpropers[i][2]] == pImproperParam->atomType3) && 
              (improperTypes[possibleImpropers[i][3]] == pImproperParam->atomType4)) {
            improperParamList.push_back(pImproperParam);
            std::vector <int> temp;
            for (unsigned int x = 0; x < 4; x++) {
              temp.push_back(possibleImpropers[i][x]);
            }
            order.push_back(temp);
            gotIt = true;
            break;
          }
        }
      }
      // A-B-C-D
      else {
        for (unsigned int i = 0; i < 11; i++) {
          if ((improperTypes[possibleImpropers[i][0]] == pImproperParam->atomType1) &&
              (improperTypes[possibleImpropers[i][1]] == pImproperParam->atomType2) &&
              (improperTypes[possibleImpropers[i][2]] == pImproperParam->atomType3) &&
              (improperTypes[possibleImpropers[i][3]] == pImproperParam->atomType4)) {
            improperParamList.push_back(pImproperParam);
            std::vector <int> temp;
            for (unsigned int x = 0; x < 4; x++) {
              temp.push_back(possibleImpropers[i][x]);
            }
            order.push_back( temp );
            gotIt = true;
            break;
          }
        }
      }
    }
    return improperParamList;
}

// ============================================================
// Function : getLJ612SE()
// ------------------------------------------------------------
// 
// ============================================================
LJ612SE* parameters::getLJ612SE(const std::string& at1, const std::string& at2)
{
    if (!bLJ612SE) {
      this->calculateSigmaEpsilon();
    }

    for (LJ612SEIterator c = this->itsLJ612SEList.begin(); c != this->itsLJ612SEList.end(); c++) {
      pLJ612SE = *c;
      if (pLJ612SE->atomType1 == at1 && pLJ612SE->atomType2 == at2) {
        return pLJ612SE;
      }
      if (pLJ612SE->atomType1 == at2 && pLJ612SE->atomType2 == at1) {
        return pLJ612SE;
      }
    }
    return 0;
}

// ============================================================
// Function : getEquivalentAtom()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<equivalentAtomsParam*> parameters::getEquivalentAtomList()
{
    return itsEquivalentAtomsList;
}

} // MTKpp namespace
