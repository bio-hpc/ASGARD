/*!
   \file molecule.cpp
   \brief Container for submolecules, bonds, angles, torsions, and impropers
   \author Martin Peters

   Container for submolecules, bonds, angles, torsions, and impropers

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.42 $

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

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "Utils/vector3d.h"

#include "bond.h"
#include "angle.h"
#include "torsion.h"
#include "improper.h"
#include "ring.h"
#include "conformer.h"
#include "protonate.h"
#include "fingerPrint.h"
#include "functionalize.h"
#include "hydrophobize.h"
#include "hybridize.h"
#include "pharmacophore.h"
#include "stdLibrary.h"
#include "stdFrag.h"
#include "parameters.h"
#include "utility.h"

#include "Diagnostics/MTKException.h"
#include "Utils/index.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : molecule()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
molecule::molecule(collection *parent):pParent(parent)
{
    // FLAGS
    inFileType = "";
    bBondsAssigned = false;
    bAnglesAssigned = false;
    bTorsionsAssigned = false;
    bImpropersAssigned = false;
    bBondTypes1Assigned = false;
    bBondTypes2Assigned = false;
    bMMAtomTypesAssigned = false;
    bAtomValencesAssigned = false;
    bAtomHybridizationsAssigned = false;
    bRingsAssigned = false;
    bRingsKekulized = false;
    bInSolution = false;
    bHydrogensAdded = false;
    bError = false;

    // variables
    itsName            = "";
    itsChain           = "";
    itsKind            = 0;
    itsMolId           = 0;
    itsNumAtoms        = 0;
    itsNumBonds        = 0;
    totalCharge        = 0;
    itsSMolIndex       = 1;
    itsAtomIndex       = 1; // must start at one for bonding
    this->itsMaxFileID = 1;

    // pointers
    pAtom              = 0;
    pBond              = 0;
    pSubMolecule       = 0;
    pAngle             = 0;
    pTorsion           = 0;
    pImproper          = 0;
    pParameters        = 0;
    pRings             = 0;
    pConformers        = 0;
    pProtonate         = 0;
    pFingerPrint       = 0;
    pFunctionalize     = 0;
    pHydrophobize      = 0;
    pPharmacophore     = 0;
    adjMatrix          = 0;
    adjMatrixSize      = 0;
    featureDistMatrix  = 0;

    //std::string errorMessage = "MAXATOMS: " + i2s(MAXATOMS);
    //errorLogger.throwError("molecule::construct", errorMessage, INFO);
    //pParameters = pParent->getParameters();
}

// ============================================================
// Function : ~molecule()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
molecule::~molecule()
{
    errorLogger.throwError(" Deleting molecule: ", this->getName(), INFO);

    for (sMolIterator c = this->itsSubMoleculeList.begin(); c != this->itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      delete pSubMolecule;
    }
    this->itsSubMoleculeList.clear();

    this->itsAtomList.erase(this->itsAtomList.begin(), this->itsAtomList.end());
    this->itsAtomList.clear();

    for (BondMapIterator b = this->itsBondMap.begin(); b != this->itsBondMap.end(); b++) {
      pBond = b->second;
      delete pBond;
    }
    this->itsBondMap.clear();

    for (AngleMapIterator a = this->itsAngleMap.begin(); a != this->itsAngleMap.end(); a++) {
      pAngle = a->second;
      delete pAngle;
    }
    this->itsAngleMap.clear();

    for (TorsionMapIterator t = this->itsTorsionMap.begin(); t != this->itsTorsionMap.end(); t++) {
      pTorsion = t->second;
      delete pTorsion;
    }
    this->itsTorsionMap.clear();

    for (ImproperMapIterator i = this->itsImproperMap.begin(); i != this->itsImproperMap.end(); i++) {
      pImproper = i->second;
      delete pImproper;
    }
    this->itsImproperMap.clear();

    for (ConformerIterator c = this->itsConformers.begin(); c != this->itsConformers.end(); c++) {
      pConformer = *c;
      delete pConformer;
    }
    this->itsConformers.clear();

    for (RingIterator r = this->itsRings.begin(); r != this->itsRings.end(); r++) {
      pRing = *r;
      delete pRing;
    }
    this->itsRings.clear();

    for (HydrophobeIterator h = this->itsHydrophobes.begin(); h != this->itsHydrophobes.end(); h++) {
      pHydrophobe = *h;
      delete pHydrophobe;
    }
    this->itsHydrophobes.clear();

    this->itsBond13List.clear();
    this->itsBond14List.clear();
    delete [] adjMatrix;
    delete [] featureDistMatrix;

    delete pRings;
    delete pConformers;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
//
// ============================================================
collection* molecule::getParent()
{
    return pParent;
}

// ============================================================
// Function : addSubMolecule()
// ------------------------------------------------------------
//
// ============================================================
submolecule* molecule::addSubMolecule()
{
    pSubMolecule = new submolecule(this);
    pSubMolecule->setIndex(this->itsSMolIndex);

    int iSmolIndex = pParent->getSubMoleculeIndex();

    pSubMolecule->setColIndex(iSmolIndex);
    pParent->setSubMoleculeIndex(iSmolIndex+1);

    itsSMolIndex++;
    this->itsSubMoleculeList.push_back(pSubMolecule);
    return pSubMolecule;
}

// ============================================================
// Function : getNumSubMolecules()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumSubMolecules()
{
    //errorLogger.throwError("molecule::getNumSubMolecules", " Getting Number of Submolecules ", INFO);
    return this->itsSubMoleculeList.size();
}

// ============================================================
// Function : addBond()
// ------------------------------------------------------------
//
// ============================================================
Bond* molecule::addBond(atom* at1, atom* at2, const int &type,
                const int &stereo, const int &topology, const double &b)
{
    if (at1 and at2) {
      // make sure at1 and at2 are in the same molecule
      if (at1->getParent()->getParent() != at1->getParent()->getParent()) {
        errorLogger.throwError("molecule::addBond", " Trying to Add a Bond between two molecules ", INFO);
        return 0;
      }
      if (!this->hasBond(at1, at2)) {
        pBond = new Bond();
        pBond->atom1 = at1;
        pBond->atom2 = at2;
        pBond->type  = type;
        pBond->stereo = stereo;
        pBond->topology = topology;
        pBond->size  = b;
        pBond->rotatable = 1;

        int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);
        itsBondMap[bondIndex] = pBond;
        return pBond;
      }
    }
    return 0;
}

// ============================================================
// Function : getBond()
// ------------------------------------------------------------
//
// ============================================================
Bond* molecule::getBond(atom* at1, atom* at2)
{
    if (at1 == at2) return 0;

    if (!at1 or !at2) {
      errorLogger.throwError("molecule::getBond", " Can't find one or both atoms ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "molecule::getBond" << " Can't find one or both atoms ";
      throw MTKException(ss.str());
    }

    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);

    BondMapIterator b = itsBondMap.find(bondIndex);

    if (b != itsBondMap.end()) {
      return itsBondMap[bondIndex];
    }
    return 0;
}

// ============================================================
// Function : getBonds()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getBonds(int bonds[])
{
    int i = 0;
    try {
      if (!this->itsBondMap.empty()) {
        for (BondMapIterator b = this->itsBondMap.begin();
             b != this->itsBondMap.end(); b++) {
          pBond = b->second;
          bonds[i] = pBond->atom1->getIndex()-1;
          bonds[i+1] = pBond->atom2->getIndex()-1;
          i+=2;
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getBonds", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getBondParams()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getBondParams(double bondParams[])
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getBondParams", " Can't find required parameters ", MTK_ERROR);
      std::stringstream ss;
      ss << "molecule::getBondParams" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    int i = 0;
    try {
      if (!this->itsBondMap.empty()) {
        for (BondMapIterator b = this->itsBondMap.begin();
             b != this->itsBondMap.end(); b++) {
          pBond = b->second;
          pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                                                 pBond->atom2->getStdAtom()->type);

          if (pBondParam) {
            //std::cout << pBond->atom1->getStdAtom()->type << " " << pBond->atom2->getStdAtom()->type << " "
            //          << pBondParam->req << " " << pBondParam->keq << std::endl;
            bondParams[i]   = pBondParam->req;
            bondParams[i+1] = pBondParam->keq;

            i+=2;
          }
          else {
            errorLogger.throwError("molecule::getBondParams", " Setting Bond Parameters to 0 ", MTK_ERROR);
            bondParams[i]   = 0.0;
            bondParams[i+1] = 0.0;
            i+=2;
          }
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getBondParams", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : add13Bond()
// ------------------------------------------------------------
//
// ============================================================
void molecule::add13Bond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);
    this->itsBond13List.push_back(bondIndex);
}

// ============================================================
// Function : add14Bond()
// ------------------------------------------------------------
//
// ============================================================
void molecule::add14Bond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);
    this->itsBond14List.push_back(bondIndex);
}

// ============================================================
// Function : addAngle()
// ------------------------------------------------------------
//
// ============================================================
Angle* molecule::addAngle(atom* at1, atom* at2, atom* at3, const double &ang)
{
    pAngle = new Angle();
    pAngle->atom1 = at1;
    pAngle->atom2 = at2;
    pAngle->atom3 = at3;
    pAngle->size = ang;

    try {
      ULONG_KIND angleIndex = indexABC_ULL(at1->getIndex(), at2->getIndex(),
                                at3->getIndex(), MAXATOMS, MAXATOMS);
/*
if (angleIndex < 0) {
std::cout << at1->getIndex() << "-"
          << at2->getIndex() << "-"
          << at3->getIndex() << "  "
          << MAXATOMS << "   "
          << angleIndex << "  new index = "
          << angleIndex2
          << std::endl;
exit(1);
}
*/
      this->itsAngleMap[angleIndex] = pAngle;
      return pAngle;
    }
    catch (MTKException& e) {
      std::cout << "Error in addAngle: " << e.message << std::endl;
    }
    return 0;
}

// ============================================================
// Function : getAngles()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getAngles(int angles[])
{
    int i = 0;
    try {
      if (!this->itsAngleMap.empty()){
        for (AngleMapIterator b = this->itsAngleMap.begin();
             b != this->itsAngleMap.end(); b++) {
          pAngle = b->second;
          angles[i]   = pAngle->atom1->getIndex()-1;
          angles[i+1] = pAngle->atom2->getIndex()-1;
          angles[i+2] = pAngle->atom3->getIndex()-1;
          i+=3;
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getAngles", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : addTorsion()
// ------------------------------------------------------------
//
// ============================================================
Torsion* molecule::addTorsion(atom* at1, atom* at2, atom* at3, atom* at4,
                              const double &tor)
{
    pTorsion = new Torsion();
    pTorsion->atom1 = at1;
    pTorsion->atom2 = at2;
    pTorsion->atom3 = at3;
    pTorsion->atom4 = at4;
    pTorsion->size = tor;
    pTorsion->parametersAssigned = false;

    try {
      ULONG_KIND torsionIndex = indexABCD_ULL(at1->getIndex(), at2->getIndex(),
                                              at3->getIndex(), at4->getIndex(),
                                              MAXATOMS, MAXATOMS, MAXATOMS);

      this->itsTorsionMap[torsionIndex] = pTorsion;
    }
    catch (MTKException& e) {
      std::cout << "Error in addTorsion: " << e.message << std::endl;
    }
    return pTorsion;
}

// ============================================================
// Function : getTorsions()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getTorsions(int torsions[])
{
    int i = 0;
    try {
      if (!this->itsTorsionMap.empty()){
        for (TorsionMapIterator b = this->itsTorsionMap.begin();
             b != this->itsTorsionMap.end(); b++) {
          pTorsion = b->second;
          torsions[i]   = pTorsion->atom1->getIndex()-1;
          torsions[i+1] = pTorsion->atom2->getIndex()-1;
          torsions[i+2] = pTorsion->atom3->getIndex()-1;
          torsions[i+3] = pTorsion->atom4->getIndex()-1;
          i += 4;
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getTorsions", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : addImproper()
// ------------------------------------------------------------
//
// ============================================================
Improper* molecule::addImproper(atom* at1, atom* at2, atom* at3, atom* at4, const double &imp)
{
    pImproper = new Improper();
    pImproper->atom1 = at1;
    pImproper->atom2 = at2;
    pImproper->atom3 = at3;
    pImproper->atom4 = at4;
    pImproper->size = imp;
    pImproper->parametersAssigned = false;

    int improperIndex = indexABCDb(at1->getIndex(), at2->getIndex(),
                                   at3->getIndex(), at4->getIndex(),
                                   MAXATOMS, MAXATOMS, MAXATOMS);
    this->itsImproperMap[improperIndex] = pImproper;
    return pImproper;
}

// ============================================================
// Function : getImpropers()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getImpropers(int impropers[])
{
    int i = 0;
    try {
      if (!this->itsImproperMap.empty()){
        for (ImproperMapIterator b = this->itsImproperMap.begin();
             b != this->itsImproperMap.end(); b++) {
          pImproper      = b->second;
          impropers[i]   = pImproper->atom1->getIndex()-1;
          impropers[i+1] = pImproper->atom2->getIndex()-1;
          impropers[i+2] = pImproper->atom3->getIndex()-1;
          impropers[i+3] = pImproper->atom4->getIndex()-1;
          i+=4;
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getImpropers", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : hasBond()
// ------------------------------------------------------------
//
// ============================================================
bool molecule::hasBond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);
    BondMapIterator b = this->itsBondMap.find(bondIndex);

    if (b != this->itsBondMap.end()) {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : has13Bond()
// ------------------------------------------------------------
//
// ============================================================
bool molecule::has13Bond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);
    std::vector<int>::iterator result;
    result = std::find(this->itsBond13List.begin(), this->itsBond13List.end(),
                  bondIndex);

    if (result != this->itsBond13List.end()) {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : has14Bond()
// ------------------------------------------------------------
//
// ============================================================
bool molecule::has14Bond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);
    std::vector<int>::iterator result;
    result = std::find(this->itsBond14List.begin(),
                  this->itsBond14List.end(), bondIndex);

    if (result != this->itsBond14List.end()) {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : hasAngle()
// ------------------------------------------------------------
//
// ============================================================
bool molecule::hasAngle(atom* at1, atom* at2, atom* at3)
{
    try {
      ULONG_KIND angleIndex = indexABC_ULL(at1->getIndex(), at2->getIndex(),
                                at3->getIndex(),
                                MAXATOMS, MAXATOMS);

      AngleMapIterator a = this->itsAngleMap.find(angleIndex);

      if (a != this->itsAngleMap.end()) {
        return 1;
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in hasAngle: " << e.message << std::endl;
    }
    return 0;
}

// ============================================================
// Function : getAngle()
// ------------------------------------------------------------
//
// ============================================================
Angle* molecule::getAngle(atom* at1, atom* at2, atom* at3)
{
    try {
      ULONG_KIND angleIndex = indexABC_ULL(at1->getIndex(), at2->getIndex(),
                                at3->getIndex(),
                                MAXATOMS,MAXATOMS);

      AngleMapIterator a = this->itsAngleMap.find(angleIndex);

      if (a != itsAngleMap.end()) {
        return this->itsAngleMap[angleIndex];
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in getAngle: " << e.message << std::endl;
    }
    return 0;
}

// ============================================================
// Function : hasTorsion()
// ------------------------------------------------------------
//
// ============================================================
bool molecule::hasTorsion(atom* at1, atom* at2, atom* at3, atom* at4)
{
    try {
      ULONG_KIND torsionIndex = indexABCD_ULL(at1->getIndex(),
                             at2->getIndex(), at3->getIndex(), at4->getIndex(),
                             MAXATOMS,MAXATOMS,MAXATOMS);

      TorsionMapIterator t = this->itsTorsionMap.find(torsionIndex);

      if (t != this->itsTorsionMap.end()) {
        return 1;
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in hasTorsion: " << e.message << std::endl;
    }
    return 0;
}

// ============================================================
// Function : getTorsion()
// ------------------------------------------------------------
//
// ============================================================
Torsion* molecule::getTorsion(atom* at1, atom* at2, atom* at3, atom* at4)
{
    try {
      ULONG_KIND torsionIndex = indexABCD_ULL(at1->getIndex(),
                             at2->getIndex(), at3->getIndex(), at4->getIndex(),
                             MAXATOMS,MAXATOMS,MAXATOMS);

      TorsionMapIterator t = this->itsTorsionMap.find(torsionIndex);

      if (t != this->itsTorsionMap.end()) {
        return this->itsTorsionMap[torsionIndex];
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in getTorsion: " << e.message << std::endl;
    }
    return 0;
}

// ============================================================
// Function : setTorsion()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setTorsion(atom* at1, atom* at2, atom* at3, atom* at4, double v)
{
    try {
      ULONG_KIND torsionIndex = indexABCD_ULL(at1->getIndex(),
                             at2->getIndex(), at3->getIndex(), at4->getIndex(),
                             MAXATOMS,MAXATOMS,MAXATOMS);

      TorsionMapIterator t = this->itsTorsionMap.find(torsionIndex);

      if (t != this->itsTorsionMap.end()) {
        std::vector<atom*> aList = this->getAtomList();
        conformers* pConf = new conformers(this);
        pConf->addRotatableTorsion(at1->getIndex(), at2->getIndex(),
                                 at3->getIndex(), at4->getIndex());
        pConf->determineAtomList();

        conformer* myConf = new conformer();
        myConf->torsions.push_back(v * DEG2RAD);
        myConf->name = "setTorsion";

        std::vector< vector3d > newCoords;
        pConf->getConformerCoordinates(myConf, newCoords);

        for (unsigned int i = 0; i < newCoords.size(); i++) {
          aList[i]->setCoords(newCoords[i][0], newCoords[i][1], newCoords[i][2]);
        }
        delete myConf;
        delete pConf;
      }
      else {
        errorLogger.throwError("molecule::setTorsion", " Can't find torsion ", MTK_ERROR);
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in setTorsion: " << e.message << std::endl;
    }
}

// ============================================================
// Function : hasImproper()
// ------------------------------------------------------------
//
// ============================================================
bool molecule::hasImproper(atom* at1, atom* at2, atom* at3, atom* at4)
{
    int improperIndex = indexABCDb(at1->getIndex(), at2->getIndex(),
                                  at3->getIndex(), at4->getIndex(),
                                  MAXATOMS,MAXATOMS,MAXATOMS);

    ImproperMapIterator i = this->itsImproperMap.find(improperIndex);

    if (i != this->itsImproperMap.end()){
      return 1;
    }
    return 0;
}

// ============================================================
// Function : delBond()
// ------------------------------------------------------------
//
// ============================================================
void molecule::delBond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getIndex(), at2->getIndex(), MAXATOMS);

    BondMapIterator b = this->itsBondMap.find(bondIndex);

    if (b != this->itsBondMap.end()) {
      this->itsBondMap.erase(b);
    }
}

// ============================================================
// Function : delAngle()
// ------------------------------------------------------------
//
// ============================================================
void molecule::delAngle(atom* at1, atom* at2, atom* at3)
{
    try {
      ULONG_KIND angleIndex = indexABC_ULL(at1->getIndex(), at2->getIndex(),
                                at3->getIndex(),
                                MAXATOMS,MAXATOMS);

      AngleMapIterator a = this->itsAngleMap.find(angleIndex);

      if (a != this->itsAngleMap.end()) {
        this->itsAngleMap.erase(a);
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in delAngle: " << e.message << std::endl;
    }
}

// ============================================================
// Function : delTorsion()
// ------------------------------------------------------------
//
// ============================================================
void molecule::delTorsion(atom* at1, atom* at2, atom* at3, atom* at4)
{
    try {
      ULONG_KIND torsionIndex = indexABCD_ULL(at1->getIndex(),
                    at2->getIndex(), at3->getIndex(), at4->getIndex(),
                    MAXATOMS, MAXATOMS,MAXATOMS);

      TorsionMapIterator t = this->itsTorsionMap.find(torsionIndex);

      if (t != this->itsTorsionMap.end()) {
        this->itsTorsionMap.erase(t);
      }
    }
    catch (MTKException& e) {
      std::cout << "Error in delTorsion: " << e.message << std::endl;
    }
}

// ============================================================
// Function : delImproper()
// ------------------------------------------------------------
//
// ============================================================
void molecule::delImproper(atom* at1, atom* at2, atom* at3, atom* at4)
{
    int improperIndex = indexABCDb(at1->getIndex(), at2->getIndex(),
                        at3->getIndex(), at4->getIndex(),
                        MAXATOMS, MAXATOMS, MAXATOMS);

    ImproperMapIterator i = this->itsImproperMap.find(improperIndex);

    if (i != this->itsImproperMap.end()) {
      this->itsImproperMap.erase(i);
    }
}

// ============================================================
// Function : numBonds()
// ------------------------------------------------------------
//
// ============================================================
int molecule::numBonds()
{
    if (!this->itsBondMap.empty()) {
      return this->itsBondMap.size();
    }
    return 0;
}

// ============================================================
// Function : numAngles()
// ------------------------------------------------------------
//
// ============================================================
int molecule::numAngles()
{
    if (!this->itsAngleMap.empty()) {
      return this->itsAngleMap.size();
    }
    return 0;

}

// ============================================================
// Function : numTorsions()
// ------------------------------------------------------------
//
// ============================================================
int molecule::numTorsions()
{
    return this->itsTorsionMap.size();
}

// ============================================================
// Function : numImpropers()
// ------------------------------------------------------------
//
// ============================================================
int molecule::numImpropers()
{
    return this->itsImproperMap.size();
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setName(const std::string &name)
{
     this->itsName = name;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
//
// ============================================================
std::string molecule::getName()
{
     return itsName;
}

// ============================================================
// Function : setChain()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setChain(const std::string &chain)
{
     this->itsChain = chain;
}

// ============================================================
// Function : getChain()
// ------------------------------------------------------------
//
// ============================================================
std::string molecule::getChain()
{
     return this->itsChain;
}

// ============================================================
// Function : setKind()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setKind(const int &n)
{
     this->itsKind = n;
}

// ============================================================
// Function : getKind()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getKind()
{
     return this->itsKind;
}

// ============================================================
// Function : setMolId()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setMolId(const int &id)
{
     itsMolId = id;
}

// ============================================================
// Function : getMolId()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMolId()
{
     return itsMolId;
}

// ============================================================
// Function : setNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setNumAtoms(const int &natoms)
{
     itsNumAtoms = natoms;
}

// ============================================================
// Function : getNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumAtoms()
{
    int nAtoms = 0;
    if (this->itsSubMoleculeList.empty()) {
      return 0;
    }
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        nAtoms++;
      }
    }
    return nAtoms;
}

// ============================================================
// Function : getNumHeavyAtoms()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumHeavyAtoms()
{
    int nHeavyAtoms = 0;
    if (this->itsSubMoleculeList.empty()) {
      return 0;
    }
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        pAtom = *d;
        if (pAtom->getElementSymbol() != "H") {
          nHeavyAtoms++;
        }
      }
    }
    return nHeavyAtoms;
}

// ============================================================
// Function : setNumBonds()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setNumBonds(const int &nbonds)
{
     itsNumBonds = nbonds;
}

// ============================================================
// Function : getNumBonds()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumBonds()
{
    if (this->itsBondMap.empty()) {
      return 0;
    }
    return this->itsBondMap.size();
}

// ============================================================
// Function : getNumElectrons()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumElectrons()
{
    int itsNumElectrons = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
       pSubMolecule = *c;
       std::vector<atom*> atList = pSubMolecule->getAtomList();
       for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
         pAtom = *d;
         itsNumElectrons += pAtom->getElement()->number;
       }
    }
    return itsNumElectrons;
}

// ============================================================
// Function : setAtomIndex()
// ------------------------------------------------------------
// Used in submolecule when an atom is added
// ============================================================
void molecule::setAtomIndex(const int &n)
{
    this->itsAtomIndex = n;
}

// ============================================================
// Function : getAtomIndex()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getAtomIndex()
{
    return this->itsAtomIndex;
}

// ============================================================
// Function : setMaxFileID()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setMaxFileID(const int &n)
{
    if (n > this->itsMaxFileID) {
      this->itsMaxFileID = n;
    }
}

// ============================================================
// Function : getMaxFileID()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMaxFileID()
{
    return itsMaxFileID;
}

// ============================================================
// Function : getSubMoleculeList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<submolecule*> molecule::getSubMoleculeList()
{
    return itsSubMoleculeList;
}

// ============================================================
// Function : getSubMoleculeList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<submolecule*> molecule::getSubMoleculeList(std::string name)
{
    std::vector<submolecule*> smolList;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
       pSubMolecule = *c;
       if (pSubMolecule->getName() == name) {
         smolList.push_back(pSubMolecule);
       }
    }
    return smolList;
}

// ============================================================
// Function : getSubMoleculeList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<submolecule*> molecule::getSubMoleculeList(std::string name, int num)
{
    std::vector<submolecule*> smolList;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
       pSubMolecule = *c;
       if ((pSubMolecule->getSubMolId() == num) and (pSubMolecule->getName() == name)) {
         smolList.push_back(pSubMolecule);
       }
    }
    return smolList;
}

// ============================================================
// Function : getSubMoleculeList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<submolecule*> molecule::getSubMoleculeList(int num)
{
    std::vector<submolecule*> smolList;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
       pSubMolecule = *c;
       if (pSubMolecule->getSubMolId() == num) {
         smolList.push_back(pSubMolecule);
       }
    }
    return smolList;
}

// ============================================================
// Function : getAtomList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> molecule::getAtomList()
{
    this->itsAtomList.clear();
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        pAtom = *d;
        this->itsAtomList.push_back(pAtom);
      }
    }
    return itsAtomList;
}

// ============================================================
// Function : getAtomList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> molecule::getHeavyAtomList()
{
    this->itsHeavyAtomList.clear();
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        pAtom = *d;
        if (pAtom->getElementSymbol() == "H") continue;
        this->itsHeavyAtomList.push_back(pAtom);
      }
    }
    return this->itsHeavyAtomList;
}

// ============================================================
// Function : getBondMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<int, Bond*> molecule::getBondMap()
{
    return itsBondMap;
}

// ============================================================
// Function : getAngleMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<ULONG_KIND, Angle*> molecule::getAngleMap()
{
    return itsAngleMap;
}

// ============================================================
// Function : getAngleParams()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getAngleParams(double angleParams[])
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getAngleParams", " Can't find required parameters ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "molecule::getAngleParams" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    int i = 0;
    try {
      if (!this->itsAngleMap.empty()){
        for (AngleMapIterator b = this->itsAngleMap.begin(); b != this->itsAngleMap.end(); b++) {
          pAngle = b->second;
          pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                   pAngle->atom2->getStdAtom()->type,
                                                   pAngle->atom3->getStdAtom()->type);
          if (pAngleParam) {
            angleParams[i]   = pAngleParam->req;
            angleParams[i+1] = pAngleParam->keq;
            i+=2;
          }
          else {
            angleParams[i]   = 0.0;
            angleParams[i+1] = 0.0;
            i+=2;
          }
        }
      }
    }
    catch (std::exception& e) {
      std::string errMessage = e.what();
      errorLogger.throwError("molecule::getAngleParams", errMessage, MTK_ERROR);
    }
    return 0;
}

// ============================================================
// Function : getTorsionMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<ULONG_KIND, Torsion*> molecule::getTorsionMap()
{
    return this->itsTorsionMap;
}

// ============================================================
// Function : getNumMMTorsions()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumMMTorsions()
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getNumMMTorsions", " Can't find required parameters ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "molecule::getNumMMTorsions" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    //int i = 0;
    int n = 0;

    if (!itsTorsionMap.empty()) {
      for (TorsionMapIterator b = itsTorsionMap.begin(); b != itsTorsionMap.end(); b++) {
        pTorsion = b->second;
        stdAtom* pStdAt1 = pTorsion->atom1->getStdAtom();
        stdAtom* pStdAt2 = pTorsion->atom2->getStdAtom();
        stdAtom* pStdAt3 = pTorsion->atom3->getStdAtom();
        stdAtom* pStdAt4 = pTorsion->atom4->getStdAtom();
        if (pStdAt1 and pStdAt2 and pStdAt3 and pStdAt4) {
          torsionParamList = pParameters->getTorsionParamList(
                             pStdAt1->type, pStdAt2->type, pStdAt3->type, pStdAt4->type);
          if (!torsionParamList.empty()) {
            //std::cout << pTorsion->atom1->getIndex()-1 << "-" << pTorsion->atom2->getIndex()-1 << "-"
            //     << pTorsion->atom3->getIndex()-1 << "-" << pTorsion->atom4->getIndex()-1 << " : "
            //     << torsionParamList.size() << std::endl;
            n += torsionParamList.size();
          }
          //else {
            //std::cout << pTorsion->atom1->getIndex()-1 << "-" << pTorsion->atom2->getIndex()-1 << "-"
            //     << pTorsion->atom3->getIndex()-1 << "-" << pTorsion->atom4->getIndex()-1 << " : 0" << std::endl;
          //}
          //i++;
        }
      }
    }
    return n;
}

// ============================================================
// Function : getMMTorsions()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMMTorsions(int torsions[], double torsionParams[])
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getMMTorsions", " Can't find required parameters ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "molecule::getMMTorsions" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    int i = 0;
    int j = 0;
    try {
      if (!itsTorsionMap.empty()) {
        for (TorsionMapIterator b = itsTorsionMap.begin(); b != itsTorsionMap.end(); b++) {
          pTorsion = b->second;
          stdAtom* pStdAt1 = pTorsion->atom1->getStdAtom();
          stdAtom* pStdAt2 = pTorsion->atom2->getStdAtom();
          stdAtom* pStdAt3 = pTorsion->atom3->getStdAtom();
          stdAtom* pStdAt4 = pTorsion->atom4->getStdAtom();
          if (pStdAt1 and pStdAt2 and pStdAt3 and pStdAt4) {
            torsionParamList = pParameters->getTorsionParamList(
              pStdAt1->type, pStdAt2->type, pStdAt3->type, pStdAt4->type);
            if (!torsionParamList.empty()) {
              for (torsionParamIterator c = torsionParamList.begin(); c != torsionParamList.end(); c++) {
                pTorsionParam = *c;
                torsions[i  ] = pTorsion->atom1->getIndex()-1;
                torsions[i+1] = pTorsion->atom2->getIndex()-1;
                torsions[i+2] = pTorsion->atom3->getIndex()-1;
                torsions[i+3] = pTorsion->atom4->getIndex()-1;

                torsionParams[j]   = pTorsionParam->Vn / pTorsionParam->npth;
                //torsionParams[i+1] = pTorsionParam->npth;
                torsionParams[j+1] = pTorsionParam->Nt;
                torsionParams[j+2] = pTorsionParam->gamma;
                i += 4;
                j += 3;
              }
            }
            else {
              errorLogger.throwError("molecule::getMMTorsions", " torsionParamList is empty", INFO);
            }
            torsionParamList.clear();
          }
          else {
            errorLogger.throwError("molecule::getMMTorsions", " Can't find stdAtoms ", INFO);
          }
        }
      }
    }
    catch (std::exception& e) {
      std::string errMessage = e.what();
      errorLogger.throwError("molecule::getMMTorsions", errMessage, INFO);
    }

    return 0;
}

// ============================================================
// Function : getImproperMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<int, Improper*> molecule::getImproperMap()
{
    return itsImproperMap;
}

// ============================================================
// Function : getNumMMImpropers()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumMMImpropers()
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getNumMMImpropers", " Can't find required parameters ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "molecule::getNumMMImpropers" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    //int i = 0;
    int n = 0;
    std::vector<std::vector<int> > order;

    if (!itsImproperMap.empty()) {
      for (ImproperMapIterator b = itsImproperMap.begin(); b != itsImproperMap.end(); b++) {
        pImproper = b->second;
        improperParamList = pParameters->getImproperParamList(pImproper->atom1->getStdAtom()->type,
                            pImproper->atom2->getStdAtom()->type, pImproper->atom3->getStdAtom()->type,
                            pImproper->atom4->getStdAtom()->type, order);
        if (!improperParamList.empty()) {
          //std::cout << pImproper->atom1->getIndex()-1 << "-" << pImproper->atom2->getIndex()-1 << "-"
          //     << pImproper->atom3->getIndex()-1 << "-" << pImproper->atom4->getIndex()-1 << " : "
          //     << improperParamList.size() << std::endl;
          n+=improperParamList.size();
        }
        else {
          //std::cout << pImproper->atom1->getIndex()-1 << "-" << pImproper->atom2->getIndex()-1 << "-"
          //     << pImproper->atom3->getIndex()-1 << "-" << pImproper->atom4->getIndex()-1 << " : 0" << std::endl;
        }
        //i++;
      }
    }
    return n;
}

// ============================================================
// Function : getMMImpropers()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMMImpropers(int impropers[], double improperParams[])
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getMMImpropers", " Can't find required parameters ", MTK_ERROR);
      std::stringstream ss;
      ss << "molecule::getMMImpropers" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    int i = 0;
    int j = 0;
    int index = 0;
    std::vector<std::vector<int> > order;

    try {
      if (!itsImproperMap.empty()){
        for (ImproperMapIterator b = itsImproperMap.begin(); b != itsImproperMap.end(); b++) {
          pImproper = b->second;
          index = 0;
          improperParamList = pParameters->getImproperParamList(pImproper->atom1->getStdAtom()->type,
                                                                pImproper->atom2->getStdAtom()->type,
                                                                pImproper->atom3->getStdAtom()->type,
                                                                pImproper->atom4->getStdAtom()->type, order);
        //std::cout << " Improper2 " << pImproper->atom1->getIndex()-1 << "-" << pImproper->atom2->getIndex()-1  << "-"
        //                      << pImproper->atom3->getIndex()-1 << "-" << pImproper->atom4->getIndex()-1  << endl;

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
    catch (std::exception& e) {
      std::string errMessage = e.what();
      errorLogger.throwError("molecule::getMMImpropers", errMessage, MTK_ERROR);
    }
    return 0;
}

// ============================================================
// Function : getNumMMnonBondedPairs()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getNumMMnonBondedPairs(double cutoff)
{
    int nBPs = 0;
    unsigned int nSubMolecules = this->itsSubMoleculeList.size();
    double cutoff2 = cutoff*cutoff;
    double subMolCutoff = (cutoff*2)*(cutoff*2);

    // Compute all submolecule centroids
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      pSubMolecule->centerOfMass();
    }

    for (unsigned int i = 0; i < nSubMolecules; i++) {
      std::vector<atom*> atList1 = itsSubMoleculeList[i]->getAtomList();
      for (unsigned int k = 0; k < atList1.size(); k++) {
        for (unsigned int j = i; j < nSubMolecules; j++) {
          // INTRA RESIDUE
          if (j == i) {
            for (unsigned int k2 = k+1; k2 < atList1.size(); k2++) {
              if (atList1[k]->getCoords()->dist2( (*atList1[k2]->getCoords())) < cutoff2) {
                if ( (!atList1[k]->hasBondedAtom(atList1[k2])) && (!atList1[k]->has13BondedAtom(atList1[k2]))) {
                  nBPs++;
                }
              }
            }
          }
          // INTER RESIDUE
          else {
            std::vector<atom*> atList2 = itsSubMoleculeList[j]->getAtomList();
            if (itsSubMoleculeList[i]->getCenterOfMass()->dist2( (*itsSubMoleculeList[j]->getCenterOfMass())) < subMolCutoff) {
              for (unsigned int l = 0; l < atList2.size(); l++) {
                if (atList1[k]->getCoords()->dist2( (*atList2[l]->getCoords())) < cutoff2) {
                  if ( (!atList1[k]->hasBondedAtom(atList2[l])) && (!atList1[k]->has13BondedAtom(atList2[l]))) {
                    nBPs++;
                  }
                }
              }
            }
          }
        }
      }
    }
    //std::cout << " nBPs = " << nBPs << std::endl;
    return nBPs;
}

// ============================================================
// Function : getMMnonBondedPairs()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMMnonBondedPairs(int nonBonded[], int nonBondedPtrs[], double nonBondedParams[],
                                  int nonBonded14Ptrs[], double cutoff)
{
    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::getMMnonBondedPairs", " Can't find required parameters ", MTK_ERROR);
      std::stringstream ss;
      ss << "molecule::getMMnonBondedPairs" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    pParameters->calculateSigmaEpsilon();

    unsigned int nSubMolecules = this->itsSubMoleculeList.size();
    double cutoff2 = cutoff*cutoff;
    double subMolCutoff = (cutoff*2)*(cutoff*2);
    int index = 0;
    int indexParam = 0;
    int nNBPtr = 0;
    int n = 0;
    try {
      for (unsigned int i = 0; i < nSubMolecules; i++) {
        std::vector<atom*> atList1 = itsSubMoleculeList[i]->getAtomList();
        for (unsigned int k = 0; k < atList1.size(); k++) {
          nNBPtr = 0;
          for (unsigned int j = i; j < nSubMolecules; j++) {
            // INTRA RESIDUE
            if (j == i) {
              for (unsigned int k2 = k+1; k2 < atList1.size(); k2++) {
                if (atList1[k]->getCoords()->dist2( (*atList1[k2]->getCoords())) < cutoff2) {
                  if ( (!atList1[k]->hasBondedAtom(atList1[k2])) && (!atList1[k]->has13BondedAtom(atList1[k2]))) {
                    pLJ612SE = pParameters->getLJ612SE(atList1[k]->getStdAtom()->type, atList1[k2]->getStdAtom()->type);
                    if (pLJ612SE) {
                      nonBonded[index] = atList1[k2]->getIndex()-1;
                      nonBondedParams[indexParam] = pLJ612SE->sigma;
                      nonBondedParams[indexParam+1] = pLJ612SE->epsilon;
                      if (atList1[k]->has14BondedAtom(atList1[k2])) {
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
                    //std::cout << "      " << index << std::endl;
                  }
                }
              }
            }
            // INTER RESIDUE
            else {
              std::vector<atom*> atList2 = itsSubMoleculeList[j]->getAtomList();
              if (itsSubMoleculeList[i]->getCenterOfMass()->dist2( (*itsSubMoleculeList[j]->getCenterOfMass())) < subMolCutoff) {
                for (unsigned int l = 0; l < atList2.size(); l++) {
                  if (atList1[k]->getCoords()->dist2( (*atList2[l]->getCoords())) < cutoff2) {
                    if ( (!atList1[k]->hasBondedAtom(atList2[l])) && (!atList1[k]->has13BondedAtom(atList2[l]))) {
                      pLJ612SE = pParameters->getLJ612SE(atList1[k]->getStdAtom()->type, atList2[l]->getStdAtom()->type);
                      if (pLJ612SE) {
                        nonBonded[index] = atList2[l]->getIndex()-1;
                        if (atList1[k]->has14BondedAtom(atList2[l])) {
                          nonBonded14Ptrs[index] = 1;
                        }
                        else {
                          nonBonded14Ptrs[index] = 0;
                        }
                        nonBondedParams[indexParam] = pLJ612SE->sigma;
                        nonBondedParams[indexParam+1] = pLJ612SE->epsilon;
                        nNBPtr++;
                        index++;
                        indexParam+=2;
                      }
                      else {
                        nonBonded[index] = -1;
                        nonBonded14Ptrs[index] = 0;
                        nonBondedParams[indexParam] = 0;
                        nonBondedParams[indexParam+1] = 0;
                        nNBPtr++;
                        index++;
                        indexParam+=2;
                      }
                      //std::cout << "      " << index << std::endl;
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
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getMMnonBondedPairs", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}


// ============================================================
// Function : getMMCharges()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMMCharges(double charges[])
{
    unsigned int nSubMolecules = this->itsSubMoleculeList.size();
    int a = 0;
    try {
      for (unsigned int i = 0; i < nSubMolecules; i++) {
        std::vector<atom*> atList = itsSubMoleculeList[i]->getAtomList();
        for (unsigned int j = 0; j < atList.size(); j++) {
          stdAtom* pStdAt = atList[j]->getStdAtom();
          if (pStdAt) {
            charges[a] = pStdAt->atmCharge * E2KCAL;
          }
          else {
            charges[a] = 0.0;
          }
          a++;
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getMMCharges", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getMMResidues()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getMMResidues(int residues[])
{
    int nAtoms = 0;
    unsigned int nSubMolecules = this->itsSubMoleculeList.size();
    try {
      for (unsigned int i = 0; i < nSubMolecules; i++) {
        residues[i] = nAtoms;
        nAtoms += itsSubMoleculeList[i]->getNumAtoms();
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getMMResidues", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : determineRings()
// ------------------------------------------------------------
// Determine all rings in the molecule
// ============================================================
void molecule::determineRings()
{
    if (!pRings) {
      pRings = new rings(this);
    }
    pRings->determine();
    this->bRingsAssigned = true;

    if (this->inFileType != "pdb") {
      this->kekulizeRings();
    }
}

// ============================================================
// Function : kekulizeRings()
// ------------------------------------------------------------
//
// ============================================================
void molecule::kekulizeRings()
{
    if (!this->bRingsAssigned) {
      this->determineRings();
    }

    if (this->itsRings.size() > 0) {
      for (unsigned int i = 0; i < itsRings.size(); i++) {
        ring* pRing = itsRings[i];
        if (pRing) pRings->kekulize(pRing);
      }

#ifdef DEBUG
      std::string errMessage = "\nmolecule: SSSR(";
      for (unsigned int i = 0; i < itsRings.size(); i++) {
        errMessage += i2s(itsRings[i]->atoms.size());
        if (i != itsRings.size()-1 ) {
          errMessage + ",";
        }
      }
      errMessage += ") \n";

      for (unsigned int i = 0; i < itsRings.size(); i++) {
        errMessage += " ring [" + i2s(i+1) + "] :";
        for (unsigned int j = 0; j < itsRings[i]->atoms.size(); j++) {
          errMessage += " " + i2s(itsRings[i]->atoms[j]->getFileID());
        }
        errMessage += " planarity = " + i2s(itsRings[i]->planar) +
                      " aromaticity = " + i2s(itsRings[i]->aromatic) + " \n";
      }
      errMessage += " ###";
      errorLogger.throwError("molecule::kekulizeRings", errMessage, INFO);
#endif
    }
    this->bRingsKekulized = true;
}

// ============================================================
// Function : getRingInfo()
// ------------------------------------------------------------
//
// ============================================================
std::string molecule::getRingInfo()
{
    std::string ringInfo;
    unsigned int nRings = this->itsRings.size();

    ringInfo = "nRings(" + i2s(nRings) + ") \n";
    ringInfo += "SSSR(";
    for (unsigned int i = 0; i < nRings; i++) {
      ringInfo += i2s(this->itsRings[i]->atoms.size());
      if (i != nRings-1 ) {
        ringInfo += ",";
      }
    }
    ringInfo += ") \n";

    for (unsigned int i = 0; i < nRings; i++) {
      ringInfo += " ring [" + i2s(i+1) + "] :";
      for (unsigned int j = 0; j < itsRings[i]->atoms.size(); j++) {
        ringInfo += " " + i2s(itsRings[i]->atoms[j]->getFileID());
      }
      ringInfo += " planarity = " + i2s(itsRings[i]->planar) +
                  " aromaticity = " + i2s(itsRings[i]->aromatic) + " \n";
    }
    return ringInfo;
}

// ============================================================
// Function : addRing()
// ------------------------------------------------------------
// Add a ring to the molecule
// ============================================================
void molecule::addRing(std::vector<atom*> r)
{
    std::vector<atom*>::iterator atomIterator;
    bool unique = true;
    unsigned int c = 0;
    // Check to see if ring is already been added
    for (unsigned int i = 0; i < this->itsRings.size(); i++) {
      for (unsigned int j = 0; j < this->itsRings[i]->atoms.size(); j++) {
        atomIterator = std::find(r.begin(), r.end(), this->itsRings[i]->atoms[j]);
        if (atomIterator != r.end()) {
          c++;
        }
      }
      if (c == r.size()) {
        unique = false;
      }
      c = 0;
    }

    bool error = false;
    if (unique) {
      ring* myRing = new ring();
      Bond* pRingBond = 0;
      for (unsigned int i = 0; i < r.size(); i++) {
        if (i != r.size()-1) {
          pRingBond = this->getBond(r[i], r[i+1]);
          if (pRingBond) {
            pRingBond->topology = 1;
          }
          else {
            //std::string errMessage = " NO BOND BETWEEN " + r[i]->getName() + "-" + r[i+1]->getName();
            //errorLogger.throwError("molecule::addRing", errMessage, MTK_ERROR);
            error = true;
          }
        }
        else {
          pRingBond = this->getBond(r[i], r[0]);
          if (pRingBond) {
            pRingBond->topology = 1;
          }
          else {
            //std::string errMessage = " NO BOND BETWEEN " + r[i]->getName() + "-" + r[0]->getName();
            //errorLogger.throwError("molecule::addRing", errMessage, MTK_ERROR);
            error = true;
          }
        }
        myRing->atoms.push_back(r[i]);
      }
      myRing->size      = myRing->atoms.size();
      myRing->planar    = 0;
      myRing->aromatic  = 0;
      myRing->hetero    = 0;
      myRing->nHetero   = 0;
      myRing->nNitrogen = 0;
      myRing->nOxygen   = 0;
      myRing->nSulfur   = 0;

      if (error) {
        delete myRing;
      }
      else {

//#ifdef DEBUG
        std::string errMessage = "\n";

        Bond* pRingBond = 0;
        for (unsigned int i = 0; i < myRing->atoms.size(); i++) {
          if (i != myRing->atoms.size()-1) {
            pRingBond = this->getBond(myRing->atoms[i], myRing->atoms[i+1]);
            if (pRingBond) {
              errMessage += "    " + i2s(pRingBond->atom1->getFileID()) + "-" + i2s(pRingBond->atom2->getFileID())
                         + " type = " + i2s(pRingBond->type) + " topology = " + i2s(pRingBond->topology) + "\n";
            }
          }
          else {
            pRingBond = this->getBond(myRing->atoms[i], myRing->atoms[0]);
            if (pRingBond) {
              errMessage += "    " + i2s(pRingBond->atom1->getFileID()) + "-" + i2s(pRingBond->atom2->getFileID())
                         + " type = " + i2s(pRingBond->type) + " topology = " + i2s(pRingBond->topology) + "\n";
            }
          }
        }
        errorLogger.throwError("molecule::addRing", errMessage, INFO);
//#endif

        for (unsigned int i = 0; i < myRing->atoms.size(); i++) {
          std::string symbol = myRing->atoms[i]->getElement()->symbol;
          if (symbol == "O") {
            myRing->hetero = 1;
            myRing->nHetero++;
            myRing->nOxygen++;
          }
          if (symbol == "S") {
            myRing->hetero = 1;
            myRing->nHetero++;
            myRing->nSulfur++;
          }
          if (symbol == "N") {
            myRing->hetero = 1;
            myRing->nHetero++;
            myRing->nNitrogen++;
          }
        }
        this->itsRings.push_back(myRing);
        myRing->index = this->itsRings.size();
      }
    }
    //else {
    //  puts("ITS NOT UNIQUE");
    //}
}

// ============================================================
// Function : getRings()
// ------------------------------------------------------------
// Get rings
// ============================================================
std::vector<ring*> molecule::getRings()
{
    return this->itsRings;
}

// ============================================================
// Function : getLastAddedRing()
// ------------------------------------------------------------
// Get rings
// ============================================================
ring* molecule::getLastAddedRing()
{
    if (!this->itsRings.size() < 1) {
      return this->itsRings.back();
    }
    return 0;
}

// ============================================================
// Function : getRotatableBonds()
// ------------------------------------------------------------
// Finds all possible rotatable bonds in the molecule
// ============================================================
std::vector<Bond*> molecule::getRotatableBonds()
{
    errorLogger.throwError("molecule::getRotatableBonds", this->getName(), INFO);

    std::vector<atom*> atomList = this->getAtomList();

    std::vector<atom*>::iterator atomIterator1;
    std::vector<atom*>::iterator atomIterator2;

    bool bOk = true;
    this->itsRotatableBonds.clear();
    if (!itsBondMap.empty()) {
      for (BondMapIterator b = itsBondMap.begin(); b != itsBondMap.end(); b++) {
        pBond = b->second;

        atomIterator1 = std::find(atomList.begin(), atomList.end(), pBond->atom1);
        atomIterator2 = std::find(atomList.begin(), atomList.end(), pBond->atom2);

        // check if bond is intermolecular
        if ((atomIterator1 != atomList.end()) and
            (atomIterator2 != atomList.end())) {
          // single bond
          if (pBond->type != 1) {
            bOk = false;
          }
          if (pBond->topology == 1) { // ring
            bOk = false;
          }
          if (((pBond->atom1->bondedAtoms.size() == 1) or
               (pBond->atom2->bondedAtoms.size() == 1))) {
            bOk = false;
          }
          // Either atom is not terminal
          if (((pBond->atom1->getType() == 2) or
               (pBond->atom2->getType() == 2))) {
            bOk = false;
          }
/*
\todo figure out a better way to deal with this
          // Aromatic-Aromatic bond
          if (((pBond->atom1->getType() == 6) and
               (pBond->atom2->getType() == 6))) {
            bOk = false;
          }
*/
          if (bOk and pBond->rotatable) {
            this->itsRotatableBonds.push_back(pBond);
          }
          else {
            pBond->rotatable = 0;
          }
          bOk = true;
        }
      }
    }
    return this->itsRotatableBonds;
}

// ============================================================
// Function : setRotatableTorsions()
// ------------------------------------------------------------
// Set rotatable torsions
// ============================================================
void molecule::setRotatableTorsions(std::vector<Torsion*> tors)
{
    this->itsRotatableTorsions.clear();

    for (unsigned int i = 0; i < tors.size(); i++) {
      this->itsRotatableTorsions.push_back(tors[i]);
    }
}

// ============================================================
// Function : getRotatableTorsions()
// ------------------------------------------------------------
// Finds all possible rotatable torsions in the molecule
// ============================================================
std::vector<Torsion*> molecule::getRotatableTorsions()
{
    this->itsRotatableTorsions.clear();
    std::vector<Bond*> rotBonds = this->getRotatableBonds();

    std::vector<atom*> atomList = this->getAtomList();

    std::vector<atom*>::iterator atomIterator1;
    std::vector<atom*>::iterator atomIterator2;
    std::vector<atom*>::iterator atomIterator3;
    std::vector<atom*>::iterator atomIterator4;

    if (!itsTorsionMap.empty()) {
      for (unsigned int i = 0; i < rotBonds.size(); i++) {
        pBond = rotBonds[i];
        for (TorsionMapIterator b = itsTorsionMap.begin(); b != itsTorsionMap.end(); b++){
          pTorsion = b->second;
          if ( ((pTorsion->atom2->getFileID() == pBond->atom1->getFileID()) and
                (pTorsion->atom3->getFileID() == pBond->atom2->getFileID())) or
               ((pTorsion->atom2->getFileID() == pBond->atom2->getFileID()) and
                (pTorsion->atom3->getFileID() == pBond->atom1->getFileID())) ) {
            // check if torsion is intermolecular
            atomIterator1 = std::find(atomList.begin(), atomList.end(), pTorsion->atom1);
            atomIterator2 = std::find(atomList.begin(), atomList.end(), pTorsion->atom2);
            atomIterator3 = std::find(atomList.begin(), atomList.end(), pTorsion->atom3);
            atomIterator4 = std::find(atomList.begin(), atomList.end(), pTorsion->atom4);
            if ((atomIterator1 != atomList.end()) and (atomIterator2 != atomList.end()) and
                (atomIterator3 != atomList.end()) and (atomIterator4 != atomList.end())) {
              this->itsRotatableTorsions.push_back(pTorsion);
              break;
            }
          }
        }
      }
    }
    return this->itsRotatableTorsions;
}

// ============================================================
// Function : setupConformers()
// ------------------------------------------------------------
// Set up conformers
// ============================================================
std::vector<Torsion*> molecule::setupConformers()
{
    this->itsConformers.clear();
    pConformers = new conformers(this);

    std::vector<Torsion*> rotTors = this->getRotatableTorsions();

    if (rotTors.size() > 0) {
      for (unsigned int i = 0; i < rotTors.size(); i++) {
        pConformers->addRotatableTorsion(rotTors[i]->atom1->getIndex(),
                                         rotTors[i]->atom2->getIndex(),
                                         rotTors[i]->atom3->getIndex(),
                                         rotTors[i]->atom4->getIndex());
/*
#ifdef DEBUG
        std::cout << " torsion = "
                  << rotTors[i]->atom1->getFileID() << "-"
                  << rotTors[i]->atom1->getHybridization() << "-"
                  << rotTors[i]->atom1->getType() << " "
                  << rotTors[i]->atom2->getFileID() << "-"
                  << rotTors[i]->atom2->getHybridization() << "-"
                  << rotTors[i]->atom2->getType() << " "
                  << rotTors[i]->atom3->getFileID() << "-"
                  << rotTors[i]->atom3->getHybridization() << "-"
                  << rotTors[i]->atom3->getType() << " "
                  << rotTors[i]->atom4->getFileID() << "-"
                  << rotTors[i]->atom4->getHybridization() << "-"
                  << rotTors[i]->atom4->getType() << " "
                  << std::endl;
#endif
*/
      }
      pConformers->determineAtomList();
    }
    return rotTors;
}

// ============================================================
// Function : setupConformers()
// ------------------------------------------------------------
// Set up conformers
// ============================================================
std::vector<Torsion*> molecule::setupConformers(std::vector<atom*> frozenAtoms)
{
    std::vector<atom*>::iterator atomIterator1;
    std::vector<atom*>::iterator atomIterator2;
    std::vector<atom*>::iterator atomIterator3;
    std::vector<atom*>::iterator atomIterator4;

    this->itsConformers.clear();
    pConformers = new conformers(this);

    std::vector<Torsion*> rotTors = this->getRotatableTorsions();

    if (rotTors.size() > 0) {
      for (unsigned int i = 0; i < rotTors.size(); i++) {
        atomIterator1 = std::find(frozenAtoms.begin(), frozenAtoms.end(), rotTors[i]->atom1);
        atomIterator2 = std::find(frozenAtoms.begin(), frozenAtoms.end(), rotTors[i]->atom2);
        atomIterator3 = std::find(frozenAtoms.begin(), frozenAtoms.end(), rotTors[i]->atom3);
        atomIterator4 = std::find(frozenAtoms.begin(), frozenAtoms.end(), rotTors[i]->atom4);
/*
#ifdef DEBUG
        std::cout << " Check torsion = " << rotTors[i]->atom1->getFileID()
                  << " " << rotTors[i]->atom2->getFileID()
                  << " " << rotTors[i]->atom3->getFileID()
                  << " " << rotTors[i]->atom4->getFileID();
#endif
*/
        if ((atomIterator1 != frozenAtoms.end()) and
            (atomIterator2 != frozenAtoms.end()) and
            (atomIterator3 != frozenAtoms.end()) and
            (atomIterator4 != frozenAtoms.end()) ) {
/*
#ifdef DEBUG
        std::cout << " Rejected " << std::endl;
#endif
*/
          continue;
        }

        if (atomIterator4 != frozenAtoms.end()) {
          if ((atomIterator1 == frozenAtoms.end()) and
              (atomIterator2 == frozenAtoms.end()) ) {
            pConformers->addRotatableTorsion(rotTors[i]->atom4->getIndex(), rotTors[i]->atom3->getIndex(),
                                             rotTors[i]->atom2->getIndex(), rotTors[i]->atom1->getIndex());
/*
#ifdef DEBUG
        std::cout << " Accepted " << std::endl;
#endif
*/
            continue;
          }
        }

        pConformers->addRotatableTorsion(rotTors[i]->atom1->getIndex(), rotTors[i]->atom2->getIndex(),
                                         rotTors[i]->atom3->getIndex(), rotTors[i]->atom4->getIndex());
/*
#ifdef DEBUG
        std::cout << " Accepted " << std::endl;
#endif
*/
      }
      pConformers->determineAtomList();
    }
    return rotTors;
}

// ============================================================
// Function : setupConformers()
// ------------------------------------------------------------
// Set up conformers
// ============================================================
std::vector<Torsion*> molecule::setupConformers(std::vector<Bond*> rotBonds)
{
    errorLogger.throwError("molecule::setupConformers", this->getName(), INFO);
    this->itsRotatableTorsions.clear();
    this->itsConformers.clear();
    pConformers = new conformers(this);

    std::vector<atom*> atomList = this->getAtomList();

    std::vector<atom*>::iterator atomIterator1;
    std::vector<atom*>::iterator atomIterator2;
    std::vector<atom*>::iterator atomIterator3;
    std::vector<atom*>::iterator atomIterator4;

    if (!itsTorsionMap.empty()) {
      for (unsigned int i = 0; i < rotBonds.size(); i++) {
        pBond = rotBonds[i];
        for (TorsionMapIterator b = itsTorsionMap.begin(); b != itsTorsionMap.end(); b++){
          pTorsion = b->second;
          if ( ((pTorsion->atom2->getFileID() == pBond->atom1->getFileID()) and
                (pTorsion->atom3->getFileID() == pBond->atom2->getFileID())) or
               ((pTorsion->atom2->getFileID() == pBond->atom2->getFileID()) and
                (pTorsion->atom3->getFileID() == pBond->atom1->getFileID())) ) {
            // check if torsion is intermolecular
            atomIterator1 = std::find(atomList.begin(), atomList.end(), pTorsion->atom1);
            atomIterator2 = std::find(atomList.begin(), atomList.end(), pTorsion->atom2);
            atomIterator3 = std::find(atomList.begin(), atomList.end(), pTorsion->atom3);
            atomIterator4 = std::find(atomList.begin(), atomList.end(), pTorsion->atom4);
            if ((atomIterator1 != atomList.end()) and (atomIterator2 != atomList.end()) and
                (atomIterator3 != atomList.end()) and (atomIterator4 != atomList.end())) {
              this->itsRotatableTorsions.push_back(pTorsion);
              break;
            }
          }
        }
      }
    }

    if (this->itsRotatableTorsions.size() > 0) {
      for (unsigned int i = 0; i < this->itsRotatableTorsions.size(); i++) {
        pConformers->addRotatableTorsion(this->itsRotatableTorsions[i]->atom1->getIndex(),
                                         this->itsRotatableTorsions[i]->atom2->getIndex(),
                                         this->itsRotatableTorsions[i]->atom3->getIndex(),
                                         this->itsRotatableTorsions[i]->atom4->getIndex());
      }
      pConformers->determineAtomList();
    }
    return this->itsRotatableTorsions;
}

// ============================================================
// Function : setupConformers()
// ------------------------------------------------------------
// Set up conformers
// ============================================================
std::vector<Torsion*> molecule::setupConformers(std::vector<Torsion*> rotTors)
{
    this->itsConformers.clear();
    this->itsRotatableTorsions.clear();
    pConformers = new conformers(this);

    if (rotTors.size() > 0) {
      for (unsigned int i = 0; i < rotTors.size(); i++) {
        this->itsRotatableTorsions.push_back(rotTors[i]);
        pConformers->addRotatableTorsion(rotTors[i]->atom1->getIndex(), rotTors[i]->atom2->getIndex(),
                                         rotTors[i]->atom3->getIndex(), rotTors[i]->atom4->getIndex());
/*
#ifdef DEBUG
        std::cout << "  torsion = " << rotTors[i]->atom1->getFileID() << " " << rotTors[i]->atom2->getFileID() << " "
                                   << rotTors[i]->atom3->getFileID() << " " << rotTors[i]->atom4->getFileID() << std::endl;
#endif
*/
      }
      pConformers->determineAtomList();
    }
    return rotTors;
}


// ============================================================
// Function : addConformer()
// ------------------------------------------------------------
// Add a conformer
// ============================================================
void molecule::addConformer(std::string name, std::vector<double> t)
{
    if (!this->pConformers) {
      this->setupConformers();
    }

    conformer* myConformer = new conformer();
    for (unsigned int i = 0; i < t.size(); i++) {
      myConformer->torsions.push_back(t[i]);
    }
    myConformer->name = name;

    this->itsConformers.push_back(myConformer);
}

// ============================================================
// Function : updateConformer()
// ------------------------------------------------------------
// Update a conformer
// ============================================================
void molecule::updateConformer(std::string name, std::vector<double> t)
{
/*
#ifdef DEBUG
    std::cout << " molecule::updateConformer -> " << name << std::endl;
#endif
*/
    conformer* curConformer = 0;
    for (unsigned int i = 0; i < this->itsConformers.size(); i++) {
      if (this->itsConformers[i]->name == name) {
        curConformer = this->itsConformers[i];
      }
    }
    if (curConformer) {
      if (t.size() == curConformer->torsions.size()) {
        for (unsigned int i = 0; i < t.size(); i++) {
          curConformer->torsions[i] = t[i];
        }
      }
    }
}

// ============================================================
// Function : delConformer()
// ------------------------------------------------------------
// Delete a conformer
// ============================================================
void molecule::delConformer(std::string name)
{
    conformer* curConformer = 0;
    for (unsigned int i = 0; i < this->itsConformers.size(); i++) {
      if (this->itsConformers[i]->name == name) {
        curConformer = this->itsConformers[i];
      }
    }

    if (curConformer) {
      ConformerIterator c = std::find(itsConformers.begin(), itsConformers.end(), curConformer);

      if (c != itsConformers.end()) {
        itsConformers.erase(c);
      }
    }
}

// ============================================================
// Function : getConformer()
// ------------------------------------------------------------
//
// ============================================================
conformer* molecule::getConformer(std::string name)
{
    conformer* curConformer = 0;
    for (unsigned int i = 0; i < this->itsConformers.size(); i++) {
      curConformer = this->itsConformers[i];
      if (curConformer->name == name) {
        return curConformer;
      }
    }
    return 0;
}

// ============================================================
// Function : getConformers()
// ------------------------------------------------------------
// Get conformers pointer
// ============================================================
conformers* molecule::getConformers()
{
    if (this->pConformers) {
      return this->pConformers;
    }
    return 0;
}

// ============================================================
// Function : getConformerNames()
// ------------------------------------------------------------
// Get conformers pointer
// ============================================================
std::vector<std::string> molecule::getConformerNames()
{
    std::vector<std::string> conformerNames;
    for (unsigned int i = 0; i < this->itsConformers.size(); i++) {
      conformerNames.push_back(this->itsConformers[i]->name);
    }
    return conformerNames;
}

// ============================================================
// Function : generateSimpleFP()
// ------------------------------------------------------------
// Generate Simple Fingerprint
// ============================================================
void molecule::generateSimpleFP()
{
    if (!this->pFingerPrint) {
      this->pFingerPrint = new fingerPrint();
      this->pFingerPrint->generateSimpleFP(this, this->itsSimpleFP);
    }
    else {
      this->pFingerPrint->generateSimpleFP(this, this->itsSimpleFP);
    }
    delete pFingerPrint;
/*
#ifdef DEBUG
    for (unsigned int i = 0; i < this->itsSimpleFP.size(); i++) {
      std::cout << this->itsSimpleFP[i] << " ";
    }
    std::cout << " " << std::endl;
#endif
*/
}

// ============================================================
// Function : getSimpleFP()
// ------------------------------------------------------------
// Get Simple Fingerprint
// ============================================================
std::vector<unsigned int> molecule::getSimpleFP()
{
    return this->itsSimpleFP;
}

// ============================================================
// Function : generateFragmentFP()
// ------------------------------------------------------------
// Generate Fragment Fingerprint
// ============================================================
/*
void molecule::generateFragmentFP()
{
    if (!this->pFingerPrint) {
      this->pFingerPrint = new fingerPrint();
      this->pFingerPrint->generateFragmentFP(this, this->itsFragmentFP);
    }
    else {
      this->pFingerPrint->generateFragmentFP(this, this->itsFragmentFP);
    }
}
*/

// ============================================================
// Function : getFragmentFP()
// ------------------------------------------------------------
// Get Fragment Fingerprint
// ============================================================
std::string molecule::getFragmentFP()
{
    return this->itsFragmentFP;
}

// ============================================================
// Function : determineHydrophobicGroups()
// ------------------------------------------------------------
//
// ============================================================
int molecule::determineHydrophobicGroups()
{
    int failure = 0;
    if (!this->pHydrophobize) {
      this->pHydrophobize = new hydrophobize(this);
      failure = this->pHydrophobize->run();
    }
    else {
      failure = this->pHydrophobize->run();
    }
    delete this->pHydrophobize;
    return failure;
}

// ============================================================
// Function : addHydrophobicGroup()
// ------------------------------------------------------------
//
// ============================================================
int molecule::addHydrophobicGroup(std::vector<atom*> ats)
{
    hydrophobe* pHydrophobe = new hydrophobe();
    for (unsigned int i = 0; i < ats.size(); i++) {
      pHydrophobe->atoms.push_back(ats[i]);
    }
    this->itsHydrophobes.push_back(pHydrophobe);
    return 0;
}

// ============================================================
// Function : getHydrophobicGroups()
// ------------------------------------------------------------
//
// ============================================================
std::vector<hydrophobe*> molecule::getHydrophobicGroups()
{
    return this->itsHydrophobes;
}

// ============================================================
// Function : determineFunctionalGroups()
// ------------------------------------------------------------
//
// ============================================================
int molecule::determineFunctionalGroups(stdLibrary* pStdLib)
{
    int success = 0;
    if (!this->pFunctionalize) {
      this->pFunctionalize = new functionalize(this);
      success = this->pFunctionalize->run(pStdLib);
    }
    else {
      success = this->pFunctionalize->run(pStdLib);
    }
    delete this->pFunctionalize;
    if (!success) return 1;
    return 0;
}

// ============================================================
// Function : addFunctionalGroup()
// ------------------------------------------------------------
//
// ============================================================
void molecule::addFunctionalGroup(funcGroup* fg)
{
    errorLogger.throwError("molecule::addFunctionalGroup", fg->pStdFrag->getName(), INFO);

    this->itsFuncGroups.push_back(fg);
    this->itsFragmentFP += fg->pStdFrag->getCode();

    // Add molecule property (printed in sdf files)
    /*
       > <FuncGroup:6MR00BNZ:1>
       1;1|2;2|3;3|4;4|5;5|6;6|
    */
    int numFuncGroup = 0;
    for (unsigned int w = 0; w < this->itsFuncGroups.size(); w++) {
      if (fg->pStdFrag->getCode() == this->itsFuncGroups[w]->pStdFrag->getCode()) {
        numFuncGroup++;
      }
    }
    std::stringstream funcGroupNumber;
    funcGroupNumber << numFuncGroup;
    std::string fName = "FuncGroup:" + fg->pStdFrag->getCode()
                        + ":" + funcGroupNumber.str();
    std::string fValue = "";

    stdAtom* pStdAtom = 0;
    typedef std::map<stdAtom*, atom*>::iterator funcGroupMapIterator;
    for (funcGroupMapIterator f = fg->atomMap.begin();
         f != fg->atomMap.end(); f++) {
      pStdAtom = f->first;
      pAtom = f->second;
      std::stringstream stdAtomIndex;
      stdAtomIndex << pStdAtom->index;
      std::stringstream atomIndex;
      atomIndex << pAtom->getIndex();

      fValue += stdAtomIndex.str() + ";" +
                atomIndex.str()    + "|";
    }
    this->addProperty(fName, fValue);
}

// ============================================================
// Function : addFunctionalGroup()
// ------------------------------------------------------------
//
// ============================================================
void molecule::addFunctionalGroup(int fragAtoms, stdFrag* pStdFrag,
                                  int molAtoms, std::vector<int> subGraph)
{
    errorLogger.throwError("molecule::addFunctionalGroup", pStdFrag->getName(), INFO);

    if (!pParameters) {
      pParameters = pParent->getParameters();
    }
    if (!pParameters) {
      errorLogger.throwError("molecule::addFunctionalGroup", " Can't find required parameters ", MTK_ERROR);
      std::stringstream ss;
      ss << "molecule::addFunctionalGroup" << " Can't find required parameters ";
      throw MTKException(ss.str());
    }

    stdAtom* pStdAtom = 0;

    typedef std::map<stdAtom*, atom*>::iterator funcGroupAtomIterator;
    int numSame = 0;
    bool different = true;

    // Remove redundant functional groups
    for (unsigned int i = 0; i < this->itsFuncGroups.size(); i++) {
      if (this->itsFuncGroups[i]->pStdFrag->getName() == pStdFrag->getName()) {
        for (int t = 0; t < fragAtoms; t++) {
          for (int t2 = 0; t2 < molAtoms; t2++) {
            if (subGraph[t*molAtoms+t2]) {
              pAtom = this->getAtom(t2+1, 1, 0);
              if (pAtom) {
                if (!this->itsFuncGroups[i]->atomMap.empty()) {
                  for (funcGroupAtomIterator b = this->itsFuncGroups[i]->atomMap.begin();
                       b != this->itsFuncGroups[i]->atomMap.end(); b++) {
                    if (pAtom == b->second) {
                      numSame++;
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (numSame == fragAtoms) {
        //std::cout << "    different1 " << std::endl;
        different = false;
      }
      numSame = 0;
    }

    // Check to see if a larger graph related to this fragment was added
    for (unsigned int i = 0; i < this->itsFuncGroups.size(); i++) {
      std::vector<std::string> subGraphs = this->itsFuncGroups[i]->pStdFrag->getSubGraphs();
      for (unsigned int j = 0; j < subGraphs.size(); j++) {
        if (pStdFrag->getSymbol() == subGraphs[j]) {
          for (int t = 0; t < fragAtoms; t++) {
            for (int t2 = 0; t2 < molAtoms; t2++) {
              if (subGraph[t*molAtoms+t2]) {
                pAtom = this->getAtom(t2+1, 1, 0);
                if (pAtom) {
                  if (!this->itsFuncGroups[i]->atomMap.empty()) {
                    for (funcGroupAtomIterator b = this->itsFuncGroups[i]->atomMap.begin();
                         b != this->itsFuncGroups[i]->atomMap.end(); b++) {
                      if (pAtom == b->second) {
                        numSame++;
                      }
                    }
                  }
                }
              }
            }
          }
          if (numSame == fragAtoms) {
            //std::cout << "    different2 " << std::endl;
            different = false;
          }
          numSame = 0;
        }
      }
    }

    if (different) {
      funcGroup* myFuncGroup = new funcGroup();

      myFuncGroup->pStdFrag = pStdFrag;

      for (int t = 0; t < fragAtoms; t++) {
        for (int t2 = 0; t2 < molAtoms; t2++) {
          if (subGraph[t*molAtoms+t2]) {
            pAtom = this->getAtom(t2+1, 1, 0);
            pStdAtom = pStdFrag->getStdAtom(t+1);
            if (pAtom && pStdAtom) {
              myFuncGroup->atomMap[pStdAtom] = pAtom;
            }
          }
        }
      }

      std::vector<int> connPts = pStdFrag->getStdConnPtsList();
      for (unsigned int x = 0; x < connPts.size(); x++) {
        pStdAtom = pStdFrag->getStdAtom(connPts[x]);
        std::string atSymbol = pParameters->getAtomTypeSymbol(pStdAtom->type);
        if (pStdAtom && (atSymbol != "H")) {
          std::vector<stdAtom*> stdAtomBondedAtoms = pStdFrag->getBondedStdAtoms(pStdAtom);
          int numStdHeavy = 0;
          for (unsigned int y = 0; y < stdAtomBondedAtoms.size(); y++) {
            std::string bdAtSymbol = pParameters->getAtomTypeSymbol(stdAtomBondedAtoms[y]->type);
            if (bdAtSymbol != "H") numStdHeavy++;
          }
          numStdHeavy++;

          pAtom = myFuncGroup->atomMap[pStdAtom];
          std::vector<atom*> bondedAtoms = pAtom->getBondedAtoms();
          int numHeavy = 0;
          for (unsigned int y = 0; y < bondedAtoms.size(); y++) {
            if (bondedAtoms[y]->getElement()->symbol != "H") numHeavy++;
          }
          if (numHeavy < numStdHeavy) {
            //std::cout << "    different3 " << std::endl;
            different = false;
          }
        }
      }

      if (different) {
        //std::cout << "      Added " << std::endl;
        this->itsFuncGroups.push_back(myFuncGroup);
        this->itsFragmentFP += pStdFrag->getCode();

        // Add molecule property (printed in sdf files)
        /*
           > <FuncGroup:6MR00BNZ:1>
           1;1|2;2|3;3|4;4|5;5|6;6|
        */
        int numFuncGroup = 0;
        for (unsigned int w = 0; w < this->itsFuncGroups.size(); w++) {
          if (myFuncGroup->pStdFrag->getCode() == this->itsFuncGroups[w]->pStdFrag->getCode()) {
            numFuncGroup++;
          }
        }
        std::stringstream funcGroupNumber;
        funcGroupNumber << numFuncGroup;
        std::string fName = "FuncGroup:" + myFuncGroup->pStdFrag->getCode()
                            + ":" + funcGroupNumber.str();
        std::string fValue = "";

        typedef std::map<stdAtom*, atom*>::iterator funcGroupMapIterator;
        for (funcGroupMapIterator f = myFuncGroup->atomMap.begin();
             f != myFuncGroup->atomMap.end(); f++) {
          pStdAtom = f->first;
          pAtom = f->second;
          std::stringstream stdAtomIndex;
          stdAtomIndex << pStdAtom->index;
          std::stringstream atomIndex;
          atomIndex << pAtom->getIndex();

          fValue += stdAtomIndex.str() + ";" +
                    atomIndex.str()    + "|";
        }
        this->addProperty(fName, fValue);
      }
      else {
        delete myFuncGroup;
      }
    }
}

// ============================================================
// Function : delFunctionalGroup()
// ------------------------------------------------------------
//
// ============================================================
void molecule::delFunctionalGroup(funcGroup* fG, int i)
{
    errorLogger.throwError("molecule::delFunctionalGroup", fG->pStdFrag->getSymbol(), INFO);

    // Delete property
    std::stringstream funcGroupNumber;
    funcGroupNumber << i;
    std::string fName = "FuncGroup:" + fG->pStdFrag->getCode()
                            + ":" + funcGroupNumber.str();

    this->delProperty(fName);

    // Delete funcGroup from vector
    FuncGroupIterator f = std::find(this->itsFuncGroups.begin(), this->itsFuncGroups.end(), fG);
    if (f != this->itsFuncGroups.end()) {
       delete *f;
       this->itsFuncGroups.erase(f);
    }

    // Clear fragment fingerprint
    this->itsFragmentFP = "";
}

// ============================================================
// Function : getFunctionalGroups()
// ------------------------------------------------------------
//
// ============================================================
std::vector<funcGroup*> molecule::getFunctionalGroups()
{
    return this->itsFuncGroups;
}

// ============================================================
// Function : getFunctionalGroup()
// ------------------------------------------------------------
//
// ============================================================
funcGroup* molecule::getFunctionalGroup(atom* pAt)
{
    typedef std::map<stdAtom*, atom*>::iterator funcGroupAtomIterator;

    if (pAt) {
      for (unsigned int i = 0; i < this->itsFuncGroups.size(); i++) {
        if (!this->itsFuncGroups[i]->atomMap.empty()) {
          for (funcGroupAtomIterator b = this->itsFuncGroups[i]->atomMap.begin();
                  b != this->itsFuncGroups[i]->atomMap.end(); b++) {
            if (pAt == b->second) {
              return this->itsFuncGroups[i];
            }
          }
        }
      }
    }
    return 0;
}

// ============================================================
// Function : generateAdjMatrix()
// ------------------------------------------------------------
//
// ============================================================
int molecule::generateAdjMatrix()
{
    std::vector<atom*> atomList = this->getAtomList();
    int nAtoms = atomList.size();

    this->adjMatrixSize = nAtoms*nAtoms;
    try {
      // can I delete adjMatrix before allocating it.
      this->adjMatrix = 0;
      this->adjMatrix   = new int [nAtoms*nAtoms];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::generateAdjMatrix", " Memory Allocation Failure", MTK_ERROR);
      return 1;
    }
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        pBond = this->getBond(atomList[i], atomList[j]);
        if (pBond) {
          this->adjMatrix[i*nAtoms+j] = pBond->type;
        }
        else {
          this->adjMatrix[i*nAtoms+j] = 0;
        }
      }
    }
/*
#ifdef DEBUG
    std::cout << "   molecule::generateAdjMatrix " << std::endl;
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        std::cout << this->adjMatrix[i*nAtoms+j] << " ";
      }
      std::cout << " " << std::endl;
    }
    std::cout << " " << std::endl;
#endif
*/
    return 0;
}

// ============================================================
// Function : generateHeavyAdjMatrix()
// ------------------------------------------------------------
//
// ============================================================
int molecule::generateHeavyAdjMatrix()
{
    std::vector<atom*> heavyAtomList = this->getHeavyAtomList();
    int nHeavyAtoms = heavyAtomList.size();

    this->heavyAdjMatrixSize = nHeavyAtoms*nHeavyAtoms;
    try {
      // can I delete adjMatrix before allocating it.
      this->heavyAdjMatrix = 0;
      this->heavyAdjMatrix   = new int [this->heavyAdjMatrixSize];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::generateHeavyAdjMatrix", " Memory Allocation Failure", MTK_ERROR);
      return 1;
    }

    for (int i = 0; i < nHeavyAtoms; i++) {
      for (int j = 0; j < nHeavyAtoms; j++) {
        pBond = this->getBond(heavyAtomList[i], heavyAtomList[j]);
        if (pBond) {
          this->heavyAdjMatrix[i*nHeavyAtoms+j] = pBond->type;
        }
        else {
          this->heavyAdjMatrix[i*nHeavyAtoms+j] = 0;
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getAdjMatrix()
// ------------------------------------------------------------
//
// ============================================================
int* molecule::getAdjMatrix()
{
    return this->adjMatrix;
}

// ============================================================
// Function : getHeavyAdjMatrix()
// ------------------------------------------------------------
//
// ============================================================
int* molecule::getHeavyAdjMatrix()
{
    return this->heavyAdjMatrix;
}

// ============================================================
// Function : getAdjMatrixSize()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getAdjMatrixSize()
{
    return this->adjMatrixSize;
}

// ============================================================
// Function : getHeavyAdjMatrixSize()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getHeavyAdjMatrixSize()
{
    return this->heavyAdjMatrixSize;
}

// ============================================================
// Function : generateFeatureDistMatrix()
// ------------------------------------------------------------
//
// ============================================================
int molecule::generateFeatureDistMatrix()
{
    std::vector< vector3d > coords;
    this->getCoordinates(coords);
    int r = this->generateFeatureDistMatrix(coords);
    return r;
}

// ============================================================
// Function : generateFeatureDistMatrix()
// ------------------------------------------------------------
//
// ============================================================
int molecule::generateFeatureDistMatrix(std::vector< vector3d > coords)
{
    errorLogger.throwError("molecule::generateFeatureDistMatrix", "", INFO);

    this->featureGroups.erase(this->featureGroups.begin(), this->featureGroups.end());
    this->featureLabels.erase(this->featureLabels.begin(), this->featureLabels.end());
    this->featureCoordinates.erase(this->featureCoordinates.begin(), this->featureCoordinates.end());

    bool oC = false;
    if (coords.size() != 0) oC = true;
    stdAtom* pStdAtom = 0;
    int nFeatures = 0;

    // Get coordinates for each feature
    for (unsigned int i = 0; i < this->itsFuncGroups.size(); i++) {
      std::vector<stdFeature*> featureList = this->itsFuncGroups[i]->pStdFrag->getStdFeatureList();
      for (unsigned int j = 0; j < featureList.size(); j++) {
        this->featureGroups.push_back(this->itsFuncGroups[i]->pStdFrag->getSymbol());
        this->featureLabels.push_back(featureList[j]->name);
        nFeatures++;
        if (featureList[j]->atoms.size() > 1) {
          // Calculate the geometric centroid of the atoms
          vector3d Coord;
          vector3d Center;
          for (unsigned int k = 0; k < featureList[j]->atoms.size(); k++) {
            pStdAtom = this->itsFuncGroups[i]->pStdFrag->getStdAtom(featureList[j]->atoms[k]);
            pAtom = this->itsFuncGroups[i]->atomMap[pStdAtom];
            if (oC) {
              Coord = coords[pAtom->getIndex()-1];
            }
            else {
              Coord = (*pAtom->getCoords());
            }
            Center = Center+Coord;
          }
          Center = Center / featureList[j]->atoms.size();
          this->featureCoordinates.push_back(Center);
        }
        else {
          pStdAtom = this->itsFuncGroups[i]->pStdFrag->getStdAtom(featureList[j]->atoms[0]);
          if (!pStdAtom) {
            errorLogger.throwError("molecule::generateFeatureDistMatrix", " Error1", MTK_ERROR);
            std::stringstream ss;
            ss << "molecule::generateFeatureDistMatrix" << " Error1";
            throw MTKException(ss.str());
          }

          pAtom = this->itsFuncGroups[i]->atomMap[pStdAtom];
          if (!pAtom) {
            errorLogger.throwError("molecule::generateFeatureDistMatrix", " Error2", MTK_ERROR);
            std::stringstream ss;
            ss << "molecule::generateFeatureDistMatrix"<< " Error2";
            throw MTKException(ss.str());
          }
          vector3d Coord;
          if (oC) {
            Coord = coords[pAtom->getIndex()-1];
          }
          else {
            Coord = (*pAtom->getCoords());
          }
          this->featureCoordinates.push_back(Coord);
        }
      }
    }

    // Loop over all rings in the molecule
    for (unsigned int i = 0; i < this->itsRings.size(); i++) {
      nFeatures+=3;
      pRing = this->itsRings[i];

      if (pRings) {
        int ringNormalResult = pRings->getPlaneNormal(pRing);
        if (ringNormalResult) {
          errorLogger.throwError("molecule::generateFeatureDistMatrix", " Error", MTK_ERROR);
          std::stringstream ss;
          ss << "molecule::generateFeatureDistMatrix" << " ring normal result error";
          throw MTKException(ss.str());
        }
      }

      if (pRing->aromatic) {
        this->featureLabels.push_back("PIC");
      }
      else {
        this->featureLabels.push_back("RNG");
      }
      this->featureGroups.push_back("RING");

      vector3d Coord;
      pRings->calcCentroid(pRing);
      Coord[0] = pRing->centroid[0];
      Coord[1] = pRing->centroid[1];
      Coord[2] = pRing->centroid[2];
      this->featureCoordinates.push_back(Coord);

      vector3d Coord2;
      for (unsigned int j = 0; j < 3; j++) {
        Coord2[j] = pRing->planeNormal(j,2) + pRing->centroid(j);
      }
      this->featureGroups.push_back("RING");
      this->featureLabels.push_back("RNL");
      this->featureCoordinates.push_back(Coord2);

      vector3d Coord3;
      for (unsigned int j = 0; j < 3; j++) {
        Coord3[j] = - pRing->planeNormal(j,2) + pRing->centroid(j);
      }
      this->featureGroups.push_back("RING");
      this->featureLabels.push_back("RNL");
      this->featureCoordinates.push_back(Coord3);
    }

    // Loop over hydrophobes
    for (unsigned int i = 0; i < this->itsHydrophobes.size(); i++) {
      nFeatures++;
      hydrophobe* pHydrophobe = this->itsHydrophobes[i];
      this->featureLabels.push_back("HPB");
      this->featureGroups.push_back("HYDROPHOBE");

      vector3d Coord;
      vector3d Center;
      for (unsigned int j = 0; j < pHydrophobe->atoms.size(); j++) {
        pAtom = pHydrophobe->atoms[j];
        Coord = (*pAtom->getCoords());
        Center = Center+Coord;
      }
      Center = Center / double(pHydrophobe->atoms.size());
      this->featureCoordinates.push_back(Center);
    }

    if (nFeatures == 0) {
      std::string errMessage = " No features found in " + this->getName();
      errorLogger.throwError("molecule::generateFeatureDistMatrix", errMessage, MTK_ERROR);
      return 1;
    }

    this->featureDistMatrixSize = nFeatures*nFeatures;

    try {
      // can I delete matrix before allocating it.
      delete [] this->featureDistMatrix;
      this->featureDistMatrix = 0;
      this->featureDistMatrix = new double [nFeatures*nFeatures];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::generateFeatureDistMatrix", " Memory Allocation Failure", MTK_ERROR);
      return 1;
    }

    for (int i = 0; i < nFeatures; i++) {
      for (int j = i; j < nFeatures; j++) {
        double distance = featureCoordinates[i].dist(featureCoordinates[j]);
        if (i != j) {
          this->featureDistMatrix[i*nFeatures+j] = distance;
          this->featureDistMatrix[j*nFeatures+i] = distance;
        }
        else {
          this->featureDistMatrix[i*nFeatures+j] = 0;
        }
      }
    }

#ifdef DEBUG
    std::string errMessage = "    Feature Distance Matrix: \n";
    for (int i = 0; i < nFeatures; i++) {
      for (int j = 0; j < nFeatures; j++) {
        errMessage += d2s(this->featureDistMatrix[i*nFeatures+j]) + " ";
      }
      errMessage += " \n";
    }

    errMessage += "    Feature Coordinates: \n";
    for (unsigned int k = 0 ; k < featureLabels.size(); k++) {
      errMessage += featureLabels[k] + " " + d2s(featureCoordinates[k].getX()) + " "
                 + d2s(featureCoordinates[k].getY()) + " "
                 + d2s(featureCoordinates[k].getZ()) + "\n";
    }

    errMessage += " \n";
    errorLogger.throwError("molecule::generateFeatureDistMatrix", errMessage, INFO);
#endif
    return 0;
}

// ============================================================
// Function : getFeatureDistMatrix()
// ------------------------------------------------------------
//
// ============================================================
double* molecule::getFeatureDistMatrix()
{
    return this->featureDistMatrix;
}

// ============================================================
// Function : getFeatureGroups()
// ------------------------------------------------------------
//
// ============================================================
std::vector<std::string> molecule::getFeatureGroups()
{
    return this->featureGroups;
}

// ============================================================
// Function : getFeatureLabels()
// ------------------------------------------------------------
//
// ============================================================
std::vector<std::string> molecule::getFeatureLabels()
{
    return this->featureLabels;
}

// ============================================================
// Function : getFeatureCoordinates()
// ------------------------------------------------------------
//
// ============================================================
std::vector<vector3d> molecule::getFeatureCoordinates()
{
    return this->featureCoordinates;
}

// ============================================================
// Function : findPharmacophore()
// ------------------------------------------------------------
//
// ============================================================
int molecule::findPharmacophore(molecule* pMol,
              std::vector<std::vector<unsigned int> > &mcp,
              std::vector<vector3d> &mcpCoords,
              double distMax)
{
    //double hbda, double picr, double pncc, double hlbc
    //this->pPharmacophore = new pharmacophore(this, distMax, hbda, picr, pncc, hlbc);

    int failure = 0;
    if (!this->pPharmacophore) {
      this->pPharmacophore = new pharmacophore(this, distMax);
      failure = this->pPharmacophore->run(pMol, mcp, mcpCoords);
    }
    else {
      failure = this->pPharmacophore->run(pMol, mcp, mcpCoords);
    }
    delete this->pPharmacophore;
    //this->pPharmacophore = 0;
    return failure;
}

// ============================================================
// Function : findPharmacophore()
// ------------------------------------------------------------
//
// ============================================================
int molecule::findPharmacophore(molecule* pMol,
              std::vector<clique*> &cliqueList,
              double distMax)//, double hbda,
//              double picr, double pncc, double hlbc)
{
#ifdef DEBUG
      //std::cout << "   molecule::findPharmacophore" << std::endl;
#endif
    //this->pPharmacophore = new pharmacophore(this, distMax, hbda, picr, pncc, hlbc);

    int failure = 0;
    if (!this->pPharmacophore) {
      this->pPharmacophore = new pharmacophore(this, distMax);
      failure = this->pPharmacophore->getCliques(pMol, cliqueList);
    }
    else {
      failure = this->pPharmacophore->getCliques(pMol, cliqueList);
    }
    delete this->pPharmacophore;
    this->pPharmacophore = 0;
    return failure;
}

// ============================================================
// Function : addProperty()
// ------------------------------------------------------------
// Add molecular property
// ============================================================
void molecule::addProperty(const std::string &name, const std::string& value)
{
    itsPropertiesMap[name] = value;
    //itsPropertiesMap.insert( make_pair(name, value) );
}

// ============================================================
// Function : getProperty()
// ------------------------------------------------------------
// Get molecular property
// ============================================================
std::string molecule::getProperty(const std::string &name)
{
    PropertyMapIterator p = itsPropertiesMap.find(name);

    if (p != itsPropertiesMap.end()) {
      return p->second;
    }
    return "";
}

// ============================================================
// Function : delProperty()
// ------------------------------------------------------------
// Delete molecular property
// ============================================================
void molecule::delProperty(const std::string &name)
{
    PropertyMapIterator p = itsPropertiesMap.find(name);

    if (p != itsPropertiesMap.end()) {
      this->itsPropertiesMap.erase(p);
    }
}


// ============================================================
// Function : getProperties()
// ------------------------------------------------------------
// Get all molecular properties
// ============================================================
std::map<std::string, std::string> molecule::getProperties()
{
    return this->itsPropertiesMap;
}

// ============================================================
// Function : hasProperty()
// ------------------------------------------------------------
// Has molecule of a certain property
// ============================================================
bool molecule::hasProperty(const std::string &name)
{
    PropertyMapIterator p1 = this->itsPropertiesMap.find(name);
    if (p1 != this->itsPropertiesMap.end()) {
      return true;
    }
    return false;
}

// ============================================================
// Function : setInSolution()
// ------------------------------------------------------------
// Set the molecule in solution
// ============================================================
void molecule::setInSolution()
{
    std::string errorMessage = " Mol: " +  this->getName()+ "\n";

    if (!bAnglesAssigned) {
      errorMessage += " Assign angles before calling this function";

      errorLogger.throwError("molecule::setInSolution", errorMessage, MTK_ERROR);
      std::stringstream ss;
      ss << "molecule::setInSolution" << errorMessage;
      throw MTKException(ss.str());
    }

    atom* pAtom2 = 0;
    Bond* pBd = 0;
    Bond* pBd1 = 0;
    Bond* pBd2 = 0;
    std::string symbol = "";
    std::vector<atom*> atomList = this->getAtomList();
    std::vector<atom*> bondedAtoms;
    std::vector<atom*> bondedAtoms2;

    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      int nBondedAtoms = pAtom->getNumBonds();
      int g = pAtom->getElement()->group;
      if (nBondedAtoms == 0) {
        if (g == 1) { // 1. Z_i = {Group 1}, Q_i = 0 --> +1
          pAtom->setFormalCharge(+1);
        }
        else if (g == 2) { // 2. Z_i = {Group 2}, Q_i = 0 --> +2
          pAtom->setFormalCharge(+2);
        }
        else if (g == 17) { // 3. Z_i = {Group 7}, Q_i = 0 --> -1
          pAtom->setFormalCharge(-1);
        }
      }
    }

    // 2. Hydroxamic acid
    /*
          HO       O
            \    //
            HN--C
                 \
                  R
    */
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      int atType = pAtom->getType();

      // 5 = Single or Double Bond
      // 6 = Single or Aromatic Bond
      int bAtType = ((atType != 5) and (atType != 6));
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "N") and (bondedAtoms.size() == 2) and bAtType) {
        if ( ((bondedAtoms[0]->getElement()->symbol == "C") and
              (bondedAtoms[1]->getElement()->symbol == "O")) or
             ((bondedAtoms[1]->getElement()->symbol == "C") and
              (bondedAtoms[0]->getElement()->symbol == "O"))) {
          atom* pAtom_C = 0;
          atom* pAtom_O = 0;
          atom* pAtom_O2 = 0;
          bool bC_OK = false;
          bool bO_OK = false;

          for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
            if (bondedAtoms[j]->getElement()->symbol == "C") {
              pAtom_C = bondedAtoms[j];
              bondedAtoms2 = bondedAtoms[j]->getBondedAtoms();
              for (unsigned int k = 0; k < bondedAtoms2.size(); k++) {
                if ((bondedAtoms2[k] != pAtom) and (bondedAtoms2[k]->getElement()->symbol == "O")) {
                  if (bondedAtoms2[k]->getType() == 2) {
                    bC_OK = true;
                    pAtom_O2 = bondedAtoms2[k];
                  }
                }
              }
            }
            else if (bondedAtoms[j]->getElement()->symbol == "O")  {
              if (bondedAtoms[j]->getType() == 2) {
                pAtom_O = bondedAtoms[j];
                bO_OK = true;
              }
            }
          }

          if (bC_OK and bO_OK) {
            Bond* pBd_N_O = this->getBond(pAtom, pAtom_O);
            Bond* pBd_N_C = this->getBond(pAtom, pAtom_C);
            Bond* pBd_C_O = this->getBond(pAtom_C, pAtom_O2);

            if (pBd_N_O and pBd_N_C and pBd_C_O) {
              if ((pBd_N_O->type != 1) or (pBd_N_C->type != 1) or (pBd_C_O->type != 2)) {
                pBd_N_O->type = 1;
                pBd_N_C->type = 1;
                pBd_C_O->type = 2;
              }
              errorMessage += " Found hydroxamic acid. \n";
              pAtom->setFormalCharge(-1);
            }
          }
        }
      }
    }

    // 3. Chain amide bond
    /*
          H      O
           \   //
            N-C
          /    \
        R       R
    */
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      int atType = pAtom->getType();
      int bAtType = ((atType != 5) and (atType != 6));
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "N") and (bondedAtoms.size() == 2) and bAtType) {
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          if (bondedAtoms[j]->getElement()->symbol == "C") {
            bondedAtoms2 = bondedAtoms[j]->getBondedAtoms();
            for (unsigned int k = 0; k < bondedAtoms2.size(); k++) {
              if ((bondedAtoms2[k] != pAtom) and (bondedAtoms2[k]->getElement()->symbol == "O")) {
                if (bondedAtoms2[k]->getType() == 2) {
                  pBd1 = this->getBond(pAtom, bondedAtoms[j]);
                  pBd2 = this->getBond(bondedAtoms[j], bondedAtoms2[k]);
                  if (pBd1 and pBd2) {
                    if ((pBd1->type == 2) and (pBd2->type == 1)) {
                      pBd1->type = 1;
                      pBd2->type = 2;
                      errorMessage += " Found non-chain amide. Resetting the bond orders. \n";
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // 4. Alphiatic amine nitrogens (sp3 not bonded to an atom with a higher order bond) --> +1
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "N") and (bondedAtoms.size() == 4)) {
        bool all = true;
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          bondedAtoms2 = bondedAtoms[j]->getBondedAtoms();
          for (unsigned int k = 0; k < bondedAtoms2.size(); k++) {
            if (bondedAtoms2[k] != pAtom) {
              pBd = this->getBond(bondedAtoms[j], bondedAtoms2[k]);
              if (pBd->type > 1) {
                all = false;
                break;
              }
            }
          }
        }
        if (all) {
          pAtom->setFormalCharge(+1);
          errorMessage += " Found alphiatic amine nitrogen. Setting charge to +1. \n";
        }
      }
    }

    // 5a. carboxylic acids (-1)
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "O") and (bondedAtoms.size() == 1)) {
        if (bondedAtoms[0]->getElement()->symbol == "C") {
          for (unsigned int j = 0; j < atomList.size(); j++) {
            pAtom2 = atomList[j];
            if (pAtom2->getElement()->symbol == "O") {
              if (pAtom->has13BondedAtom(pAtom2)) {
                pBond = this->getBond(bondedAtoms[0], pAtom2);
                if (pBond->type == 2) {
                  pAtom->setFormalCharge(-1);
                  errorMessage += " Found carboxylic acid. Setting charge to -1. \n";
                }
              }
            }
          }
        }
      }
    }

    // 5b. sulfonic acid (-1)
    /*
           O
           \\
            S==0
          / \
        R    O-
    */
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "S") and (bondedAtoms.size() == 4)) {
        int nOxygens = 0;
        for (unsigned h = 0; h < bondedAtoms.size(); h++) {
          if ((bondedAtoms[h]->getElement()->symbol == "O") and
              (bondedAtoms[h]->getType() == 2)) {
            nOxygens++;
          }
        }
        if (nOxygens == 3) {
          for (unsigned h = 0; h < bondedAtoms.size(); h++) {
            pBond = this->getBond(bondedAtoms[h], pAtom2);
            if ((bondedAtoms[h]->getElement()->symbol == "O") and
                (pBond->type == 1)) {
              bondedAtoms[h]->setFormalCharge(-1);
              errorMessage += " Found sulfonic acid. Setting charge to -1. \n";
              break;
            }
          }
        }
      }
    }

    // 5c. phosphonic acid (-1)
    /*
           O   O-
           \\ /
            P
          /  \
        R     O-R
    */
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "O") and (bondedAtoms.size() == 1)) {
        if (bondedAtoms[0]->getElement()->symbol == "P") {
          //Bond* pB1 = this->getBond(bondedAtoms[0], pAtom);
          for (unsigned int j = 0; j < atomList.size(); j++) {
            pAtom2 = atomList[j];
            if (pAtom2->getElement()->symbol == "O") {
              if (pAtom->has13BondedAtom(pAtom2)) {
                pBond = this->getBond(bondedAtoms[0], pAtom2);
                if (pBond->type == 2) {
                  pAtom->setFormalCharge(-1);
                  errorMessage += " Found phosphonic acid. Setting charge to -1. \n";
                  break;
                }
              }
            }
          }
        }
      }
    }

    // 6. amidine      guanidine (+1)
    /*
              +NH2          NH2
             //            /
        R---C        +HN==C
             \        /    \
              NH2    R      NH2
    */
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      bondedAtoms = pAtom->getBondedAtoms();
      pAtom2 = 0;
      if ((pAtom->getElement()->symbol == "N") and
          (pAtom->getType() == 2 or pAtom->getType() == 3)) { // terminal or open chain
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          if (bondedAtoms[j]->getElement()->symbol == "C") {
            pBd = this->getBond(bondedAtoms[j], pAtom);
            if (pBd->type == 2) {
              pAtom2 = bondedAtoms[j];
            }
          }
        }
        if (pAtom2) {
          bondedAtoms2 = pAtom2->getBondedAtoms();
          for (unsigned int k = 0; k < bondedAtoms2.size(); k++) {
            if (bondedAtoms2[k] != pAtom) {
              int at2Type = bondedAtoms2[k]->getType();
              if ((bondedAtoms2[k]->getElement()->symbol == "N") and
                  (at2Type < 5 or at2Type > 6)) {
                pBd = this->getBond(bondedAtoms2[k], pAtom2);
                if (pBd->type == 1) {
                  pAtom->setFormalCharge(1);
                  bondedAtoms[k]->setHybridization(3); // sp2
                  errorMessage += " Found amidine/guanidine. \n";
                }
              }
            }
          }
        }
      }
    }

    // 7. sulfonamide (-1)
    /*
           O     O
           \\  //
        HN-- S
         -    \
              R
    */
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      bondedAtoms = pAtom->getBondedAtoms();
      if ((pAtom->getElement()->symbol == "N") and (bondedAtoms.size() == 1)) {
        if (bondedAtoms[0]->getElement()->symbol == "S") {
          atom* pSAtom = bondedAtoms[0];
          std::vector<atom*> SbondedAtoms = pSAtom->getBondedAtoms();
          int nOxygens = 0;
          for (unsigned h = 0; h < SbondedAtoms.size(); h++) {
            if (SbondedAtoms[h] == pAtom) continue;
            pBond = this->getBond(SbondedAtoms[h], pSAtom);
            if ((SbondedAtoms[h]->getElement()->symbol == "O") and
                (pBond->type == 2)) {
              nOxygens++;
            }
          }
          if (nOxygens == 2) {
            errorMessage += " Found sulfonamide. \n";
            pAtom->setFormalCharge(-1);
            break;
          }
        }
      }
    }

    errorLogger.throwError("molecule::setInSolution", errorMessage, WARNING);

    this->bInSolution = true;
}

// ============================================================
// Function : determineValences()
// ------------------------------------------------------------
// Set all valences
// ============================================================
void molecule::determineValences()
{
    double v = 0;
    atom* pAtom2 = 0;
    std::string symbol = "";
    std::vector<atom*> atomList = this->getAtomList();
    std::vector<atom*> bondedAtoms;

//std::cout << atomList.size() << std::endl;

    std::string errMessage = " Mol: " + this->getName() + "\n";
    // Determine valence for each atom in the molecule
    for (unsigned int i = 0; i < atomList.size(); i++) {

      pAtom = atomList[i];
      if (!pAtom) {
        errMessage = " Mol: " + this->getName() + " error finding atom";
        errorLogger.throwError("molecule::determineValences", errMessage, MTK_ERROR);
        return;
      }

      if (pAtom->getStdAtom()) continue;

      if (!pAtom->getElement()) {
        errMessage = " Mol: " + this->getName() + " error finding element";
        errorLogger.throwError("molecule::determineValences", errMessage, MTK_ERROR);
        return;
      }

      v = pAtom->getElement()->valence;

//std::cout << pAtom->getIndex() << ":" << pAtom->getElement()->symbol << std::endl;

/*
      errMessage += i2s(pAtom->getElement()->valence) + " ";
*/
      bondedAtoms = pAtom->getBondedAtoms();
      symbol = pAtom->getElement()->symbol;
      for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
        pAtom2 = bondedAtoms[j];
        if (!pAtom2) {
          errMessage = " Mol: " + this->getName() + " error finding atom";
          errorLogger.throwError("molecule::determineValences", errMessage, MTK_ERROR);
          return;
        }

// std::cout <<"     " << pAtom2->getIndex() << ":" << pAtom2->getElement()->symbol << std::endl;

        pBond = this->getBond(pAtom,pAtom2);
        if (!pBond) {
          errorLogger.throwError("molecule::determineValences",
           " Intermolecule bond found assuming its single ", INFO);
          v += 1;
          continue;
        }

        if (pBond->type == 0) {
          errorLogger.throwError("molecule::determineValences",
           " A bond type is undefined ", INFO);
          v += 1;
        }
        else if (pBond->type < 4) { // single, double, triple
          v += pBond->type;
/*
          errMessage += i2s(pBond->type) + " ";
*/
        }
        else if (pBond->type == 4) { // aromatic
          if (symbol == "C") {
            v += 1.5;
          }
          else {
            v += 1.5;
          }
        }
        else if (pBond->type == 5) { // single or double bond
          v += 2;
        }
        else if (pBond->type == 6) { // single or aromatic
          v += 1;
        }
        else if (pBond->type == 7) { // double or aromatic
          v += 2;
        }
      }

      int fc = pAtom->getFormalCharge();
      for (int j = 0; j < abs(fc); j++) {
        if (fc < 0) {
          v++;
        }
        else {
          v--;
        }
      }

      /// Account for neighbors formal charge, check this code: (pAtom->getFormalCharge() == 0)
      if ((pAtom->getType() == 2) and (pAtom->getFormalCharge() == 0)) {
        fc = pAtom2->getFormalCharge();
        if (pAtom->getElement()->symbol != "O" and
            pAtom2->getElement()->symbol != "N" and
            pBond->type != 2) { // not terminal
          for (int j = 0; j < abs(fc); j++) {
            if (fc < 0) {
              v++;
            }
            else {
              v--;
            }
          }
        }
      }
      /////
      errMessage += "\n ";

      if (v > 8) v = 8;
      pAtom->setValence(static_cast<int>(v));

      errMessage += i2s(pAtom->getFileID())
                + " | type = " + i2s(pAtom->getType())
                + " | valence = " + i2s(pAtom->getValence())
                + " | hybrid = " + i2s(pAtom->getHybridization())
                + " | charge = " + i2s(pAtom->getFormalCharge()) + "\n";
      v = 0;
//std::cout << " ------------------------ " << std::endl;
    }
    errorLogger.throwError("molecule::determineValences", errMessage, INFO);

    this->bAtomValencesAssigned = true;
}

// ============================================================
// Function : determineHybridizations()
// ------------------------------------------------------------
// Set all heavy atom hybridizations
// ============================================================
void molecule::determineHybridizations(int method, std::string paramFile)
{
    hybridize* pHybridizer = new hybridize(this, paramFile);

    if (method == 0) { // mdl, mol2 file formats where bond types are present
      if (inFileType == "mol" or inFileType == "sdf" or inFileType == "mol2") {
        pHybridizer->run();
      }
    }
    else if (method == 1) { // Meng Algorithm
      pHybridizer->runMeng();
    }
    else if (method == 2) { // Labute Algorithm
      pHybridizer->runLabute();
    }
    delete pHybridizer;
}

// ============================================================
// Function : addHydrogens()
// ------------------------------------------------------------
// add hydrogens
// ============================================================
void molecule::addHydrogens()
{
    errorLogger.throwError("molecule::addHydrogens", this->getName(), INFO);

    this->determineValences();
    if (!pProtonate) {
      pProtonate = new protonate(this);
    }
    pProtonate->run();
    this->bHydrogensAdded = true;
    delete pProtonate;
}

// ============================================================
// Function : removeHydrogens()
// ------------------------------------------------------------
// remove hydrogens
// ============================================================
void molecule::removeHydrogens()
{
    errorLogger.throwError("molecule::removeHydrogens", this->getName(), INFO);

    std::vector<atom*> Hs;
    std::vector<atom*> atomList = this->getAtomList();
    for (unsigned int i = 0; i < atomList.size(); i++) {
      if (atomList[i]->getElementSymbol() == "H") {
        Hs.push_back(atomList[i]);
      }
    }
    for (unsigned int i = 0; i < Hs.size(); i++) {
      this->delAtom(Hs[i]);
    }

    std::string errMessage = this->getName() + ": " + i2s(Hs.size()) + " Hydrogens to be deleted ";
    errorLogger.throwError("molecule::removeHydrogens", errMessage, INFO);
}

// ============================================================
// Function : getSubMolecule()
// ------------------------------------------------------------
//
// ============================================================
submolecule* molecule::getSubMolecule(int id)
{
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      if (pSubMolecule->getSubMolId() == id) {
        return pSubMolecule;
      }
    }
    return 0;
}

// ============================================================
// Function : getSubMolecule()
// ------------------------------------------------------------
//
// ============================================================
submolecule* molecule::getSubMolecule(int id, bool smolIndex, bool fileId)
{
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      if (smolIndex && (pSubMolecule->getIndex() == id)) {
        return pSubMolecule;
      }
      if (fileId && (pSubMolecule->getSubMolId() == id)) {
        return pSubMolecule;
      }
    }
    return 0;
}

// ============================================================
// Function : getSubMolecule()
// ------------------------------------------------------------
//
// ============================================================
submolecule* molecule::getSubMolecule(std::string name)
{
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      if (pSubMolecule->getName() == name) {
        return pSubMolecule;
      }
    }
    return 0;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* molecule::getAtom(int number, bool atomIndex, bool fileId, bool atomColIndex)
{
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        pAtom = *d;
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
    return 0;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* molecule::getAtom(stdAtom* pStdAtom)
{
    stdAtom* pStdAt = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        pAtom = *d;
        pStdAt = pAtom->getStdAtom();
        if (pStdAt == pStdAtom) {
          return pAtom;
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* molecule::getAtom(std::string name)
{
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (AtomIterator d = atList.begin(); d != atList.end(); d++) {
        pAtom = *d;
        if (pAtom->getName() == name) {
          return pAtom;
        }
      }
    }
    return 0;
}

// ============================================================
// Function : delAtom()
// ------------------------------------------------------------
//
// ============================================================
void molecule::delAtom(atom* pAt)
{
    if (!pAt) return;

    pSubMolecule = pAt->getParent();
    if (!pSubMolecule) return;

    errorLogger.throwError("molecule::delAtom", i2s(pAt->getFileID()), INFO);

    // Delete all bonds
    if (!this->itsBondMap.empty()) {
      BondMapIterator b;
      for (b = this->itsBondMap.begin(); b != this->itsBondMap.end(); b++) {
        pBond = b->second;
        if (pBond->atom1 == pAt or pBond->atom2 == pAt) {
          delete b->second;
          this->itsBondMap.erase(b);
        }
      }
    }

    // Delete all angles
    if (!this->itsAngleMap.empty()) {
      AngleMapIterator a;
      for (a = this->itsAngleMap.begin(); a != this->itsAngleMap.end(); a++) {
        pAngle = a->second;
        if (pAngle->atom1 == pAt or pAngle->atom2 == pAt or pAngle->atom3 == pAt) {
          delete a->second;
          this->itsAngleMap.erase(a);
        }
      }
    }

    // Delete all torsions
    if (!this->itsTorsionMap.empty()) {
      TorsionMapIterator t;
      for (t = this->itsTorsionMap.begin(); t != this->itsTorsionMap.end(); t++) {
        pTorsion = t->second;
        if (pTorsion->atom1 == pAt or pTorsion->atom2 == pAt or
            pTorsion->atom3 == pAt or pTorsion->atom4 == pAt) {
          delete t->second;
          this->itsTorsionMap.erase(t);
        }
      }
    }

    // Delete all impropers
    if (!this->itsImproperMap.empty()) {
      ImproperMapIterator i;
      for (i= this->itsImproperMap.begin(); i != this->itsImproperMap.end(); i++) {
        pImproper = i->second;
        if (pImproper->atom1 == pAt or pImproper->atom2 == pAt or
            pImproper->atom3 == pAt or pImproper->atom4 == pAt) {
          delete i->second;
          this->itsImproperMap.erase(i);
        }
      }
    }

    // If atom is in a ring, delete ring
    if (!this->itsConformers.empty()) {
      RingIterator r;
      for (r = this->itsRings.begin(); r != this->itsRings.end(); r++) {
        ring* cRing = *r;
        AtomIterator a = std::find(cRing->atoms.begin(), cRing->atoms.end(), pAt);
        if (a != cRing->atoms.end()) {
          delete *r;
          this->itsRings.erase(r);
        }
      }
    }

    // Delete all conformers
    if (!this->itsConformers.empty()) {
      ConformerIterator c;
      for (c = this->itsConformers.begin(); c != this->itsConformers.end(); c++) {
        delete *c;
        this->itsConformers.erase(c);
      }
    }

    // Clear fingerprint
    for (unsigned int i = 0; i < this->itsSimpleFP.size(); i++) {
      this->itsSimpleFP[i] = 0;
    }

    // Delete atom in submolecule
    AtomIterator a = std::find(pSubMolecule->itsAtomList.begin(), pSubMolecule->itsAtomList.end(), pAt);
    if (a != pSubMolecule->itsAtomList.end()) {
       delete *a;
       pSubMolecule->itsAtomList.erase(a);
    }
}

// ============================================================
// Function : moveCenterOfMass()
// ------------------------------------------------------------
//
// ============================================================
void molecule::moveCenterOfMass(vector3d* c)
{
    vector3d* Coord;

    vector3d  Center;
    std::vector<atom*> atomList = this->getAtomList();
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      Coord = pAtom->getCoords();
      Center = Center+(*Coord);
    }

    Center = Center / atomList.size();

    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      Coord = pAtom->getCoords();

      double x = Coord->getX() - Center.getX();
      double y = Coord->getY() - Center.getY();
      double z = Coord->getZ() - Center.getZ();

      pAtom->setCoords(x, y, z);
    }
}

// ============================================================
// Function : centerOfMass()
// ------------------------------------------------------------
//
// ============================================================
void molecule::centerOfMass(vector3d* c)
{
    vector3d* Coord;
    vector3d  Center;

    std::vector<atom*> atomList = this->getAtomList();
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      Coord = pAtom->getCoords();
      Center = Center+(*Coord);
    }

    Center = Center / atomList.size();
    *c =  Center;
}

// ============================================================
// Function : centerOfMass()
// ------------------------------------------------------------
//
// ============================================================
void molecule::centerOfMass(double c[3])
{
    vector3d* Coord;
    vector3d  Center;

    std::vector<atom*> atomList = this->getAtomList();
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      Coord = pAtom->getCoords();
      Center = Center+(*Coord);
    }

    Center = Center / atomList.size();
    c[0] =  Center[0];
    c[1] =  Center[1];
    c[2] =  Center[2];
}

// ============================================================
// Function : setCoordinates()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setCoordinates(std::vector< vector3d > &coords)
{
    std::vector<atom*> atomList = this->getAtomList();
    if (coords.size() == atomList.size()) {
      for (unsigned int i = 0; i < atomList.size(); i++) {
        pAtom = atomList[i];
        pAtom->setCoords(coords[i].getX(), coords[i].getY(), coords[i].getZ());
      }
    }
}

// ============================================================
// Function : getCoordinates()
// ------------------------------------------------------------
//
// ============================================================
void molecule::getCoordinates(double coords[][3])
{
    int i = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        coords[i][0] = atList[j]->getX();
        coords[i][1] = atList[j]->getY();
        coords[i][2] = atList[j]->getZ();
        i++;
      }
    }
}

// ============================================================
// Function : getCoordinates()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getCoordinates(double coords[])
{
    int i = 0;
    try {
      for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
        pSubMolecule = *c;
        std::vector<atom*> atList = pSubMolecule->getAtomList();
        for (unsigned int j = 0; j < atList.size(); j++) {
          coords[i  ] = atList[j]->getX();
          coords[i+1] = atList[j]->getY();
          coords[i+2] = atList[j]->getZ();
          i+=3;
        }
      }
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getCoordinates", " Memory Out of bounds Failure ", MTK_ERROR);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getCoordinates()
// ------------------------------------------------------------
//
// ============================================================
void molecule::getCoordinates(std::vector< vector3d > &coords)
{
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        coords.push_back((*atList[j]->getCoords()));
      }
    }
}

// ============================================================
// Function : getHeavyAtomIndices()
// ------------------------------------------------------------
//
// ============================================================
int* molecule::getHeavyAtomIndices()
{
    try {
      this->heavyAtomIndices = new int [this->getNumHeavyAtoms()];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getHeavyAtomIndices", " Memory Allocation Failure", MTK_ERROR);
      return 0;
    }

    int i = 0;
    int iIndex = 0;
    std::string currentSymbol = "";
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        if (!atList[j]->getElement()) {
          errorLogger.throwError("molecule::getHeavyAtomIndices", " Can't find element ", MTK_ERROR);
          return 0;
        }
        currentSymbol = atList[j]->getElement()->symbol;
        if (currentSymbol != "H") {
          this->heavyAtomIndices[iIndex] = i;
          iIndex++;
        }
        i++;
      }
    }
    return this->heavyAtomIndices;
}

// ============================================================
// Function : getAtomSymbols()
// ------------------------------------------------------------
//
// ============================================================
/*
void molecule::getAtomSymbols(char symbols[][2])
{
    int i = 0;
    string currentSymbol = "";
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        currentSymbol = atList[j]->getElement()->symbol;
        symbols[i][0] = currentSymbol[0];
        symbols[i][1] = currentSymbol[1];
        i++;
      }
    }
}
*/

// ============================================================
// Function : getRes1LSymbols()
// ------------------------------------------------------------
//
// ============================================================
void molecule::getRes1LSymbols(char codes1l[])
{
    int i = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::string oneLetterCode = pSubMolecule->get1LName();
      if (oneLetterCode.size() == 1) {
        codes1l[i] = oneLetterCode[0];
      }
      else {
        //codes1l[i] = ' ';
        codes1l[i] = 'X';
      }
      i++;
    }
}

// ============================================================
// Function : getRes3LSymbols()
// ------------------------------------------------------------
//
// ============================================================
void molecule::getRes3LSymbols(char codes3l[])
{
    int i = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::string resName = pSubMolecule->getName();
      if (resName.size() == 1) {
        codes3l[i] = resName[0];
        codes3l[i+1] = resName[0];
        codes3l[i+2] = resName[0];
      }
      else {
        codes3l[i] = ' ';
        codes3l[i+1] = ' ';
        codes3l[i+2] = ' ';
      }
      i+=3;
    }
}

// ============================================================
// Function : getAtomSymbols()
// ------------------------------------------------------------
//
// ============================================================
char* molecule::getAtomSymbols()
{
    try {
      this->atomSymbols = new char [this->getNumAtoms()*2];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getAtomSymbols", " Memory Allocation Failure", MTK_ERROR);
      return 0;
    }

    int i = 0;
    std::string currentSymbol = "";
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        currentSymbol = atList[j]->getElement()->symbol;
        this->atomSymbols[i  ] = currentSymbol[0];
        this->atomSymbols[i+1] = currentSymbol[1];
        i+=2;
      }
    }
    return this->atomSymbols;
}

// ============================================================
// Function : getHeavyAtomSymbols()
// ------------------------------------------------------------
//
// ============================================================
char* molecule::getHeavyAtomSymbols()
{
    try {
      this->heavyAtomSymbols = new char [this->getNumHeavyAtoms()*2];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getHeavyAtomSymbols", " Memory Allocation Failure", MTK_ERROR);
      return 0;
    }

    int i = 0;
    std::string currentSymbol = "";
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        currentSymbol = atList[j]->getElement()->symbol;
        if (currentSymbol == "H") continue;
        this->heavyAtomSymbols[i  ] = currentSymbol[0];
        this->heavyAtomSymbols[i+1] = currentSymbol[1];
        i+=2;
      }
    }
    return this->heavyAtomSymbols;
}

// ============================================================
// Function : getAtomTypes()
// ------------------------------------------------------------
//
// ============================================================
char* molecule::getAtomTypes()
{
    // Assuming all type are of size 2
    try {
      this->atomTypes = new char [this->getNumAtoms()*2];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getAtomTypes", " Memory Allocation Failure", MTK_ERROR);
      return 0;
    }

    int i = 0;
    std::string currentType = "";
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        stdAtom* pStdAtom = atList[j]->getStdAtom();
        if (pStdAtom) {
          currentType = pStdAtom->type;
          if (currentType.size() == 1) currentType+=" ";
          this->atomTypes[i  ] = currentType[0];
          this->atomTypes[i+1] = currentType[1];
          i+=2;
        }
        else {
          errorLogger.throwError("molecule::getAtomTypes", this->getName(), MTK_ERROR);
          delete this->atomTypes;
          this->atomTypes = 0;
          return 0;
        }
      }
    }
    return this->atomTypes;
}

// ============================================================
// Function : getHeavyAtomTypes()
// ------------------------------------------------------------
//
// ============================================================
char* molecule::getHeavyAtomTypes()
{
    // Assuming all type are of size 2
    try {
      this->heavyAtomTypes = new char [this->getNumHeavyAtoms()*2];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getHeavyAtomTypes", " Memory Allocation Failure", MTK_ERROR);
      return 0;
    }

    int i = 0;
    std::string currentType = "";
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        if (atList[j]->getElementSymbol() == "H") continue;
        stdAtom* pStdAtom = atList[j]->getStdAtom();
        if (pStdAtom) {
          currentType = pStdAtom->type;
          if (currentType.size() == 1) currentType+=" ";
          this->heavyAtomTypes[i  ] = currentType[0];
          this->heavyAtomTypes[i+1] = currentType[1];
          i+=2;
        }
        else {
          errorLogger.throwError("molecule::getHeavyAtomTypes", this->getName(), MTK_ERROR);
          delete this->heavyAtomTypes;
          this->heavyAtomTypes = 0;
          return 0;
        }
      }
    }
    return this->heavyAtomTypes;
}

// ============================================================
// Function : getAtomKinds()
// ------------------------------------------------------------
//
// ============================================================
int* molecule::getAtomKinds()
{
    try {
      atomKinds = new int [this->getNumAtoms()];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("molecule::getAtomKinds", " Memory Allocation Failure", MTK_ERROR);
      return 0;
    }

    int i = 0;
    int currentKind = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        currentKind = atList[j]->getType();
        atomKinds[i] = currentKind;
        currentKind = 0;
        i++;
      }
    }
    return this->atomKinds;
}

// ============================================================
// Function : getAtomCharges()
// ------------------------------------------------------------
//
// ============================================================
void molecule::getAtomCharges(double charges[])
{
    int i = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      std::vector<atom*> atList = pSubMolecule->getAtomList();
      for (unsigned int j = 0; j < atList.size(); j++) {
        charges[i] = double(atList[j]->getAtomicNum());
        i++;
      }
    }
}

// ============================================================
// Function : getFormalCharge()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getFormalCharge()
{
    int formalCharge = 0;
    if (bInSolution) {
      for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
        pSubMolecule = *c;
        std::vector<atom*> atList = pSubMolecule->getAtomList();
        for (unsigned int j = 0; j < atList.size(); j++) {
          formalCharge += atList[j]->getFormalCharge();
        }
      }
    }
    else {
      return this->getTotalCharge();
    }
    return formalCharge;
}

// ============================================================
// Function : getStdCharge()
// ------------------------------------------------------------
//
// ============================================================
double molecule::getStdCharge()
{
    double stdCharge = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      if (pSubMolecule->hasStdFrag()) {
        stdCharge += pSubMolecule->getStdFrag()->getCharge();
      }
    }
    return stdCharge;
}

// ============================================================
// Function : setTotalCharge()
// ------------------------------------------------------------
//
// ============================================================
void molecule::setTotalCharge(int d)
{
    this->totalCharge = d;
}

// ============================================================
// Function : getTotalCharge()
// ------------------------------------------------------------
//
// ============================================================
int molecule::getTotalCharge()
{
    return this->totalCharge;
}

// ============================================================
// Function : getMolecularWeight()
// ------------------------------------------------------------
//
// ============================================================
double molecule::getMolecularWeight()
{
    double w = 0;
    for (sMolIterator c = itsSubMoleculeList.begin(); c != itsSubMoleculeList.end(); c++) {
      pSubMolecule = *c;
      w += pSubMolecule->getMolecularWeight();
    }
    return w;
}

} // MTKpp namespace
