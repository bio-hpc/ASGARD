/*!
   \file connections.cpp
   \brief Assigns connectivity
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

#include "connections.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"
#include "angle.h"
#include "torsion.h"
#include "improper.h"
#include "metalCenter.h"
#include "utility.h"

#include "Utils/vector3d.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "parameters.h"

#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : connections()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
connections::connections()
{
    missingBondParams = 0;
    missingAngleParams = 0;
    missingTorsionParams = 0;
}
// ============================================================
// Function : connections()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
connections::connections(collection* pCol):pCollection(pCol)
{
    missingBondParams = 0;
    missingAngleParams = 0;
    missingTorsionParams = 0;
}

// ============================================================
// Function : ~connections()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
connections::~connections() {}


// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
void connections::run(collection* pCol)
{
    if (!pCol) {
      errorLogger.throwError("connections::run", " Error finding collection ", MTK_ERROR);
      return;
    }

    pCollection = pCol;
    this->run();
}

// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
void connections::run()
{
    if (!pCollection) {
      errorLogger.throwError("connections::run", " Error finding collection ", MTK_ERROR);
      return;
    }

    errorLogger.throwError("connections::run", " Collection Begin ", INFO);
    this->assignBonds();
    this->assignDisulfideBonds();
    this->assignAngles();
    this->assignTorsions();
    this->assignImpropers();

    this->assignStdBonds();
    this->assignStdAngles();
    this->assignStdTorsions();
    this->assignStdImpropers();
}

// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
void connections::run(molecule* pMol)
{
    errorLogger.throwError("connections::run", " Molecule Begin ", INFO);
    this->assignBonds(pMol);
    this->assignDisulfideBonds(pMol);
    this->assignAngles(pMol);
    this->assignTorsions(pMol);
    this->assignImpropers(pMol);

    // assign parameters
    this->assignStd(pMol);
}

// ============================================================
// Function : assignStd()
// ------------------------------------------------------------
//
// ============================================================
void connections::assignStd()
{
    errorLogger.throwError("connections::assignStd", " Collection Begin ", INFO);
    this->assignStdBonds();
    this->assignStdAngles();
    this->assignStdTorsions();
    this->assignStdImpropers();
}

// ============================================================
// Function : assignStd()
// ------------------------------------------------------------
//
// ============================================================
void connections::assignStd(molecule* pMol)
{
    errorLogger.throwError("connections::assignStd", " Molecule Begin ", INFO);
    this->assignStdBonds(pMol);
    this->assignStdAngles(pMol);
    this->assignStdTorsions(pMol);
    this->assignStdImpropers(pMol);
}

// ============================================================
// Function : assignStdBondsAngles()
// ------------------------------------------------------------
//
// ============================================================
void connections::assignStdBondsAngles(molecule* pMol)
{
    errorLogger.throwError("connections::assignStdBondsAngles", " Molecule Begin ", INFO);
    this->assignStdBonds(pMol);
    this->assignStdAngles(pMol);
}

// ============================================================
// Function : assignBonds()
// ------------------------------------------------------------
// Assigns all Bonds for Every Molecule in the Collection
// ============================================================
void connections::assignBonds()
{
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignBonds(pMolecule);
    }
}

// ============================================================
// Function : assignBonds()
// ------------------------------------------------------------
// Assigns all Bonds in a Molecule
// ============================================================
void connections::assignBonds(molecule* pMol)
{
    errorMessage = " Molecule Begin # " + i2s(pMol->getMolId()) + "\n";
    pStdFragMinus1 = 0;
    pSubMoleculeMinus1 = 0;
    subMoleculeList = pMol->getSubMoleculeList();
    int nBondsBefore = 0;
    int nBondsAfter = 0;
    int nBondsCreated = 0;

    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];
      // Assigning bonds for current smol and
      // bonds between it and the previous residue
      pStdFrag = pSubMolecule->getStdFrag();
      nBondsBefore = pMol->getNumBonds();
      if (pStdFrag) {
        std::string stdFragType = pStdFrag->getType();
        if ((stdFragType == "s") or (stdFragType == "m")) {
          pSubMoleculeMinus1 = 0;
        }
        if (pStdFragMinus1) {
          if (pStdFragMinus1->getType() == "e") {
            pSubMoleculeMinus1 = 0;
            pStdFragMinus1 = 0;
          }
        }

        if (pStdFragMinus1) {
          this->bondByLibrary(pMol, pSubMolecule, pStdFrag,
                              pSubMoleculeMinus1, pStdFragMinus1);
        }
        else {
          this->bondByLibrary(pMol, pSubMolecule, pStdFrag, 0, 0);
        }
      }
      else {
        for (unsigned int j2 = 0; j2 < subMoleculeList.size(); j2++) {
          //if (subMoleculeList[j] == subMoleculeList[j2]) continue;
          this->bondByDistance(pMol, pSubMolecule, subMoleculeList[j2]);
        }
      }
      nBondsAfter = pMol->getNumBonds();
      nBondsCreated = nBondsAfter - nBondsBefore;
	    errorMessage += pSubMolecule->getName() + ": " + i2s(nBondsCreated);
	    if (nBondsCreated == 1) {
		    errorMessage += " bond was added.\n";
	    } else {
		    errorMessage += " bonds were added.\n";
	    }
      pSubMoleculeMinus1 = pSubMolecule;
      pStdFragMinus1     = pStdFrag;
    }
    pSubMoleculeMinus1 = 0;
    pStdFragMinus1 = 0;
    pMol->bBondsAssigned = true;
    errorLogger.throwError("connections::assignBonds", errorMessage, INFO);
}

// ============================================================
// Function : bondByLibrary()
// ------------------------------------------------------------
// Creates bonds in a molecule using bond-by-library
// ============================================================
void connections::bondByLibrary(molecule* pMolecule,
                                submolecule* pSubMolecule, stdFrag* pStdFrag,
                                submolecule* pSubMoleculeMinus1, stdFrag* pStdFragMinus1)
{
    //std::cout << " connections::bondByLibrary" <<std::endl;
    int nBondsAdded = 0;
    atomList = pSubMolecule->getAtomList();
    double distance = 0.0;
    stdBond* pStdBond = 0;

    // intra-residue bonds
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom1 = atomList[i];
      pStdAtom1 = pAtom1->getStdAtom();
      if (pStdAtom1) {
        for (unsigned int j = 0; j < atomList.size(); j++) {
          pAtom2 = atomList[j];
          if (pAtom1 == pAtom2) continue;
          pStdAtom2 = pAtom2->getStdAtom();
          if (pStdAtom2) {
            pStdBond = pStdFrag->getStdBond(pStdAtom1->index, pStdAtom2->index);
            if  (pStdBond) {
              if (!pMolecule->hasBond(pAtom1, pAtom2)) {
                coord1 = pAtom1->getCoords();
                coord2 = pAtom2->getCoords();
                distance = coord1->dist(*coord2);
                pBond = pMolecule->addBond(pAtom1, pAtom2, pStdBond->type,
                        pStdBond->stereo, pStdBond->topology, distance);
                        //std::cout << " " << pSubMolecule->getName() << " " << pAtom1->getName() << " " 
                        //    << pAtom2->getName() << std::endl;
                nBondsAdded++;
                pAtom1->addBondedAtom(pAtom2);
                pAtom2->addBondedAtom(pAtom1);
                pStdBond = 0;
              }
            }
          }
          else {
            errorMessage = " Missing standard atom info for '" + pAtom2->getName() + "' in " + pStdFrag->getName();
            errorLogger.throwError("connections::bondByLibrary", errorMessage, MTK_ERROR);
          }
        }
      }
      else {
        errorMessage = " Missing standard atom info for '" + pAtom2->getElement()->symbol + "' in " + pStdFrag->getName();
        errorLogger.throwError("connections::bondByLibrary", errorMessage, MTK_ERROR);
      }
    }

    // inter-residue bonds
    int tempNumber = 0;
    atom* pBondAtom = 0;
    if (pSubMoleculeMinus1) {
      pStdFragMinus1 = pSubMoleculeMinus1->getStdFrag();
      atomList2      = pSubMoleculeMinus1->getAtomList();
      if (pStdFragMinus1) {
        for (unsigned int i = 0; i < atomList.size(); i++) {
          pAtom1 = atomList[i];
          pStdAtom1 = pAtom1->getStdAtom();
          if (pStdAtom1) {
            if (pStdAtom1->bond12 < 0) {
              for (unsigned int j = 0; j < atomList2.size(); j++) {
                pAtom2 = atomList2[j];
                pStdAtom2 = pAtom2->getStdAtom();
                if (pStdAtom2) {
                  if (pStdAtom2->chain == "M") {
                    if (pStdAtom2->index > tempNumber) {
                      tempNumber = pStdAtom2->index;
                      pBondAtom = pAtom2;
                    }
                  }
                }
              }

              if (!pMolecule->hasBond(pAtom1, pBondAtom)){
                coord1 = pAtom1->getCoords();
                coord2 = pBondAtom->getCoords();
                distance = coord1->dist(*coord2);
                if (pMolecule->getName() == "Reference") {
                  //std::cout << " INTER BOND " << pAtom1->getParent()->getName() << " " << pAtom1->getName() << " " 
                  //          << pAtom2->getParent()->getName() << " " << pAtom2->getName() << std::endl;
                  continue;
                }
                pMolecule->addBond(pAtom1, pBondAtom, 0, 0, 0, distance);
                pAtom1->addBondedAtom(pBondAtom);
                pBondAtom->addBondedAtom(pAtom1);
                nBondsAdded++;
              }
            }
          }
        }
      }
    }

    // Loop Atoms
    stdLoopList = pStdFrag->getStdLoopList();
    for (unsigned int k = 0; k < stdLoopList.size(); k++) {
      pStdLoop = stdLoopList[k];
      pStdAtom1 = pStdFrag->getStdAtom(pStdLoop->atom1);
      pStdAtom2 = pStdFrag->getStdAtom(pStdLoop->atom2);
      if (pStdAtom1 && pStdAtom2) {
        pAtom1 = pSubMolecule->getAtom(pStdAtom1->identity);
        pAtom2 = pSubMolecule->getAtom(pStdAtom2->identity);

//std::cout << pStdFrag->getName() <<std::endl;
//pStdFrag->print();
        if (pAtom1 and pAtom2) {
//std::cout << "   " << pAtom1->getName() << " " << pAtom2->getName() << std::endl;

          if (!pMolecule->hasBond(pAtom1, pAtom2)) {
            coord1 = pAtom1->getCoords();
            coord2 = pAtom2->getCoords();
            distance = coord1->dist(*coord2);
            pMolecule->addBond(pAtom1, pAtom2, pStdLoop->type, pStdLoop->stereo, 1, distance);
//std::cout << "   " << pAtom1->getName() << " " << pAtom2->getName() << std::endl;
            pAtom1->addBondedAtom(pAtom2);
            pAtom2->addBondedAtom(pAtom1);
            nBondsAdded++;
          }
        }
        else {
          //pSubMolecule->print();
        }
      }
    }
}

// ============================================================
// Function : bondByDistance()
// ------------------------------------------------------------
// Creates bonds in a molecule using bond-by-distance
// ============================================================
void connections::bondByDistance(molecule* pMolecule, submolecule* pSubMolecule,
                                 submolecule* pSubMoleculeMinus1)
{
    //std::cout << "connections::bondByDistance" <<std::endl;
    bool checkBond = false;
    double distance = 0.0;
    int bondType = 1;
    if (pMolecule && !(pSubMolecule or pSubMoleculeMinus1)) {
      atomList = pMolecule->getAtomList();
      for (unsigned int i = 0; i < atomList.size(); i++) {
        pAtom1 = atomList[i];
        coord1 = pAtom1->getCoords();
        for (unsigned int j = i+1; j < atomList.size(); j++) {
          pAtom2 = atomList[j];
          coord2 = pAtom2->getCoords();
          distance = coord1->dist(*coord2);
          checkBond = BondExists(distance, pAtom1, pAtom2);
          if (checkBond) {
            if (!pMolecule->hasBond(pAtom1, pAtom2)) {
              pMolecule->addBond(pAtom1, pAtom2, bondType, 0, 0, distance);
              pAtom1->addBondedAtom(pAtom2);
              pAtom2->addBondedAtom(pAtom1);
            }
          }
        }
      }
    }

    if (pSubMolecule) {
      atomList = pSubMolecule->getAtomList();
      for (unsigned int i = 0; i < atomList.size(); i++) {
        pAtom1 = atomList[i];
        coord1 = pAtom1->getCoords();
        for (unsigned int j = i+1; j < atomList.size(); j++) {
          pAtom2 = atomList[j];

          coord2 = pAtom2->getCoords();
          distance = coord1->dist(*coord2);
          checkBond = BondExists(distance, pAtom1, pAtom2);

          //std::cout << pAtom1->getFileID() << "-" << pAtom2->getFileID() <<  " " << distance << " " << checkBond << std::endl;

          if (checkBond) {
            if (!pMolecule->hasBond(pAtom1, pAtom2)) {
              pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, distance);
              pAtom1->addBondedAtom(pAtom2);
              pAtom2->addBondedAtom(pAtom1);
            }
          }
        }
      }
    }

    if (pSubMolecule && pSubMoleculeMinus1) {
      if (pSubMolecule != pSubMoleculeMinus1) {
        atomList = pSubMolecule->getAtomList();
        atomList2 = pSubMoleculeMinus1->getAtomList();
        for (unsigned int i = 0; i < atomList.size(); i++) {
          pAtom1 = atomList[i];
          coord1 = pAtom1->getCoords();
          for (unsigned int j = 0; j < atomList2.size(); j++) {
            pAtom2 = atomList2[j];
            coord2 = pAtom2->getCoords();
            distance = coord1->dist(*coord2);
            checkBond = BondExists(distance, pAtom1, pAtom2);
            if (checkBond){
              if (!pMolecule->hasBond(pAtom1, pAtom2)) {
                pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, distance);
                pAtom1->addBondedAtom(pAtom2);
                pAtom2->addBondedAtom(pAtom1);
              }
            }
          }
        }
      }
    }
}

// ============================================================
// Function : BondExists()
// ------------------------------------------------------------
// Algorithm: J. Comp. Chem. 12, 891-898, 1991
// ============================================================
bool connections::BondExists(double distance, atom* pAtom1, atom* pAtom2)
{
    double covalentBondRadii1 = pAtom1->getElement()->covalentRadius;
    double covalentBondRadii2 = pAtom2->getElement()->covalentRadius;
    double d = covalentBondRadii1 + covalentBondRadii2 + BONDTOLERANCE;
/*
std::cout << "connections::BondExists " 
          << pAtom1->getElement()->symbol << "-" << pAtom2->getElement()->symbol << " rad1="
          << covalentBondRadii1 << " rad2=" << covalentBondRadii2 << " d=" << d << " dist=" << distance << std::endl;
*/
    if (distance < 0.5) {
      std::string atSymbol1 = pAtom1->getElement()->symbol;
      std::string atSymbol2 = pAtom2->getElement()->symbol;

      errorMessage = " Atom clash. The distance between "
                   + i2s(pAtom1->getFileID()) + "-" + atSymbol1 + " and " + i2s(pAtom2->getFileID())
                   + "-" + atSymbol2 + " is " + d2s(distance);
      errorLogger.throwError("connections::BondExists", errorMessage, INFO);

      return false;
    }
    if (distance <= d) {
      return true;
    }
    return false;
}

// ============================================================
// Function : assignDisulfideBonds()
// ------------------------------------------------------------
// Determine disulfide bonds
//
// Original Algorithm and S-S Parameters from Andrew Wollacott
//
// ============================================================
void connections::assignDisulfideBonds()
{
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignDisulfideBonds(pMolecule);
    }
}

// ============================================================
// Function : assignDisulfideBonds()
// ------------------------------------------------------------
// Determine disulfide bonds
//
// Original Algorithm and S-S Parameters from Andrew Wollacott
//
// ============================================================
void connections::assignDisulfideBonds(molecule* pMol)
{
    double Distcutoff  = 2.5;
    double SSdist      = 2.038;
    double SSdistK     = 166.0;
    double CT_SSangle  = 103.7;
    double CT_SSangleK = 68.0;
    double Esscutoff   = 30.0;

    double dist = 0;
    double cb1_sg1_sg2 = 0;
    double cb2_sg2_sg1 = 0;
    double ssEnergy = 0;

    submolecule* cys1;
    submolecule* cys2;
    atom* sg1;
    atom* sg2;
    atom* cb1;
    atom* cb2;

    pMol->bDisulfideBondsAssigned = 1;
    subMoleculeList = pMol->getSubMoleculeList();

    std::vector<submolecule*> allCys;
    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];
      if ((pSubMolecule->getName() == "CYS") or (pSubMolecule->getName() == "CYX")) {
        allCys.push_back(pSubMolecule);
      }
    }

    for (unsigned int k = 0; k < allCys.size(); k++) {
      cys1 = allCys[k];
      sg1 = cys1->getAtom(" SG ");
      cb1 = cys1->getAtom(" CB ");

      if (!sg1 or !cb1) {
        errorMessage  = "\n ### CYS Residue: " + i2s(allCys[k]->getSubMolId());
        errorMessage += "\n###  Can't find SG and/or CB ";
        errorLogger.throwError("connections::assignDisulfideBonds", errorMessage, WARNING);
        continue;
      }

      coord1 = sg1->getCoords();
      for ( unsigned int l = k+1; l < allCys.size(); l++) {
        cys2 = allCys[l];
        sg2 = cys2->getAtom(" SG ");
        cb2 = cys2->getAtom(" CB ");

        if (!sg2 or !cb2) {
          continue;
        }

        coord2 = sg2->getCoords();

        // Check distance between CYS's
        dist = coord1->dist(*coord2);
        if (dist < Distcutoff) {
          // Check Angles
          cb1_sg1_sg2 = RAD2DEG * angle(*(cb1->getCoords()), *(sg1->getCoords()), *(sg2->getCoords()));
          cb2_sg2_sg1 = RAD2DEG * angle(*(cb2->getCoords()), *(sg2->getCoords()), *(sg1->getCoords()));
          ssEnergy = 0.0;
          ssEnergy = SSdistK * pow(dist-SSdist,2);
          ssEnergy = ssEnergy + CT_SSangleK * pow((cb1_sg1_sg2-CT_SSangle)/RAD2DEG,2);
          ssEnergy = ssEnergy + CT_SSangleK * pow((cb2_sg2_sg1-CT_SSangle)/RAD2DEG,2);
          if (ssEnergy < Esscutoff) {
            if (!pMol->hasBond(sg1, sg2)) {
              Bond* pSSBond = pMol->addBond(sg1, sg2, 1, 0, 2, dist);
              pSSBond->kind = 3;
              sg1->addBondedAtom(sg2);
              sg2->addBondedAtom(sg1);

              errorMessage = " Molecule Name: " + i2s(pMol->getMolId());
              errorMessage += " Assigning disulfide bond: " +
                 i2s(sg1->getIndex()) + "-" +
                 i2s(sg2->getIndex()) + "\n";
              cys1->setName("CYX");
              cys2->setName("CYX");
              errorLogger.throwError("connections::assignDisulfideBonds", errorMessage, INFO);
            }
          }
        }
      }
    }
}

// ============================================================
// Function : assignStdBonds()
// ------------------------------------------------------------
//
// ============================================================
void connections::assignStdBonds()
{
    errorLogger.throwError("connections::assignStdBonds", " Collection Begin ", INFO);
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignStdBonds(pMolecule);
    }
}

// ============================================================
// Function : assignStdBonds()
// ------------------------------------------------------------
//
// ============================================================
void connections::assignStdBonds(molecule* pMol)
{
    errorMessage = " Molecule # " + i2s(pMol->getMolId()) + "\n";

    typedef std::map<int, Bond*>::iterator BondMapIterator;
    std::map<int, Bond*> moleculeBondMap;

    pParameters = pCollection->getParameters();
    bondParam* pBondParam;

    moleculeBondMap = pMol->getBondMap();

    if (pParameters and (!moleculeBondMap.empty()) and (pMol->getNumAtoms() > 1)) {
      for (BondMapIterator b = moleculeBondMap.begin(); b != moleculeBondMap.end(); b++) {
        pBond = b->second;

        bool errorOccured = false;
        if (!pBond->atom1->getStdAtom()) {
          errorOccured = true;
        }
        if (!pBond->atom2->getStdAtom()) {
          errorOccured = true;
        }

        if (errorOccured) {
          errorMessage += " Unable to assign parameters for bond : "
                       + pBond->atom1->getName() + "-" + pBond->atom2->getName() + "\n";
          missingBondParams++;
          continue;
        }

        pBondParam = pParameters->getBondParam(pBond->atom1->getStdAtom()->type,
                               pBond->atom2->getStdAtom()->type);
        if (pBondParam) {
          pBond->pBondParam = pBondParam;
        }
        else {
          errorMessage += " No Bond Parameter for the bond containing: "
                       + pBond->atom1->getName() + ":" + pBond->atom1->getStdAtom()->type + " <--> "
                       + pBond->atom2->getName() + ":" + pBond->atom2->getStdAtom()->type + "\n";
          missingBondParams++;
        }
      }
    }
    else {
      if (pMol->getNumAtoms() > 1) {
        errorMessage += " Please assign bonds before assigning MM parameters for molecule: "
                     + pMol->getName() + "\n";
      }
    }
    if (missingBondParams) {
      errorMessage += " Total Number of Missing Bond Parameters = "
                   + i2s(missingBondParams) + "\n";
    }
    errorLogger.throwError("connections::assignStdBonds", errorMessage, INFO);
}

// ============================================================
// Function : assignAngles()
// ------------------------------------------------------------
// Assigns all Angles for Every Molecule in the Collection
// ============================================================
void connections::assignAngles()
{
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignAngles(pMolecule);
    }

    std::vector<metalCenter*> metalCenterList = pCollection->getMetalCenters();
    for (unsigned int i = 0; i < metalCenterList.size(); i++) {
      metalCenter* pMetalCenter = metalCenterList[i];
      this->assignAngles(pMetalCenter);
    }
}

// ============================================================
// Function : assignAngles()
// ------------------------------------------------------------
// Assigns all Angles in a Molecule
// ============================================================
void connections::assignAngles(molecule* pMol)
{
    errorMessage = " Molecule # " + i2s(pMol->getMolId()) + " ";

    int nAnglesAdded = 0;
    bool gotAngle;
    double dAngle = 0.0;

    pMol->bAnglesAssigned = true;
    subMoleculeList = pMol->getSubMoleculeList();
    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];

      // Get atoms in submolecule
      atomList = pSubMolecule->getAtomList();
      for (unsigned int k = 0; k < atomList.size(); k++) {
        pAtom1 = atomList[k];
        for (unsigned int l = 0; l < pAtom1->bondedAtoms.size(); l++) {
          pAtom2 = pAtom1->bondedAtoms[l];
          if (pAtom1->getParent()->getParent() != pAtom2->getParent()->getParent()) continue;
          for (unsigned int o = 0; o < pAtom2->bondedAtoms.size(); o++) {
            pAtom3 = pAtom2->bondedAtoms[o];
            if (pAtom1->getParent()->getParent() != pAtom3->getParent()->getParent()) continue;
            if (pAtom1 != pAtom3) {
              gotAngle = pMol->hasAngle(pAtom1, pAtom2, pAtom3);
              if (!gotAngle) {
                dAngle = angle(*(pAtom1->getCoords()), *(pAtom2->getCoords()), *(pAtom3->getCoords()));
                pAngle = pMol->addAngle(pAtom1, pAtom2, pAtom3, dAngle);
                pAtom1->addBonded13Atom(pAtom3);
                pAtom3->addBonded13Atom(pAtom1);
                pMol->add13Bond(pAtom1, pAtom3);
                nAnglesAdded++;
              }
            }
          }
        }
      }
    }

	errorMessage += i2s(nAnglesAdded);
	if (nAnglesAdded == 1) {
		errorMessage += " angle was added";
	} else {
		errorMessage += " angles were added";
	}
    errorLogger.throwError("connections::assignAngles", errorMessage, INFO);
}

// ============================================================
// Function : assignAngles()
// ------------------------------------------------------------
// Assigns all Angles in a metal center
// ============================================================
void connections::assignAngles(metalCenter* pMetalCenter)
{
    errorMessage = " Metal center: " + pMetalCenter->getName();
    errorLogger.throwError("connections::assignAngles", errorMessage, INFO);

    bool gotAngle;
    double dAngle = 0.0;

    // Get atoms in metal center
    pAtom1 = pMetalCenter->getMetalAtom();
    for (unsigned int l = 0; l < pAtom1->bondedAtoms.size(); l++) {
      pAtom2 = pAtom1->bondedAtoms[l];
      for (unsigned int o = 0; o < pAtom2->bondedAtoms.size(); o++) {
        pAtom3 = pAtom2->bondedAtoms[o];
        if (pAtom1 != pAtom3) {
          gotAngle = pMetalCenter->hasAngle(pAtom1, pAtom2, pAtom3);
          if (!gotAngle) {
            errorMessage = " Adding angle " + pAtom1->getName() + " "
                         +  pAtom2->getName() + " " +  pAtom3->getName();
            errorLogger.throwError("connections::assignAngles", errorMessage, INFO);

            dAngle = angle(*(pAtom1->getCoords()),
                           *(pAtom2->getCoords()),
                           *(pAtom3->getCoords()));
            pAngle = pMetalCenter->addAngle(pAtom1, pAtom2, pAtom3, dAngle);
            pAtom1->addBonded13Atom(pAtom3);
            pAtom3->addBonded13Atom(pAtom1);
            //pMol->add13Bond(pAtom1, pAtom3);
          }
        }
      }
    }
}

// ============================================================
// Function : assignStdAngles()
// ------------------------------------------------------------
// Assigns parameters to all Angles for Every Molecule in the Collection
// ============================================================
void connections::assignStdAngles()
{
    errorLogger.throwError("connections::assignStdAngles", " Collection Begin ", INFO);
    moleculeList = pCollection->getMoleculeList();

    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignStdAngles(pMolecule);
    }
}

// ============================================================
// Function : assignStdAngles()
// ------------------------------------------------------------
// Assigns parameters to all Angles for Every Molecule in the Collection
// ============================================================
void connections::assignStdAngles(molecule* pMol)
{
    errorMessage = " Molecule # " + i2s(pMol->getMolId());
    errorLogger.throwError("connections::assignStdAngles", errorMessage, INFO);

    typedef std::map<ULONG_KIND, Angle*>::iterator AngleMapIterator;
    std::map<ULONG_KIND, Angle*>      moleculeAngleMap;

    angleParam* pAngleParam;
    pParameters = pCollection->getParameters();
    moleculeAngleMap = pMol->getAngleMap();

    std::string angleErrorMessage = "\n";

    if ((!moleculeAngleMap.empty()) and (pMol->getNumAtoms() > 2)) {
      for (AngleMapIterator a = moleculeAngleMap.begin(); a != moleculeAngleMap.end();a++) {
        pAngle = a->second;

        bool errorOccured = false;
        if (!pAngle->atom1->getStdAtom()) {
          errorOccured = true;
        }
        if (!pAngle->atom2->getStdAtom()) {
          errorOccured = true;
        }

        if (!pAngle->atom3->getStdAtom()) {
          errorOccured = true;
        }

        if (errorOccured) {
          errorMessage = " Unable to assign parameters for angle: " +
                         pAngle->atom1->getName() + "-" +
                         pAngle->atom2->getName() + "-" +
                         pAngle->atom3->getName() + "\n";
          angleErrorMessage += errorMessage;
          missingAngleParams++;
          continue;
        }

        pAngleParam = pParameters->getAngleParam(pAngle->atom1->getStdAtom()->type,
                                                 pAngle->atom2->getStdAtom()->type,
                                                 pAngle->atom3->getStdAtom()->type);
        if (pAngleParam) {
          pAngle->pAngleParam = pAngleParam;
        }
        else {
          errorMessage = " No angle parameter for the angle containing: |" +
                         pAngle->atom1->getName() + "|-|" +
                         pAngle->atom2->getName() + "|-|" +
                         pAngle->atom3->getName() + "|#|" +
                         pAngle->atom1->getStdAtom()->type + "|-|" +
                         pAngle->atom2->getStdAtom()->type + "|-|" +
                         pAngle->atom3->getStdAtom()->type + "|\n";
          angleErrorMessage += errorMessage;
          missingAngleParams++;
        }
      }
    }
    else {
      if (pMol->getNumAtoms() > 2) {
        errorMessage = " Please assign angles before assigning MM parameters for molecule: " + pMol->getName();
        errorLogger.throwError("connections::assignStdAngles", errorMessage, MTK_ERROR);
      }
    }
    if (missingAngleParams) {
      errorMessage = " Total Number of Missing Angle Parameters = " +
                     i2s(missingAngleParams) + "\n";
      angleErrorMessage += errorMessage;
      errorLogger.throwError("connections::assignStdAngles", angleErrorMessage, INFO);
    }
}

// ============================================================
// Function : assignTorsions()
// ------------------------------------------------------------
// Assigns all Torsions for Every Molecule in the Collection
// ============================================================
void connections::assignTorsions()
{
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignTorsions(pMolecule);
    }
}

// ============================================================
// Function : assignTorsions()
// ------------------------------------------------------------
// Assigns all Torsions for Every Molecule in the Collection
// ============================================================
void connections::assignTorsions(molecule* pMol)
{
    errorMessage = " Molecule # " + i2s(pMol->getMolId());
    errorLogger.throwError("connections::assignTorsions", errorMessage, INFO);

    int nTorsionsAdded = 0;
    bool gotTorsion;
    double dTorsion;

    pMol->bTorsionsAssigned = true;

    subMoleculeList = pMol->getSubMoleculeList();
    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];
      atomList = pSubMolecule->getAtomList();
      for (unsigned int k = 0; k < atomList.size(); k++) {
        pAtom1 = atomList[k];
        for (unsigned int l = 0; l < pAtom1->bondedAtoms.size(); l++) {
          pAtom2 = pAtom1->bondedAtoms[l];
          for (unsigned int o = 0; o < pAtom2->bondedAtoms.size(); o++) {
            pAtom3 = pAtom2->bondedAtoms[o];
            if (pAtom1 != pAtom3) {
              for (unsigned int p = 0; p < pAtom3->bondedAtoms.size(); p++) {
                pAtom4 = pAtom3->bondedAtoms[p];
                if (pAtom4 == pAtom2) continue;
                if (pAtom4 == pAtom1) continue;
                gotTorsion = pMol->hasTorsion(pAtom1, pAtom2, pAtom3, pAtom4);
                if (!gotTorsion) {
                  dTorsion = torsion(*(pAtom1->getCoords()), *(pAtom2->getCoords()),
                                     *(pAtom3->getCoords()), *(pAtom4->getCoords()));
                  pMol->addTorsion(pAtom1, pAtom2, pAtom3, pAtom4, dTorsion);
                  pAtom1->addBonded14Atom(pAtom4);
                  pAtom4->addBonded14Atom(pAtom1);
                  pMol->add14Bond(pAtom1, pAtom4);
                  nTorsionsAdded++;
                }
              }
            }
          }
        }
      }
    }
	errorMessage = i2s(nTorsionsAdded);
	if (nTorsionsAdded == 1) {
		errorMessage += " torsion was added";
	} else {
		errorMessage += " torsions were added";
	}
    errorLogger.throwError("connections::assignTorsions", errorMessage, INFO);
}

// ============================================================
// Function : assignStdTorsions()
// ------------------------------------------------------------
// Assigns parameters to all Torsions for Every Molecule in the Collection
// ============================================================
void connections::assignStdTorsions()
{
    errorLogger.throwError("connections::assignStdTorsions", " Collection Begin ", INFO);
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->assignStdTorsions(pMolecule);
    }
}

// ============================================================
// Function : assignStdTorsions()
// ------------------------------------------------------------
// Assigns parameters to all Torsions for Every Molecule in the Collection
// ============================================================
void connections::assignStdTorsions(molecule* pMol)
{
    errorMessage = " Molecule # " + i2s(pMol->getMolId());
    errorLogger.throwError("connections::assignStdTorsions", errorMessage, INFO);

    typedef std::map<ULONG_KIND, Torsion*>::iterator TorsionMapIterator;
    std::map<ULONG_KIND, Torsion*> moleculeTorsionMap;
    std::vector<torsionParam*> torsionParamList;

    torsionParam* pTorsionParam;

    pParameters = pCollection->getParameters();

    moleculeTorsionMap = pMol->getTorsionMap();

    std::string torsionErrorMessage = "\n";

    if ((!moleculeTorsionMap.empty()) and (pMol->getNumAtoms() > 3)) {
      for (TorsionMapIterator a = moleculeTorsionMap.begin();
           a != moleculeTorsionMap.end(); a++) {
        pTorsion = a->second;

        bool errorOccured = false;
        if (!pTorsion->atom1->getStdAtom()) {
          errorOccured = true;
        }
        if (!pTorsion->atom2->getStdAtom()) {
          errorOccured = true;
        }
        if (!pTorsion->atom3->getStdAtom()) {
          errorOccured = true;
        }
        if (!pTorsion->atom4->getStdAtom()) {
          errorOccured = true;
        }

        if (errorOccured) {
          errorMessage = " Unable to assign parameters for torsion: "
                       + pTorsion->atom1->getName() + "-"
                       + pTorsion->atom2->getName() + "-"
                       + pTorsion->atom3->getName() + "-"
                       + pTorsion->atom4->getName() + "\n";
          torsionErrorMessage += errorMessage;
          missingTorsionParams++;
          continue;
        }

        torsionParamList = pParameters->getTorsionParamList(
         pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
         pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

        if (!torsionParamList.empty()) {
          pTorsion->pTorsionParamList = torsionParamList;
          pTorsion->parametersAssigned = true;
        }
        else {
          errorMessage = "No Torsion Parameter for the torsion containing: "
                       + pTorsion->atom1->getParent()->getName() + "-"
                       + pTorsion->atom2->getParent()->getName() + "-"
                       + pTorsion->atom3->getParent()->getName() + "-"
                       + pTorsion->atom4->getParent()->getName() + " ::: "
                       + pTorsion->atom1->getName() + " - " + pTorsion->atom2->getName() + " - "
                       + pTorsion->atom3->getName() + " - " + pTorsion->atom4->getName() + " ::: "
                       + pTorsion->atom1->getStdAtom()->type + " - " + pTorsion->atom2->getStdAtom()->type + " - "
                       + pTorsion->atom3->getStdAtom()->type + " - " + pTorsion->atom4->getStdAtom()->type
                       + " " + i2s(pTorsion->atom1->getParent()->getSubMolId()) + "\n";
          torsionErrorMessage += errorMessage;
          missingTorsionParams++;
        }
      }
    }
    else {
      if (pMol->getNumAtoms() > 3) {
        errorMessage = "Please assign torsions before assigning MM parameters for molecule: " + pMol->getName();
        errorLogger.throwError("connections::assignStdTorsions", errorMessage, MTK_ERROR);
      }
    }

    std::vector<torsionParam*> toBeDeleted;
    std::vector<torsionParam*>::iterator vecIterator;

    moleculeTorsionMap = pMol->getTorsionMap();
    if (!moleculeTorsionMap.empty()) {
      // LOOP OVER ALL TORSIONS
      for (TorsionMapIterator a = moleculeTorsionMap.begin();
                          a != moleculeTorsionMap.end(); a++) {
        pTorsion = a->second;

        // DEAL WITH PRO
        if (pTorsion->atom1->getParent()->getName() == "PRO") {
          for (unsigned int k = 0; k < pTorsion->pTorsionParamList.size(); k++ ) {
            pTorsionParam = pTorsion->pTorsionParamList[k];
            if ((pTorsionParam->atomType1 == "CT") and (pTorsionParam->atomType2 == "CT") and
                (pTorsionParam->atomType3 == "N" ) and (pTorsionParam->atomType4 == "C" ) ) {
              if ((pTorsionParam->Nt == -4) or (pTorsionParam->Nt == 1)) {
                toBeDeleted.push_back(pTorsionParam);
              }
            }
          }
        }

        if (!toBeDeleted.empty()) {
          std::sort( toBeDeleted.begin(), toBeDeleted.end() );
          for (unsigned int n = 0; n < toBeDeleted.size(); n++) {
            vecIterator = std::find(pTorsion->pTorsionParamList.begin(),
                          pTorsion->pTorsionParamList.end(),toBeDeleted[n]);
            pTorsion->pTorsionParamList.erase(vecIterator);
          }
        }
        toBeDeleted.clear();
      }
    }
    if (missingTorsionParams) {
      errorMessage = " Total Number of Missing Torsion Parameters = " +
                     i2s(missingTorsionParams) + "\n";
      torsionErrorMessage += errorMessage;
      errorLogger.throwError("connections::assignStdTorsions", torsionErrorMessage, INFO);
    }
}

// ============================================================
// Function : assignImpropers()
// ------------------------------------------------------------
// Assigns all Impropers for Every Molecule in the Collection
// ============================================================
void connections::assignImpropers()
{
    moleculeList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < moleculeList.size(); i++){
      pMolecule = moleculeList[i];
      this->assignImpropers(pMolecule);
    }
}

// ============================================================
// Function : assignImpropers()
// ------------------------------------------------------------
// Assigns all Impropers in a Molecule
// ============================================================
void connections::assignImpropers(molecule* pMol)
{
    errorMessage = " Molecule # " + i2s(pMol->getMolId());
    errorLogger.throwError("connections::assignImpropers", errorMessage, INFO);

    bool gotImproper;
    double dImproper;

    pMol->bImpropersAssigned = true;
    subMoleculeList = pMol->getSubMoleculeList();

    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];
      atomList = pSubMolecule->getAtomList();

      for (unsigned int k = 0; k < atomList.size(); k++) {
        gotImproper = false;
        pAtom1 = atomList[k];
        for (unsigned int l = 0; l < pAtom1->bondedAtoms.size(); l++) {
          pAtom2 = pAtom1->bondedAtoms[l];
          if (pAtom2->bondedAtoms.size() == 3) {
            for (unsigned int o = 0; o < 3; o++) {
              if (gotImproper) continue;
              pAtom3 = pAtom2->bondedAtoms[o];
              if (pAtom1 != pAtom3) {
                for (unsigned int p = 0; p < pAtom2->bondedAtoms.size(); p++) {
                  pAtom4 = pAtom2->bondedAtoms[p];
                  if ((pAtom1 != pAtom4) && (pAtom3 != pAtom4)){
                    gotImproper = pMol->hasImproper(pAtom1, pAtom2, pAtom3, pAtom4);
                    dImproper = torsion(*(pAtom1->getCoords()), *(pAtom2->getCoords()),
                                        *(pAtom3->getCoords()), *(pAtom4->getCoords()));
                    if (!gotImproper) {
                      pMol->addImproper(pAtom1, pAtom2, pAtom3, pAtom4, dImproper);
                      gotImproper = true;
                      continue;
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

// ============================================================
// Function : assignStdImpropers()
// ------------------------------------------------------------
// Assigns parameters to all impropers for Every Molecule in the Collection
// ============================================================
void connections::assignStdImpropers()
{
/*
#ifdef DEBUG
    std::cout << "   connections::assignStdImpropers" << std::endl;
#endif
    typedef std::map<int, Improper*>::iterator ImproperMapIterator;
    std::map<int, Improper*>      moleculeImproperMap;
    std::vector<improperParam*> improperParamList;

    moleculeList = pCollection->getMoleculeList();
    pParameters = pCollection->getParameters();

    for (unsigned int i = 0; i < moleculeList.size(); i++){
      pMolecule = moleculeList[i];
      moleculeImproperMap = pMolecule->getImproperMap();
      if (!moleculeImproperMap.empty()){
        // Loop over all impropers
        for (ImproperMapIterator a = moleculeImproperMap.begin(); a!=moleculeImproperMap.end();a++){
          pImproper = a->second;
          improperParamList = pParameters->getImproperParamList(
                                               pTorsion->atom1->getStdAtom()->type,
                                               pTorsion->atom2->getStdAtom()->type,
                                               pTorsion->atom3->getStdAtom()->type,
                                               pTorsion->atom4->getStdAtom()->type);
          if (!improperParamList.empty()) {
            pImproper->pImproperParamList = improperParamList;
            for (unsigned int k = 0; k < improperParamList.size(); k++ ){
               std::cout << improperParamList[k]->atomType1 << " "
                    << improperParamList[k]->atomType2 << " "
                    << improperParamList[k]->atomType3 << " "
                    << improperParamList[k]->atomType4 << " "
                    << improperParamList[k]->Nt << " "
                    << improperParamList[k]->Vn << " "
                    << improperParamList[k]->gamma << std::endl;
            }
          }
          else{
            std::cout << "No Improper Parameter for the improper containing: " <<
                      pImproper->atom1->getName() << " - " << pImproper->atom2->getName() << " - " <<
                      pImproper->atom3->getName() << " - " << pImproper->atom4->getName() << " - " <<
                      pImproper->atom1->getStdAtom()->type << " - " << pImproper->atom2->getStdAtom()->type <<
                      " - " <<
                      pImproper->atom3->getStdAtom()->type << " - " << pImproper->atom4->getStdAtom()->type
                      << std::endl;
          }
        }
      }
      else{
        std::cout << "Please assign impropers before assigning MM parameters for molecule: "
                  << pMolecule->getName() << std::endl;
      }
    }*/
}

// ============================================================
// Function : assignStdImpropers()
// ------------------------------------------------------------
// Assigns parameters to all impropers for Every Molecule in the Collection
// ============================================================
void connections::assignStdImpropers(molecule* pMol)
{
//#ifdef DEBUG
//    std::cout << "   connections::assignStdImpropers: molecule " << pMol->getName() << std::endl;
//#endif
}

} // MTKpp namespace

