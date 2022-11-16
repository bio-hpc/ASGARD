/*!
   \file atomTyper.cpp
   \brief Atom Types molecules
   \author Martin Peters

   Assigns standard fragments to each residue/fragment

   Assigns standard atoms to each atom

   $Date: 2010/03/29 20:42:27 $
   $Revision: 1.20 $

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

#include "atomTyper.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "ring.h"
#include "functionalize.h"
#include "utility.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"

#include "parameters.h"
#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : atomTyper()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
atomTyper::atomTyper()
{
    pCollection  = 0;
    pMolecule    = 0;
    pSubMolecule = 0;
    pAtom        = 0;

    pStdLibrary  = 0;
    pStdGroup    = 0;
    pStdFrag     = 0;
    pStdAtom     = 0;

    bPerceiveHistidines = 1;
    bPerceiveCysteines = 1;
}

// ============================================================
// Function : atomTyper()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
atomTyper::atomTyper(int b)
{
    pCollection  = 0;
    pMolecule    = 0;
    pSubMolecule = 0;
    pAtom        = 0;

    pStdLibrary  = 0;
    pStdGroup    = 0;
    pStdFrag     = 0;
    pStdAtom     = 0;

    bPerceiveHistidines = b;
    bPerceiveCysteines = b;
}

// ============================================================
// Function : ~atomTyper()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
atomTyper::~atomTyper() {}

// ============================================================
// Function : atomTypeByLib()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::atomTypeByLib(collection* pCol)
{
    pCollection = pCol;
    moleculeList = pCollection->getMoleculeList();
    pStdLibrary = pCollection->getStdLibrary();
    pParameters = pCollection->getParameters();

    if (bPerceiveHistidines) {
      this->perceiveHistidines();
    }

    if (bPerceiveCysteines) {
      this->perceiveCysteines();
    }

    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      pMolecule = moleculeList[i];
      this->atomTypeByLib(pMolecule);
    }
}

// ============================================================
// Function : atomTypeByLib()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::atomTypeByLib(molecule* pMolecule)
{
    std::string errorMessage = "\n";

    collection* pCollection = pMolecule->getParent();
    pStdLibrary = pCollection->getStdLibrary();
    pParameters = pCollection->getParameters();

    subMoleculeList = pMolecule->getSubMoleculeList();
    int nGotStdAtom = 0;

    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];
      if (!pSubMolecule->getStdFrag()) {
        // - First Residue
        if (j == 0) {
          pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName(), std::string("s"));

          if (!pStdFrag) {
            pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName(), std::string("l"));

            if (!pStdFrag) {
              pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName());
            }
          }
        }
        // - Last Residue
        else if (j == subMoleculeList.size()-1) {
          pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName(), std::string("e"));

          if (!pStdFrag) {
            pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName(), std::string("l"));

            if (!pStdFrag) {
              pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName());
            }
          }
        }
        // - All Other Residues
        else {
          pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName(), std::string("l"));

          if (!pStdFrag) {
            pStdFrag = pStdLibrary->getStdFrag(pSubMolecule->getName());
          }
        }
      }
      else {
        pStdFrag = pSubMolecule->getStdFrag();
      }

      if (pStdFrag) {
        errorMessage += " Residue: " + i2s(pSubMolecule->getSubMolId()) + " Assigned " 
                      + pSubMolecule->getName() + " (" + pStdFrag->getType() + ")\n";

        // - Keep a pointer to the stdFragment
        pSubMolecule->setStdFrag(pStdFrag);
        // - Atom type all atoms in submolecule
        atomList = pSubMolecule->getAtomList();
        for (unsigned int k = 0; k < atomList.size(); k++) {
          pAtom = atomList[k];
          pStdAtom = pStdFrag->getStdAtom(pAtom->getName());
          if (pStdAtom) {
            pAtom->setStdAtom(pStdAtom);
            pAtom->setName(pStdAtom->identity);
            nGotStdAtom++;

            atomType* pAtomType = pParameters->getAtomType(pStdAtom->type);
            if (!pAtomType) {
              std::string errMessage = "|" + pAtom->getName() + "|-" + pStdAtom->type + " Can't find atomType";
              errorLogger.throwError("atomTyper::atomTypeByLib", errMessage, 1);
            }

            /*
              atom kind
            */
            pAtom->setType(pStdAtom->kind);

            /*
              atom hybridization definitions
              - 0 = undefined
              - 1 = s
              - 2 = sp
              - 3 = sp2
              - 4 = sp3
              - 5 = sp3d
              - 6 = sp3d2
              - 7 = other
            */
            if (pAtomType->hybridization == "s") {
              pAtom->setHybridization(1);
            }
            else if (pAtomType->hybridization == "sp") {
              pAtom->setHybridization(2);
            }
            else if (pAtomType->hybridization == "sp2") {
              pAtom->setHybridization(3);
            }
            else if (pAtomType->hybridization == "sp3") {
              pAtom->setHybridization(4);
            }
            else if (pAtomType->hybridization == "sp3d") {
              pAtom->setHybridization(5);
            }
            else if (pAtomType->hybridization == "sp3d2") {
              pAtom->setHybridization(6);
            }
            else {
              pAtom->setHybridization(0);
            }
          }
          else {
            std::stringstream smolN;
            smolN << pSubMolecule->getSubMolId();
            //std::cout << j << " " << subMoleculeList.size() << std::endl;
            std::string errMessage = "\nUnable to assign stdAtom: Unknown atom: |" + pAtom->getName()
                      + "|\n in standard residue: " + pSubMolecule->getName() + "("  + pStdFrag->getType() + ")-" +
                      smolN.str() + "\n";

            std::vector<stdAtom*> s_atomList = pStdFrag->getStdAtomList();
            for (unsigned int d = 0; d < s_atomList.size(); d++) {
              errMessage += " AVAILABLE ATOM NAME : |" + s_atomList[d]->identity + "|\n";
            }
            errorLogger.throwError("atomTyper::atomTypeByLib", errMessage, 1);
            //exit(1);
            std::stringstream ss;
            ss << " atomTyper::atomTypeByLib" << errMessage << std::endl;
            std::cout << ss.str();
            throw MTKException(ss.str());
          }
        }
      }
      else {
        std::string errMessage = "Unable to assign stdFrag: unknown residue |" +
                    pSubMolecule->getName() + "|";
        errorLogger.throwError("atomTyper::atomTypeByLib", errMessage, WARNING);
      }
      pStdFrag = 0;
    }

    if (pMolecule->getNumAtoms() == nGotStdAtom) {
      pMolecule->bMMAtomTypesAssigned = true;
      pMolecule->bAtomHybridizationsAssigned = true;
    }

    errorLogger.throwError("atomTyper::atomTypeByLib", errorMessage, INFO);
}

// ============================================================
// Function : assignStdProperties()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::assignStdProperties(molecule* pMolecule)
{
    if (!pStdLibrary or !pParameters) {
      pCollection = pMolecule->getParent();
      pStdLibrary = pCollection->getStdLibrary();
      pParameters = pCollection->getParameters();
    }

    if (pMolecule->bMMAtomTypesAssigned and pMolecule->bAtomHybridizationsAssigned) {
      if (pMolecule->getNumSubMolecules() == 1) {
        this->assignStdRings(pMolecule);
        this->assignStdFunctionalGroups(pMolecule);
        this->assignStdHydrophobicGroups(pMolecule);
      }
    }
}

// ============================================================
// Function : atomTypeByCon()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::atomTypeByCon(collection* pCollection)
{
    std::cout << pCollection->getName() << std::endl;
}

// ============================================================
// Function : perceiveHistidines()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::perceiveHistidines()
{
    std::string errMessage = " Perceiving Histidine Residue Types \n";
    // Find all histidines
    std::vector<submolecule*> histidines;
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      std::vector<submolecule*> subMolList = moleculeList[i]->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        if (subMolList[j]->getName() == "HIS") {
          if (subMolList[j]->getStdFrag()) continue;
          histidines.push_back(subMolList[j]);
        }
      }
    }

    for (unsigned int i = 0; i < histidines.size(); i++) {
      std::stringstream smolN;
      smolN << histidines[i]->getSubMolId();
      errMessage += "  " + histidines[i]->getName() + ":" + smolN.str() + "\n";
      std::vector<atom*> atomList = histidines[i]->getAtomList();

      if (histidines[i]->getAtom("1H  ") or histidines[i]->getAtom(" H1 ")) {
        histidines[i]->setName("HIP");
        errMessage += "        Setting to HIP \n";
        continue;
      }

      bool bHIE = false;
      bool bHID = false;

      for (unsigned int j = 0; j < atomList.size(); j++) {
        if (atomList[j]->getName() == " NE2" or atomList[j]->getName() == " ND1") {

          int numHeavyNeighbors = pCollection->getNumHeavyNeighbors(atomList[j], 2.35);
          if (numHeavyNeighbors > 0) {
            if (atomList[j]->getName() == " NE2") {
              histidines[i]->setName("HID");
              errMessage += "        Setting to HID \n";
            }
            else {
              histidines[i]->setName("HIE");
              errMessage += "        Setting to HIE \n";
            }
          }

          int numNeighbors = 0;
          for (unsigned int k = 0; k < atomList.size(); k++) {
            if (k == j) continue;
            if (atomList[j]->getCoords()->dist((*atomList[k]->getCoords())) < 1.5) {
              numNeighbors++;
              //std::cout << "    " << atomList[j]->getName() << " " << atomList[k]->getName() << std::endl;
            }
          }
          if (numNeighbors == 3) {
            if (atomList[j]->getName() == " NE2") {
              bHIE = true;
            }
            else {
              bHID = true;
            }
          }
/*
          std::cout << histidines[i]->getName() <<  ":" << smolN.str() << " " << numHeavyNeighbors 
          << " " << numNeighbors << " "
          << atomList[j]->getNumBonds()
          << "\n";
*/
        }
      }

      if (bHIE and bHID) {
        histidines[i]->setName("HIP");
        errMessage += "        Setting to HIP \n";
      }

      // ND1 - HD1 --> HID
      // NE2 - HE2 --> HIE
      if (histidines[i]->getName() == "HIS" and bHIE) {
        histidines[i]->setName("HIE");
        errMessage += "        Setting to HIE \n";
      }
      else if (histidines[i]->getName() == "HIS" and bHID) {
        histidines[i]->setName("HID");
        errMessage += "        Setting to HID \n";
      }
    }
    errorLogger.throwError("atomTyper::perceiveHistidines", errMessage, INFO);
}

// ============================================================
// Function : perceiveCysteines()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::perceiveCysteines()
{
    std::string errMessage = " Perceiving Cysteines Residue Types \n";
    // Find all Cysteines
    std::vector<submolecule*> cysteines;
    for (unsigned int i = 0; i < moleculeList.size(); i++) {
      std::vector<submolecule*> subMolList = moleculeList[i]->getSubMoleculeList();
      for (unsigned int j = 0; j < subMolList.size(); j++) {
        if (subMolList[j]->getName() == "CYS") {
          if (subMolList[j]->getStdFrag()) continue;
          cysteines.push_back(subMolList[j]);
        }
      }
    }

    for (unsigned int i = 0; i < cysteines.size(); i++) {
      std::stringstream smolN;
      smolN << cysteines[i]->getSubMolId();
      errMessage += "  " + cysteines[i]->getName() + ":" + smolN.str() + "\n";
      std::vector<atom*> atomList = cysteines[i]->getAtomList();

      if (cysteines[i]->getAtom("HG  ") or cysteines[i]->getAtom(" HG ")) {
        cysteines[i]->setName("CYS");
        errMessage += "        Setting to CYS \n";
        continue;
      }

      bool bCYM = false;

      for (unsigned int j = 0; j < atomList.size(); j++) {
        if (atomList[j]->getName() == " SG " or atomList[j]->getName() == "SG  ") {
          int numHeavyNeighbors = pCollection->getNumHeavyNeighbors(atomList[j], 2.90);
          if (numHeavyNeighbors) {
            bCYM = true;
          }
          //std::cout << cysteines[i]->getName() <<  ":" << smolN.str() << " " << numHeavyNeighbors << "\n";
        }
      }
      if (bCYM) {
        cysteines[i]->setName("CYM");
        errMessage += "        Setting to CYM \n";
      }
    }
    errorLogger.throwError("atomTyper::perceiveCysteines ", errMessage, INFO);
}

// ============================================================
// Function : assignStdRings()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::assignStdRings(molecule* pMol)
{
    errorLogger.throwError("atomTyper::assignStdRings", pMol->getName(), INFO);

    subMoleculeList = pMol->getSubMoleculeList();
    pSubMolecule = subMoleculeList[0];
    pStdFrag = pSubMolecule->getStdFrag();

    rings* pRings = new rings(pMol);
    std::vector<stdRing*> ringList = pStdFrag->getStdRingList();
    for (unsigned int i = 0; i < ringList.size(); i++) {
      std::vector<atom*> cRing;
      for (unsigned int j = 0; j < ringList[i]->atoms.size(); j++) {
        pStdAtom = pStdFrag->getStdAtom(ringList[i]->atoms[j]);
        if (!pStdAtom) {
          std::string errMessage = " Can't find stdAtom " + ringList[i]->atoms[j];
          errorLogger.throwError("atomTyper::assignStdRings", errMessage, 1);
          //exit(1);
          std::stringstream ss;
          ss << "atomTyper::assignStdRings"<< errMessage;
          throw MTKException(ss.str());
        }
        pAtom = pMol->getAtom(pStdAtom);
        if (pAtom) {
          cRing.push_back(pAtom);
        }
        else {
          errorLogger.throwError("atomTyper::assignStdRings", " Failed to find atom ", MTK_ERROR);
          //exit(1);
          std::stringstream ss;
          ss << "atomTyper::assignStdRings"<< " Failed to find atom ";
          throw MTKException(ss.str());
        }
      }
      pMol->addRing(cRing);
      ring* pRing = pMol->getLastAddedRing();
      if (!pRing) {
        errorLogger.throwError("atomTyper::assignStdRings", " Failed to find ring ", MTK_ERROR);
        //exit(1);
        std::stringstream ss;
        ss << "atomTyper::assignStdRings"<< " Failed to find ring ";
        throw MTKException(ss.str());
      }
      pRing->size      = ringList[i]->size;
      pRing->planar    = ringList[i]->planar;
      pRing->aromatic  = ringList[i]->aromatic;
      pRing->hetero    = ringList[i]->hetero;
      pRing->nHetero   = ringList[i]->nHetero;
      pRing->nNitrogen = ringList[i]->nNitrogen;
      pRing->nOxygen   = ringList[i]->nOxygen;
      pRing->nSulfur   = ringList[i]->nSulfur;
      pRings->getPlaneNormal(pRing);
    }
}

// ============================================================
// Function : assignStdFunctionalGroups()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::assignStdFunctionalGroups(molecule* pMol)
{
    errorLogger.throwError("atomTyper::assignStdFunctionalGroups", pMol->getName(), INFO);

    subMoleculeList = pMol->getSubMoleculeList();
    if (subMoleculeList.size() < 1) {
      errorLogger.throwError("atomTyper::assignStdFunctionalGroups", " No residues in molecule ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "atomTyper::assignStdFunctionalGroups"<< " No residues in molecule ";
      throw MTKException(ss.str());
    }

    pSubMolecule = subMoleculeList[0];
    if (!pSubMolecule) {
      errorLogger.throwError("atomTyper::assignStdFunctionalGroups", " No residues in molecule ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "atomTyper::assignStdFunctionalGroups"<< " No residues in molecule ";
      throw MTKException(ss.str());
    }

    pStdFrag = pSubMolecule->getStdFrag();
    if (!pStdFrag) {
      errorLogger.throwError("atomTyper::assignStdFunctionalGroups", " Assigned standard residues ", MTK_ERROR);
      std::stringstream ss;
      ss << "atomTyper::assignStdFunctionalGroups"<< " Assigned standard residues ";
      throw MTKException(ss.str());
    }

    std::vector<stdFuncGroup*> funcGroupList = pStdFrag->getStdFuncGroupList();
    for (unsigned int i = 0; i < funcGroupList.size(); i++) {

      // name:terminal/CAR
      std::string gname = funcGroupList[i]->groupName;
      std::string fname = funcGroupList[i]->fragName;
      stdGroup* pStdGrp = pStdLibrary->getStdGroup(gname);
      if (!pStdGrp) {
        std::string errMessage = " Can't find stdGroup " + gname;
        errorLogger.throwError("atomTyper::assignStdFunctionalGroups", errMessage, MTK_ERROR);
        //exit(1);
        std::stringstream ss;
        ss << "atomTyper::assignStdFunctionalGroups"<< " Can't find stdGroup ";
        throw MTKException(ss.str());
      }

      stdFrag* pFuncFrag = pStdGrp->getStdFrag(fname);
      if (!pFuncFrag) {
        std::string errMessage = " Can't find stdFrag " + fname;
        errorLogger.throwError("atomTyper::assignStdFunctionalGroups", errMessage, MTK_ERROR);
        //exit(1);
        std::stringstream ss;
        ss << "atomTyper::assignStdFunctionalGroups"<< " Can't find stdFrag ";
        throw MTKException(ss.str());
      }

      // <funcGroup atoms="2 1 3 " fragment="CAR" group="terminal"/>
      funcGroup* fgroup = new funcGroup();
      for (unsigned int j = 0; j < funcGroupList[i]->atoms.size(); j++) {

        pStdAtom = pFuncFrag->getStdAtom(j+1);
        if (!pStdAtom) {
          std::string errMessage = " Can't find stdAtom " + funcGroupList[i]->atoms[j];
          errorLogger.throwError("atomTyper::assignStdFunctionalGroups", errMessage, MTK_ERROR);
          //exit(1);
          std::stringstream ss;
          ss << "atomTyper::assignStdFunctionalGroups" << errMessage;
          throw MTKException(ss.str());
        }

        pAtom = pMol->getAtom(funcGroupList[i]->atoms[j], 1, 0);
        if (!pAtom) {
          std::string errMessage = " Can't find atom |" + pStdAtom->identity + "|";
          errorLogger.throwError("atomTyper::assignStdFunctionalGroups", errMessage, MTK_ERROR);
          //exit(1);
          std::stringstream ss;
          ss << "atomTyper::assignStdFunctionalGroups"<< errMessage;
          throw MTKException(ss.str());
        }

        fgroup->atomMap[pStdAtom] = pAtom;
      }

      fgroup->pStdFrag = pFuncFrag;
      pMol->addFunctionalGroup(fgroup);
    }
}

// ============================================================
// Function : assignStdHydrophobicGroups()
// ------------------------------------------------------------
//
// ============================================================
void atomTyper::assignStdHydrophobicGroups(molecule* pMol)
{
    errorLogger.throwError("atomTyper::assignStdHydrophobicGroups", pMol->getName(), INFO);

    subMoleculeList = pMol->getSubMoleculeList();
    pSubMolecule = subMoleculeList[0];
    pStdFrag = pSubMolecule->getStdFrag();

    std::vector<stdFeature*> featureList = pStdFrag->getStdFeatureList();
    for (unsigned int i = 0; i < featureList.size(); i++) {
      if (featureList[i]->name == "HPB") {
        std::vector<atom*> fatoms;
        for (unsigned int j = 0; j < featureList[i]->atoms.size(); j++) {

          pStdAtom = pStdFrag->getStdAtom(featureList[i]->atoms[j]);
          if (!pStdAtom) {
            std::string errMessage = " Can't find stdAtom " + featureList[i]->atoms[j];
            errorLogger.throwError("atomTyper::assignStdHydrophobicGroups", errMessage, MTK_ERROR);
            //exit(1);
            std::stringstream ss;
            ss << "atomTyper::assignStdHydrophobicGroups"<< errMessage;
            throw MTKException(ss.str());
          }

          pAtom = pMol->getAtom(pStdAtom);
          if (!pAtom) {
            std::string errMessage = " Can't find atom " + pStdAtom->identity;
            errorLogger.throwError("atomTyper::assignStdHydrophobicGroups", errMessage, MTK_ERROR);
            //exit(1);
            std::stringstream ss;
            ss << "atomTyper::assignStdHydrophobicGroups"<< errMessage;
            throw MTKException(ss.str());
          }
          fatoms.push_back(pAtom);
        }
        pMol->addHydrophobicGroup(fatoms);
      }
    }
}


} // MTKpp namespace
