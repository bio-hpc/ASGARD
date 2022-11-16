/*!
   \file proProtonate.cpp
   \brief Protonates proteins
   \author Martin Peters
   \author Andrew Wollacott

   $Date: 2010/03/29 20:44:27 $
   $Revision: 1.13 $

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

#include "proProtonate.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"
#include "ring.h"
#include "angle.h"
#include "torsion.h"
#include "improper.h"
#include "connections.h"
#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "parameters.h"
#include "utility.h"

#include "Log/errorHandler.h"

// Utils
#include "Utils/idObject.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

#include <math.h>

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : proProtonate()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
proProtonate::proProtonate()
{
    dummyAtom1 = 0;
    dummyAtom2 = 0;
    dummyAtom3 = 0;
}

// ============================================================
// Function : ~proProtonate()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
proProtonate::~proProtonate()
{
    delete dummyAtom1;
    delete dummyAtom2;
    delete dummyAtom3;
}

// ============================================================
// Function : buildMissingAtoms()
// ------------------------------------------------------------
//  Build Dummy atoms needed for the first residue
// ============================================================
void proProtonate::buildMissingAtoms(submolecule* pSubMolecule, std::vector<stdAtom*> missingAtoms,
                                     std::vector<atom*> first3Atoms)
{
    std::string errMessage = "\n";

    errMessage += pSubMolecule->getName() + " First 3 atom = " + i2s(first3Atoms[0]->getFileID())
               + "-" + i2s(first3Atoms[1]->getFileID()) + "-"
               + i2s(first3Atoms[2]->getFileID()) + "\n";

    pMol = pSubMolecule->getParent();
    pCol = pMol->getParent();
    pParam = pCol->getParameters();
    if (!pMol or !pCol or !pParam) {
      std::cout << " Error in proProtonate " << std::endl;
      //exit(1);
      throw MTKException(" Error in proProtonate ");
    }

    std::vector<atom*> tmpAtoms = first3Atoms;
    int c = missingAtoms.size()-1;

    for (unsigned int i = 0; i < missingAtoms.size(); i++) {
      pAtom1 = pSubMolecule->addAtom();
      std::string atSymbol = pParam->getAtomTypeSymbol(missingAtoms[c]->type);
      pAtom1->setElement(pCol->pElements->getElement(atSymbol));

      if (!pAtom1->getElement()) {
        std::cout << " Error ... exiting " << std::endl;
        //exit(0);
        throw MTKException(" Error in proProtonate ");
      }

      //pAtom1->setFileID(pMol->getNumAtoms());
      int ifileID = pMol->getMaxFileID();
      pAtom1->setFileID(ifileID+1);
      pMol->setMaxFileID(ifileID+1);

      pAtom1->setCoords(0,0,0);
      pAtom1->setStdAtom(missingAtoms[c]);
      pAtom1->setName(missingAtoms[c]->identity);

      double bond12 = tmpAtoms[0]->getStdAtom()->bondLength;
      double angle13 = tmpAtoms[1]->getStdAtom()->bondAngle/RAD2DEG;
      double torsion14 = -(tmpAtoms[2]->getStdAtom()->bondTorsion/RAD2DEG);
      errMessage += " Adding" + missingAtoms[c]->identity + " " +  d2s(bond12) + " "
                 + d2s(angle13) + " " + d2s(torsion14) + "\n";
      buildCoord(*(pAtom1->getCoords()), *(tmpAtoms[0]->getCoords()),
                 *(tmpAtoms[1]->getCoords()), *(tmpAtoms[2]->getCoords()),
                 bond12, angle13, torsion14);

      pBond = pMol->addBond(pAtom1, tmpAtoms[0], 1, 0, 1, 0.0);
      pAtom1->addBondedAtom(tmpAtoms[0]);
      tmpAtoms[0]->addBondedAtom(pAtom1);
      if (atSymbol == "H") {
        pAtom1->setType(1);
      }
      pAtom1->setValence(pAtom1->getValence()+1);

      tmpAtoms[0] = pAtom1;
      tmpAtoms[1] = tmpAtoms[0];
      tmpAtoms[2] = tmpAtoms[1];
      c--;
    }
    errorLogger.throwError("proProtonate::buildMissingAtoms", errMessage, INFO);
}

// ============================================================
// Function : buildDummyAtoms()
// ------------------------------------------------------------
//  Build Dummy atoms needed for the first residue
// ============================================================
void proProtonate::buildDummyAtoms(std::vector<atom*> first3Atoms)
{
    std::string errMessage = " First 3 atom = " + i2s(first3Atoms[0]->getFileID())
              + "-" + i2s(first3Atoms[1]->getFileID()) + "-" 
              + i2s(first3Atoms[2]->getFileID()) + "\n";

    delete dummyAtom1;
    delete dummyAtom2;
    delete dummyAtom3;

    dummyAtom1 = new atom();
    dummyAtom2 = new atom();
    dummyAtom3 = new atom();

    prev3Atoms.push_back(dummyAtom1);
    prev3Atoms.push_back(dummyAtom2);
    prev3Atoms.push_back(dummyAtom3);

    dummyAtom1->setCoords(0,0,0);
    dummyAtom2->setCoords(0,0,0);
    dummyAtom3->setCoords(0,0,0);

    dummyAtom1->setName("XXXX");
    dummyAtom2->setName("XXXX");
    dummyAtom3->setName("XXXX");

    stdAtom* pStdAt1 = first3Atoms[0]->getStdAtom();
    stdAtom* pStdAt2 = first3Atoms[1]->getStdAtom();
    stdAtom* pStdAt3 = first3Atoms[2]->getStdAtom();

/*
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000
   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000   0.00000
   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000   0.00000
   4  N     N3    M    3   2   1     1.335   116.600   180.000   0.14140
   8  CA    CT    M    4   3   2     1.449   121.900   180.000   0.09620
  14  C     C     M    8   4   3     1.522   111.100   180.000   0.61630
*/
    double bond12 = pStdAt1->bondLength;
    double angle13 = pStdAt2->bondAngle/RAD2DEG;
    double torsion14 = -(pStdAt3->bondTorsion/RAD2DEG);

    buildCoord(*(dummyAtom3->getCoords()), *(first3Atoms[0]->getCoords()),
               *(first3Atoms[1]->getCoords()), *(first3Atoms[2]->getCoords()),
               bond12, angle13, torsion14);

    angle13 = pStdAt1->bondAngle/RAD2DEG;
    torsion14 = -(pStdAt2->bondTorsion/RAD2DEG);

    buildCoord(*(dummyAtom2->getCoords()), *(dummyAtom3->getCoords()),
               *(first3Atoms[0]->getCoords()), *(first3Atoms[1]->getCoords()),
               1.522, angle13, torsion14);

    torsion14 = -(pStdAt1->bondTorsion/RAD2DEG);

    buildCoord(*(dummyAtom1->getCoords()), *(dummyAtom2->getCoords()),
               *(dummyAtom3->getCoords()), *(first3Atoms[0]->getCoords()),
               1.449, 111.1, torsion14);
/*
#ifdef DEBUG
    std::cout << dummyAtom1->getCoords() << std::endl;
    std::cout << dummyAtom2->getCoords() << std::endl;
    std::cout << dummyAtom3->getCoords() << std::endl;
#endif
*/
    errorLogger.throwError("proProtonate::buildDummyAtoms", errMessage, INFO);
}

// ============================================================
// Function : getBestTorsion()
// ------------------------------------------------------------
//
// ============================================================
std::vector<double> proProtonate::getBestTorsions(stdFrag* pStdFrag,
                    stdAtom* stdAt1, atom* at2, atom* at3, atom* at4)
{
    std::string errMessage = "";
    std::vector<double> newTorsions;
    if (at3->getName() == "XXXX" and at4->getName() == "XXXX") {
      newTorsions.push_back(stdAt1->bondTorsion);
      return newTorsions;
    }
    if (at4->getName() == "XXXX") {
      newTorsions.push_back(stdAt1->bondTorsion);
      return newTorsions;
    }

    pStdAtom2 = at2->getStdAtom();
    int nStdBonds = 2;
    if (pStdAtom2) {
      nStdBonds = pStdFrag->numStdBonds(pStdAtom2);
    }
    else {
      newTorsions.push_back(stdAt1->bondTorsion);
      return newTorsions;
    }

    std::vector<atom*> bondedAtoms = at2->getBondedAtoms();
    atom* bondedAtom = 0;

    double defaultTor = 0.0;
    if (nStdBonds == 4) defaultTor = 120.0;
    if (nStdBonds == 3) defaultTor = 180.0;
    if (nStdBonds == 2) defaultTor = stdAt1->bondTorsion;

    errMessage = i2s(at2->getFileID()) + "-" + i2s(at3->getFileID()) + "-" + i2s(at4->getFileID()) + "  " +
                 pStdAtom2->identity + " has " + i2s(nStdBonds) + " ";
	if (nStdBonds == 1) {
		errMessage += "bond";
	} else {
		errMessage += "bonds";
	}
	errMessage += ", default torsion = " + d2s(defaultTor) + "\n";

    vector3d* a;
    vector3d* b;
    vector3d* c;
    vector3d* d;
    b = at2->getCoords();
    c = at3->getCoords();
    d = at4->getCoords();

    double tor = 0.0;
    double prevTor = 0.0;
    double diffTor = 0.0;
    double hTor = -360.1;
    double retTor = defaultTor;

    for (unsigned int i = 0; i < bondedAtoms.size(); i++) {
      bondedAtom = bondedAtoms[i];
      if (bondedAtom != at3) {
        a = bondedAtom->getCoords();
        //double dTorsion = torsion((*a), (*b), (*c), (*d));
        double dTorsion2 = torsion2((*a), (*b), (*c), (*d)); // -PI to PI
        if (dTorsion2 < 0) {
            dTorsion2 = 2*PI + dTorsion2;
        }

        tor = RAD2DEG * dTorsion2;
        if (tor > hTor) hTor = tor;
        diffTor = fabs(tor - prevTor);
        prevTor = tor;
      }
    }

    if (nStdBonds == 4) {
      diffTor /= 120.0;
      if ((diffTor > 1.7) && (diffTor < 2.1)) { // andrews values: 1.9 - 2.1
        retTor = hTor - 120.0;
      }
      else {
        retTor = hTor + 120.0;
      }
    }

    if (nStdBonds == 3) {
      retTor = hTor + 180.0;
    }
    newTorsions.push_back(retTor);

    errorLogger.throwError("proProtonate::getBestTorsions", errMessage, INFO);
    return newTorsions;
}

// ============================================================
// Function : addHydrogens()
// ------------------------------------------------------------
//  Add Hydrogens in a molecule using libraries
// ============================================================
void proProtonate::addHydrogens(submolecule* pSubMolecule, stdFrag* pStdFrag)
{
    std::string errMessage = "\nUsing libraries for " + pStdFrag->getSymbol() + "\n";

    molecule* cMolecule = pSubMolecule->getParent();
    pCol = cMolecule->getParent();
    pParam = pCol->getParameters();

    std::vector<atom*> atomList = pSubMolecule->getAtomList();
    std::vector<stdAtom*> stdAtomList = pStdFrag->getStdAtomList();

    int nPrevChainAts = prev3Atoms.size();
    pAtom1 = 0; pAtom2 = 0; pAtom3 = 0; pAtom4 = 0;
    pStdAtom1 = 0; pStdAtom2 = 0; pStdAtom3 = 0; pStdAtom4 = 0;

    // check if all atoms are present, update prev3Atoms list and return
    if (atomList.size() == stdAtomList.size()) {
      for (unsigned int i = 0; i < stdAtomList.size(); i++) {
        pStdAtom1 = stdAtomList[i];
        pAtom1 = pSubMolecule->getAtom(pStdAtom1);
        if (pAtom1 and (pStdAtom1->chain == "M")) {
          prev3Atoms.push_back(pAtom1);
        }
      }
      return;
    }

    int bond12 = 0;
    int bond13 = 0;
    int bond14 = 0;

    double distance = 0.0;
    double angle = 0.0;
    std::vector<double> torsions;

    if (pStdFrag) {
      // make sure all atoms are in the stdFrag
      bool fragOk = true;
      for (unsigned int i = 0; i < atomList.size(); i++) {
        if (!pStdFrag->hasStdAtom(atomList[i]->getStdAtom()->identity)) {
          fragOk = false;
          break;
        }
      }
      if (fragOk) {
        // Add atoms
        for (unsigned int i = 0; i < stdAtomList.size(); i++) {
          pStdAtom1 = stdAtomList[i];
          pAtom1 = pSubMolecule->getAtom(pStdAtom1);

          if (pAtom1 and (pStdAtom1->chain == "M")) {
            prev3Atoms.push_back(pAtom1);
          }

          bond12 = pStdAtom1->bond12;
          if (bond12 > 0) {
            pStdAtom2 = pStdFrag->getStdAtom(bond12);
            pAtom2 = pSubMolecule->getAtom(pStdAtom2);
          }
          else {
            if (std::abs(bond12) == 1) {
              pAtom2 = prev3Atoms[nPrevChainAts-1];
            }
            if (std::abs(bond12) == 2) {
              pAtom2 = prev3Atoms[nPrevChainAts-2];
            }
            if (std::abs(bond12) == 3) {
              pAtom2 = prev3Atoms[nPrevChainAts-3];
            }
          }

          bond13 = pStdAtom1->bond13;
          if (bond13 > 0) {
            pStdAtom3 = pStdFrag->getStdAtom(bond13);
            pAtom3 = pSubMolecule->getAtom(pStdAtom3);
          }
          else {
            if (std::abs(bond13) == 1) {
              pAtom3 = prev3Atoms[nPrevChainAts-1];
            }
            if (std::abs(bond13) == 2) {
              pAtom3 = prev3Atoms[nPrevChainAts-2];
            }
            if (std::abs(bond13) == 3) {
              pAtom3 = prev3Atoms[nPrevChainAts-3];
            }
          }

          bond14 = pStdAtom1->bond14;
          if (bond14 > 0) {
            pStdAtom4 = pStdFrag->getStdAtom(bond14);
            pAtom4 = pSubMolecule->getAtom(pStdAtom4);
          }
          else {
            if (std::abs(bond14) == 1) {
              pAtom4 = prev3Atoms[nPrevChainAts-1];
            }
            if (std::abs(bond14) == 2) {
              pAtom4 = prev3Atoms[nPrevChainAts-2];
            }
            if (std::abs(bond14) == 3) {
              pAtom4 = prev3Atoms[nPrevChainAts-3];
            }
          }

          if (!pAtom1) {
            if (pParam == 0) {
              std::cout << " parameter is null in proProtonate " << std::endl;
              return;
            }
            std::string atSymbol = pParam->getAtomTypeSymbol(pStdAtom1->type);

            if (atSymbol == "") {
              std::cout << " pStdAtom symbol is null " << std::endl;
              return;
            }

            if (!pAtom2) {
              std::cout << " pAtom2 == 0 " << std::endl;
              return;
            }
            if (!pAtom3) {
              std::cout << " pAtom3 == 0 " << std::endl;
              return;
            }
            if (!pAtom4) {
              std::cout << pStdFrag->getSymbol() << "-" << pSubMolecule->getSubMolId()
                        << " --> " << pStdAtom1->identity << std::endl;
              std::cout << " pAtom4 == 0 " << std::endl;
              //exit(0);
              throw MTKException(" pAtom4 == 0 ");
            }

            pAtom1 = pSubMolecule->addAtom();
            pAtom1->setElement(pCol->pElements->getElement(atSymbol));

            //pAtom1->setFileID(cMolecule->getNumAtoms());
            int ifileID = cMolecule->getMaxFileID();
            pAtom1->setFileID(ifileID+1);
            cMolecule->setMaxFileID(ifileID+1);

            pAtom1->setCoords(0,0,0);

            pAtom1->setStdAtom(pStdAtom1);
            pAtom1->setName(pStdAtom1->identity);

            distance = pStdAtom1->bondLength;
            angle = pStdAtom1->bondAngle/RAD2DEG;

            torsions = this->getBestTorsions(pStdFrag, pStdAtom1, pAtom2, pAtom3, pAtom4);

            errMessage += pAtom1->getName() + "-" + pAtom2->getName()
                        + "-" + pAtom3->getName() + "-" + pAtom4->getName()
                        + " distance = " + d2s(distance) + " angle = " + d2s(angle)
                        + " torsion = " + d2s(torsions[0]) + "\n";

            buildCoord(*(pAtom1->getCoords()), *(pAtom2->getCoords()),
                       *(pAtom3->getCoords()), *(pAtom4->getCoords()),
                       distance, angle, torsions[0]/RAD2DEG);

            pBond = cMolecule->addBond(pAtom1, pAtom2, 1, 0, 1, 0.0);
            pAtom1->addBondedAtom(pAtom2);
            pAtom2->addBondedAtom(pAtom1);
            if (atSymbol == "H") {
              pAtom1->setType(1);
            }
            pAtom1->setValence(pAtom1->getValence()+1);
          }
        }
      }
    }
    errorLogger.throwError("proProtonate::addHydrogens", errMessage, INFO);
}

// ============================================================
// Function : optimizePolarHs()
// ------------------------------------------------------------
//
// ============================================================
void proProtonate::optimizePolarHs()
{
    std::string errMessage = "\n" + pMol->getName() + "\n";
 
    if (!pMol) {
      std::cout << " pMol is zero " << std::endl;
      //exit(0);
      throw MTKException(" pMol is zero ");
    }

    // Get donors (O, N attach to H by a polar bond)
    std::map<int, Bond*> bondMap =  pMol->getBondMap();
    typedef std::map<int, Bond*>::iterator BondMapIterator;
    std::vector<atom*> polarHs;
    std::vector<atom*> bondAtoms;
    std::vector<atom*> angleAtoms;

    std::vector<std::vector<atom*> > polarTorsions;

    if (!bondMap.empty()) {
      for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
        pBond = b->second;
        if (!pBond) {
          std::cout << " pBond is zero " << std::endl;
          //exit(0);
          throw MTKException(" pBond is zero ");
        }

        submolecule* pSubMol1 = pBond->atom1->getParent();
        if (!pSubMol1) {
          std::cout << " pSubMol1 is zero " << std::endl;
          //exit(0);
          throw MTKException(" pSubMol1 is zero ");
        }

        submolecule* pSubMol2 = pBond->atom2->getParent();
        if (!pSubMol2) {
          std::cout << " pSubMol2 is zero " << std::endl;
          //exit(0);
          throw MTKException(" pSubMol2 is zero ");
        }
        if (pSubMol1 != pSubMol2) {
          continue;
        }

        stdFrag* pStdFrag = pSubMol1->getStdFrag();
        if (!pStdFrag) {
          std::cout << " pStdFrag is zero " << std::endl;
          //exit(0);
          throw MTKException(" pStdFrag is zero ");
        }

        stdBond* pStdBd = pStdFrag->getStdBond(pBond->atom1->getName(), pBond->atom2->getName());

        if (pStdBd) {
          if (pStdBd->kind == 1) { // polar bond
            if (pBond->atom1->getElementSymbol() == "H") {
              polarHs.push_back(pBond->atom1);
              //std::cout << " polarH found: " << pBond->atom1->getFileID() << " " << pBond->atom1->getName() << std::endl;
            }
            else if (pBond->atom2->getElementSymbol() == "H") {
              polarHs.push_back(pBond->atom2);
              //std::cout << " polarH found: " << pBond->atom1->getFileID() << " " << pBond->atom1->getName() << std::endl;
            }
          }
        }
      }
    }
    unsigned int nPolarHs = polarHs.size();
    if (nPolarHs < 1) {
      return;
    }
    std::map<ULONG_KIND, Torsion*> torsionMap =  pMol->getTorsionMap();
    typedef std::map<ULONG_KIND, Torsion*>::iterator TorsionMapIterator;

    errMessage += "Number of polar Hs = " + i2s(nPolarHs) + "\n";

    std::vector<std::vector<Torsion*> > polarHsTorsions;
    Torsion* pTorsion = 0;
    torsionParam* pTorsionParam;
    std::vector<torsionParam*> torsionParamList;
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    if (!torsionMap.empty()) {
      for (unsigned int i = 0; i < nPolarHs; i++) {
        std::vector<Torsion*> polarHsTorsion;
        for (TorsionMapIterator t = torsionMap.begin(); t != torsionMap.end(); t++) {
          pTorsion = t->second;
          torsionParamList = pTorsion->pTorsionParamList;
          if (torsionParamList.empty()) {
#ifdef DEBUG
            std::cout << " proProtonate: "
                      << pTorsion->atom1->getFileID() << " " << pTorsion->atom2->getFileID() << " "
                      << pTorsion->atom3->getFileID() << " " << pTorsion->atom4->getFileID() << std::endl;
            std::cout << pTorsion->atom1->getName() << " " << pTorsion->atom2->getName() << " "
                      << pTorsion->atom3->getName() << " " << pTorsion->atom4->getName()
                      << " torsion param list is empty "
                      << std::endl;
#endif
          }

          if (pTorsion->atom1 == polarHs[i]) {
            polarHsTorsion.push_back(pTorsion);
          }
          else if (pTorsion->atom4 == polarHs[i]) {
            polarHsTorsion.push_back(pTorsion);
          }
        }
        polarHsTorsions.push_back(polarHsTorsion);
      }
    }

    errMessage += " Polar Torsions: \n";
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      errMessage += i2s(polarHs[i]->getFileID()) + ":" + polarHs[i]->getName() + "\n";
      for (unsigned int j = 0; j < polarHsTorsions[i].size(); j++) {
        errMessage += "    " + i2s(polarHsTorsions[i][j]->atom1->getFileID()) + " "
                    + i2s(polarHsTorsions[i][j]->atom2->getFileID()) + " "
                    + i2s(polarHsTorsions[i][j]->atom3->getFileID()) + " "
                    + i2s(polarHsTorsions[i][j]->atom4->getFileID()) + "\n";
      }
    }

    // Get sphere of acceptors for each donor
    //std::vector<atom*> atomList =  pMol->getAtomList();
    typedef std::vector<atom*>::iterator atomIterator;
    //Torsion* rotTor = 0;
    std::vector<molecule*> molList = pCol->getMoleculeList();

    std::vector<atom*> donors;
    std::vector<atom*> atoms13;
    std::vector<atom*> atoms14;

    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      if (polarHsTorsions[i].size() > 0) {
        //rotTor = polarHsTorsions[i][0];
        if (polarHsTorsions[i][0]->atom1 == polarHs[i]) {
          donors.push_back(polarHsTorsions[i][0]->atom2);
          atoms13.push_back(polarHsTorsions[i][0]->atom3);
          atoms14.push_back(polarHsTorsions[i][0]->atom4);
        }
        else if (polarHsTorsions[i][0]->atom4 == polarHs[i]) {
          donors.push_back(polarHsTorsions[i][0]->atom3);
          atoms13.push_back(polarHsTorsions[i][0]->atom2);
          atoms14.push_back(polarHsTorsions[i][0]->atom1);
        }
      }
    }

    std::vector<std::vector<atom*> > acceptors;
    std::vector<atom*>::iterator atIt;
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      std::vector<atom*> accs;
      for (unsigned int j = 0; j < molList.size(); j++) {
        std::vector<submolecule*> subMolList = molList[j]->getSubMoleculeList();
        for (unsigned int k = 0; k < subMolList.size(); k++) {
          stdFrag* pStdFg = subMolList[k]->getStdFrag();
          if (pStdFg) {
            std::vector<atom*> atomList = subMolList[k]->getAtomList();
            for (unsigned int l = 0; l < atomList.size(); l++) {
              if ((atomList[l] == polarHs[i]) or (atomList[l] == donors[i]) or (atomList[l] == atoms13[i])) continue;
              if (donors[i]->getCoords()->dist(*atomList[l]->getCoords()) < 4.0) {
                stdAtom* pStdAt = atomList[l]->getStdAtom();
                if (pStdFg->hasStdFeature(pStdAt, "HBA")) {
                  atIt = std::find(donors.begin(), donors.end(), atomList[l]);
                  if (atIt == donors.end()) {
                    //std::cout << donors[i]->getFileID() << ":" << donors[i]->getName()
                    //          << " is close to "
                    //          << atomList[l]->getFileID() << ":" << atomList[l]->getName() << std::endl;
                    accs.push_back(atomList[l]);
                  }
                  else {
                    // Add donor and polar hydrogen to list of nonbonded interactions
                    accs.push_back(atomList[l]);
                    //for (unsigned int p = 0; p < donors.size(); p++) {
                    //  if (donors[p] == atomList[l]) {
                    //    accs.push_back(polarHs[p]);
                    // }
                    //}
                  }
                }
              }
            }
          }
        }
      }
      acceptors.push_back(accs);
    }

    errMessage += "\n\nNumber of Acceptors: \n";
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      errMessage += "H:" + i2s(polarHs[i]->getFileID()) + " : " + i2s(acceptors[i].size()) + "\n";
    }

    int segments = 36;
    double *torEnergies;
    double segmentSize = static_cast<double>(360/segments);
    try {
      torEnergies = new double [nPolarHs*segments];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      //exit(0);
      throw MTKException(" Memory Allocation Failure ");
    }

    std::vector<std::vector<vector3d> > newCoordinates;

    // Calculate torsional energies
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
/*
      std::cout << polarHs[i]->getParent()->getName() << ": " << polarHs[i]->getFileID() << " "
                << donors[i]->getFileID() << " "
                << atoms13[i]->getFileID() << " "
                << atoms14[i]->getFileID() << std::endl;
*/
      std::vector<vector3d> itsCoordinates;
      //stdAtom* pStdH = polarHs[i]->getStdAtom();

      if (polarHsTorsions[i].size() > 0) {
        double dist = polarHs[i]->getCoords()->dist(*(donors[i]->getCoords()));
        double ang = angle(*(polarHs[i]->getCoords()),
                           *(donors[i]->getCoords()), *(atoms13[i]->getCoords()));

        for (int r = 0; r < segments; r++) {
          // Build atom coordinates at new dihedral angle
          buildCoord(*(polarHs[i]->getCoords()), *(donors[i]->getCoords()),
                     *(atoms13[i]->getCoords()), *(atoms14[i]->getCoords()),
                     dist, ang, static_cast<double>(segmentSize*r));

          vector3d itsCoord = *(polarHs[i]->getCoords());
          itsCoordinates.push_back(itsCoord);

          double tor_energy = 0;
          // Calculate torsional energy
          for (unsigned int j = 0; j < polarHsTorsions[i].size(); j++) {
            Torsion* pTor = polarHsTorsions[i][j];
            torsionParamList = pTor->pTorsionParamList;

            double torSize = torsion(*(pTor->atom1->getCoords()), *(pTor->atom2->getCoords()),
                                     *(pTor->atom3->getCoords()), *(pTor->atom4->getCoords()));
            if (torSize != torSize) torSize = 0.0;

            if (!torsionParamList.empty()) {
              for (torsionParamIterator c = torsionParamList.begin(); c != torsionParamList.end(); c++) {
                pTorsionParam = *c;
                if (pTorsionParam->Vn > 0.0) {
                  tor_energy += ( (pTorsionParam->Vn / pTorsionParam->npth) *
                                (1 + cos(pTorsionParam->Nt * torSize - pTorsionParam->gamma)));
                  if (tor_energy != tor_energy) {
                    //std::cout << pTorsionParam->Nt << " " <<  torSize << " " << pTorsionParam->gamma << " "
                    //          << (1 + cos(pTorsionParam->Nt * torSize - pTorsionParam->gamma) ) << std::endl;
                    //exit(0);
                    std::stringstream ss;
                    ss << pTorsionParam->Nt << " " <<  torSize << " " << pTorsionParam->gamma << " "
                              << (1 + cos(pTorsionParam->Nt * torSize - pTorsionParam->gamma) ) << std::endl;
                    std::cout << ss.str();
                    throw MTKException(ss.str());
                  }
                }
              }
            }
          }
          torEnergies[i*segments+r] = tor_energy;
          //std::cout << "     " << r*segmentSize << " " << tor_energy << std::endl;
        }
        newCoordinates.push_back(itsCoordinates);
      }
    }

/*
    errMessage += "\n\nEnergies: \n";
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      errMessage += "H: " + i2s(i) + "\n";
      for (int r = 0; r < segments; r++) {
        errMessage += "   " + d2s(torEnergies[i*segments+r]) + "\n";
      }
    }

    for (unsigned int i = 0; i < newCoordinates.size(); i++) {
      for (unsigned int j = 0; j < newCoordinates[i].size(); j++) {
        std::cout << newCoordinates[i][j] << std::endl;
      }
    }
*/

    double *eleEnergies;
    int *bestTors;

    std::vector<idObject*>  sortPolarHs;
    for (unsigned int i = 0; i < nPolarHs; i++) {
      idObject* idO = new idObject(i, 0.0);
      sortPolarHs.push_back(idO);
    }

    try {
      eleEnergies = new double [segments];
      bestTors = new int [nPolarHs];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      //exit(0);
      throw MTKException(" Memory Allocation Failure ");
    }

    double ele_energy = 0;
    for (unsigned int c = 0; c < nPolarHs; c++) {
      for (unsigned int i = 0; i < nPolarHs; i++) {
        int curH = sortPolarHs[i]->getI();
        stdAtom* pStdH = polarHs[curH]->getStdAtom();
        double stdHCharge = pStdH->atmCharge;

        if (acceptors[curH].size() > 0) {
          for (int r = 0; r < segments; r++) {
            ele_energy = 0;
            for (unsigned int j = 0; j < acceptors[curH].size(); j++) {
              stdAtom* pStdAcc = acceptors[curH][j]->getStdAtom();
              double stdAccCharge = pStdAcc->atmCharge;
              double hAccepetorDist = newCoordinates[curH][r].dist(*acceptors[curH][j]->getCoords());
              ele_energy += (stdHCharge * stdAccCharge)/hAccepetorDist;
              // Add vdW term
            }

            if ((i != 0) and (c != 0)) {
              for (unsigned int ii = 0; ii < i; ii++) {
                int curH2 = sortPolarHs[ii]->getI();
                stdAtom* pStdH2 = polarHs[curH2]->getStdAtom();
                double stdH2Charge = pStdH2->atmCharge;
                double hH2Dist = newCoordinates[curH][r].dist(newCoordinates[curH2][bestTors[curH2]]);
                ele_energy += (stdHCharge * stdH2Charge)/hH2Dist;
                // add vdW term
              }
            }
            eleEnergies[r] = ele_energy;
          }
        }
        else {
          for (int r = 0; r < segments; r++) {
            eleEnergies[r] = 0.0;
          }
        }
        double eMin = BIGNUM;
        int bestTor = 0;
        for (int x = 0; x < segments; x++) {
          if (eMin > eleEnergies[x] + torEnergies[i*segments+x]) {
            eMin = eleEnergies[x] + torEnergies[i*segments+x];
            bestTor = x;
          }
        }
        sortPolarHs[i]->setD(eMin);
        bestTors[curH] = bestTor;
      }
      // sort sortPolarHs in ascending order
      std::sort(sortPolarHs.begin(), sortPolarHs.end(), idObject::less);

      double eTotal = 0.0;
      for (unsigned int e = 0; e < sortPolarHs.size(); e++) {
        eTotal+=sortPolarHs[e]->getD();
      }
      //std::cout << " Energy " << eTotal << std::endl;
    }

    //errMessage += "Best Torsion: \n";
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      if (polarHsTorsions[i].size() > 0) {
        double dist = polarHs[i]->getCoords()->dist(*(donors[i]->getCoords()));
        double ang = angle(*(polarHs[i]->getCoords()),
                           *(donors[i]->getCoords()), *(atoms13[i]->getCoords()));
        // Build atom coordinates at new dihedral angle
        buildCoord(*(polarHs[i]->getCoords()), *(donors[i]->getCoords()),
                     *(atoms13[i]->getCoords()), *(atoms14[i]->getCoords()),
                     dist, ang, static_cast<double>(segmentSize*bestTors[i]));
        //errMessage += "H:" + i2s(i) + d2s(bestTors[i]) + "\n";
      }
    }
    errorLogger.throwError("proProtonate::optimizePolarHs", errMessage, INFO);
}

} // MTKpp namespace
