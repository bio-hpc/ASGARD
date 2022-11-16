/*!
   \file protonate.cpp
   \brief Protonates a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:44:27 $
   $Revision: 1.23 $

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

#include "protonate.h"
#include "proProtonate.h"
#include "ligProtonate.h"
#include "watProtonate.h"

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
#include "Utils/idObject.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

#include "Log/errorHandler.h"
#include "utility.h"

#include <math.h>

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : protonate()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
protonate::protonate(collection *col):pCol(col)
{
    bCol = 1;
    bMol = 0;
    pParam = 0;
    pPro = 0;
    pLig = 0;
    pWat = 0;
    this->initialize();
}

// ============================================================
// Function : protonate()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
protonate::protonate(molecule *mol):pMol(mol)
{
    pCol = mol->getParent();
    pParam = 0;
    bCol = 0;
    bMol = 1;
    pPro = 0;
    pLig = 0;
    pWat = 0;
    this->initialize();
}

// ============================================================
// Function : ~protonate()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
protonate::~protonate()
{
    if (pPro) delete pPro;
    if (pLig) delete pLig;
    if (pWat) delete pWat;
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
void protonate::initialize()
{
    errorLogger.throwError("protonate::initialize", "", INFO);

    pPro = new proProtonate();
    pLig = new ligProtonate();

    bool bInitializeWatProtonate = false;
    std::vector<molecule*> molList = pCol->getMoleculeList();
    for (unsigned int i = 0; i < molList.size(); i++) {
      molecule* pLMol = molList[i];
      std::string molName = pLMol->getName();
      if ((molName == "HOH") or (molName == "WAT") or (molName == "MOH")) {
        bInitializeWatProtonate = true;
      }
    }

    if (bInitializeWatProtonate) {
      pWat = new watProtonate();
    }
    else {
      pWat = 0;
    }
}

// ============================================================
// Function : run()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
void protonate::run()
{
    if (bCol) {
      this->runCol();
    }
    else if (bMol) {
      this->runMol(pMol);
    }
}

// ============================================================
// Function : runCol()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
void protonate::runCol()
{
    errorLogger.throwError("protonate::runCol", pCol->getName(), INFO);

    bool bPDB = false;
    std::vector<molecule*> molList = pCol->getMoleculeList();
    for (unsigned int i = 0; i < molList.size(); i++) {
      pMol = molList[i];
      if (pMol->inFileType == "pdb") bPDB = true;
      std::string molName = pMol->getName();
      if ((molName != "HOH") and (molName != "WAT") and (molName != "MOH")) {
        this->runMol(pMol);
      }
    }
    if (bPDB and pWat) {
      pWat->run(pCol);
    }
}

// ============================================================
// Function : runMol()
// ------------------------------------------------------------
//
// ============================================================
void protonate::runMol(molecule* pMolecule)
{
    errorLogger.throwError("protonate::runMol", pMolecule->getName(), INFO);

    pStdFragMinus1 = 0;
    subMoleculeList = pMolecule->getSubMoleculeList();
    bool ok = true;
    for (unsigned int j = 0; j < subMoleculeList.size(); j++) {
      pSubMolecule = subMoleculeList[j];
      pStdFrag = pSubMolecule->getStdFrag();

      if (pStdFrag) {
        if (pSubMolecule->numHeavyAtoms() != pStdFrag->numStdHeavyAtoms()) {
          int err = buildMissingHeavyAtoms(pSubMolecule);
          if (err) {
            std::string em = " Residue "
                      + pSubMolecule->getName() + ":"
                      + i2s(pSubMolecule->getSubMolId()) + " contains "
                      + i2s(pSubMolecule->numHeavyAtoms())
                      + " heavy atoms, should contain "
                      + i2s(pStdFrag->numStdHeavyAtoms());
            errorLogger.throwError("protonate::runMol", em, MTK_ERROR);
            continue;
          }
        }
        if ((pStdFrag->getType() == "s") or (pStdFrag->getType() == "m")) {
          pSubMoleculeMinus1 = 0;
          pStdFragMinus1 = 0;
        }

        if (pStdFrag->numStdAtoms() == 1) { // HACK2 (for metals and ions)
          prev3Atoms.push_back(pSubMolecule->getAtomList()[0]);
          pSubMoleculeMinus1 = pSubMolecule;
          pStdFragMinus1 = pStdFrag;
          continue;
        }

        if (pStdFragMinus1) {
          pPro->addHydrogens(pSubMolecule, pStdFrag);
        }
        else { // First residue
          std::vector<atom*> first3Atoms1;

          // Get first 3 main chain atoms
          std::vector<atom*> atomList = pMolecule->getAtomList();
          std::vector<stdAtom*> missingChainAtoms;

          for (unsigned int t = j; t < subMoleculeList.size(); t++) {
            stdFrag* curStdFrag = subMoleculeList[t]->getStdFrag();
            if (!curStdFrag) {
              std::cout << " Error in protonate::runMol " << std::endl;
              return;
            }
            std::vector<stdAtom*> stdAtomList = curStdFrag->getStdAtomList();
            for (unsigned int x = 0; x < stdAtomList.size(); x++) {
              if (stdAtomList[x]->chain == "M") {
                atom* pAtom = subMoleculeList[t]->getAtom(stdAtomList[x]);
                if (pAtom) {
                  first3Atoms1.push_back(pAtom);
                }
                else {
                  missingChainAtoms.push_back(stdAtomList[x]);
                }
                if (first3Atoms1.size() == 3) {
                  break;
                }
              }
            }
            if (first3Atoms1.size() == 3) {
              break;
            }
          }

          if (first3Atoms1.size() == 3) {
            pPro->buildMissingAtoms(pSubMolecule, missingChainAtoms, first3Atoms1);
          }
          else {
            ok = false;
          }

          std::vector<atom*> first3Atoms2;
          for (unsigned int t = j; t < subMoleculeList.size(); t++) {
            stdFrag* curStdFrag = subMoleculeList[t]->getStdFrag();
            std::vector<stdAtom*> stdAtomList = curStdFrag->getStdAtomList();
            for (unsigned int x = 0; x < stdAtomList.size(); x++) {
              if (stdAtomList[x]->chain == "M") {
                atom* pAtom = subMoleculeList[t]->getAtom(stdAtomList[x]);
                if (pAtom) {
                  first3Atoms2.push_back(pAtom);
                }
                if (first3Atoms2.size() == 3) break;
              }
            }
            if (first3Atoms2.size() == 3) break;
          }

          if (first3Atoms2.size() == 3) {
            pPro->buildDummyAtoms(first3Atoms2);
          }
          else {
            ok = false;
          }
          first3Atoms2.clear();

          if (ok) {
            pPro->addHydrogens(pSubMolecule, pStdFrag);
          }
          else {
            std::string errMess =  "\n First Residue:" +
                      pSubMolecule->getName() + " in " + pMolecule->getName()
                      + ", Failed to Build Atoms ";
            errorLogger.throwError("protonate::runMol", errMess, WARNING);
          }
        }
      }
      else { // Non standard residue, atom and bond types must be present
        pLig->addHydrogens(pSubMolecule);
      }

      pSubMoleculeMinus1 = pSubMolecule;
      pStdFragMinus1     = pStdFrag;
    }

    // Update the connectivity
    connections* pConnections = new connections(pCol);
    pConnections->assignAngles(pMolecule);
    pConnections->assignTorsions(pMolecule);
    pConnections->assignImpropers(pMolecule);

    if ((pMolecule->bMMAtomTypesAssigned) and (pMolecule->getNumAtoms() > 3)) {
      pConnections->assignStd(pMolecule);
      pPro->optimizePolarHs();
    }
    delete pConnections;
    pMolecule->bHydrogensAdded = true;
}

// ============================================================
// Function : optimizePolarHs()
// ------------------------------------------------------------
//
// ============================================================
void protonate::optimizePolarHs()
{
    std::string eMessage = "\n";

    errorLogger.throwError("protonate::optimizePolarHs", "", INFO);

    std::vector<molecule*> molList = pCol->getMoleculeList();
    typedef std::vector<atom*>::iterator atomIterator;
    typedef std::map<int, Bond*>::iterator BondMapIterator;
    std::vector<atom*> polarHs;
    std::vector<atom*> bondAtoms;
    std::vector<atom*> angleAtoms;
    std::vector<std::vector<atom*> > polarTorsions;


    for (unsigned int i = 0; i < molList.size(); i++) {
      pMol = molList[i];
      std::string molName = pMol->getName();
      if ((molName == "HOH") and (molName == "WAT") and (molName == "MOH")) {
        continue;
      }

      // Get donors (O, N attach to H by a polar bond)
      std::map<int, Bond*> bondMap =  pMol->getBondMap();
      if (!bondMap.empty()) {
        for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
          Bond* pBond = b->second;
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
                if (pBond->atom2->getElementSymbol() == "N" and pBond->atom2->getType() != 2) continue;
                polarHs.push_back(pBond->atom1);

                eMessage += " polarH: " + pSubMol1->getName() + ":" + i2s(pBond->atom1->getFileID())
                          + "@|" + pBond->atom1->getName() + "|\n";
              }
              else if (pBond->atom2->getElementSymbol() == "H") {
                if (pBond->atom1->getElementSymbol() == "N" and pBond->atom1->getType() != 2) continue;
                polarHs.push_back(pBond->atom2);

                eMessage += " polarH: " + pSubMol1->getName() + ":" + i2s(pBond->atom1->getFileID())
                          + "@|" + pBond->atom1->getName() + "|\n";
              }
            }
          }
        }
      }
    }

    unsigned int nPolarHs = polarHs.size();

    eMessage += " Total Number found = " + i2s(polarHs.size());
    errorLogger.throwError("protonate::optimizePolarHs", eMessage, INFO);
    eMessage = "";

    typedef std::map<ULONG_KIND, Torsion*>::iterator TorsionMapIterator;
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    std::vector<std::vector<Torsion*> > polarHsTorsions;
    Torsion* pTorsion = 0;
    std::vector<torsionParam*> torsionParamList;

    for (unsigned int i = 0; i < molList.size(); i++) {
      pMol = molList[i];
      std::string molName = pMol->getName();
      if ((molName == "HOH") and (molName == "WAT") and (molName == "MOH")) {
        continue;
      }

      std::map<ULONG_KIND, Torsion*> torsionMap =  pMol->getTorsionMap();
      if (!torsionMap.empty()) {
        for (unsigned int i = 0; i < nPolarHs; i++) {
          if (pMol != polarHs[i]->getParent()->getParent()) continue;
          std::vector<Torsion*> polarHsTorsion;
          for (TorsionMapIterator t = torsionMap.begin(); t != torsionMap.end(); t++) {
            pTorsion = t->second;

            torsionParamList = pTorsion->pTorsionParamList;
            if (torsionParamList.empty()) {
#ifdef DEBUG
              std::cout << " protonate: "
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
          if (polarHsTorsion.size() == 0) {
            std::cout << " Error in protonate::optimizePolarHs ... exiting " << std::endl;
            //exit(0);
            throw MTKException(" Error in protonate::optimizePolarHs ... exiting ");
          }
          polarHsTorsions.push_back(polarHsTorsion);
        }
      }
    }

#ifdef DEBUG
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      std::cout << polarHs[i]->getFileID() << ":" << polarHs[i]->getName() << std::endl;
      for (unsigned int j = 0; j < polarHsTorsions[i].size(); j++) {
        std::cout << polarHsTorsions[i][j]->atom1->getFileID() << " "
                  << polarHsTorsions[i][j]->atom2->getFileID() << " "
                  << polarHsTorsions[i][j]->atom3->getFileID() << " "
                  << polarHsTorsions[i][j]->atom4->getFileID() << std::endl;
      }
    }
#endif

    // Get sphere of acceptors for each donor
    std::vector<atom*> donors;
    std::vector<atom*> atoms13;
    std::vector<atom*> atoms14;

    for (unsigned int i = 0; i < nPolarHs; i++) {
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

    eMessage += " Donor - Acceptor List \n";
    std::vector<std::vector<atom*> > acceptors;
    std::vector<atom*>::iterator atIt;
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      std::vector<atom*> accs;
      eMessage += i2s(donors[i]->getFileID()) + ":" + donors[i]->getName() + " \n";
      for (unsigned int j = 0; j < molList.size(); j++) {
        std::vector<submolecule*> subMolList = molList[j]->getSubMoleculeList();
        for (unsigned int k = 0; k < subMolList.size(); k++) {
          stdFrag* pStdFg = subMolList[k]->getStdFrag();
          if (pStdFg) {
            stdGroup* pStdGp = pStdFg->getParent();
            std::vector<atom*> atomList = subMolList[k]->getAtomList();
            for (unsigned int l = 0; l < atomList.size(); l++) {

              if ((atomList[l] == polarHs[i]) or (atomList[l] == donors[i])) continue;

              if (polarHs[i]->has13BondedAtom(atomList[l]) or
                  polarHs[i]->has14BondedAtom(atomList[l])) continue;

              if (donors[i]->getCoords()->dist(*atomList[l]->getCoords()) < 6.0) {
                stdAtom* pStdAt = atomList[l]->getStdAtom();
                if (pStdFg->hasStdFeature(pStdAt, "HBA")) {
                  atIt = std::find(donors.begin(), donors.end(), atomList[l]);
                  if (atIt == donors.end()) {
                    // Add to list of acceptors
                    eMessage += "    " + i2s(atomList[l]->getFileID()) + ":" + atomList[l]->getName() + "\n";
                    accs.push_back(atomList[l]);
                  }
                  else {
                    // Add donor and polar hydrogen to list of nonbonded interactions
                    accs.push_back(atomList[l]);
                  }
                }
                if (pStdGp->getName() == "metals") {
                  accs.push_back(atomList[l]);
                  eMessage += "    " + i2s(atomList[l]->getFileID()) + ":" + atomList[l]->getName() + " (Metal) \n";
                }
              }

            }
          }
        }
      }

      acceptors.push_back(accs);
    }
    errorLogger.throwError("protonate::optimizePolarHs", eMessage, INFO);

    eMessage = " Number of Acceptor(s) per Donor \n";
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      eMessage +=  "   " + i2s(polarHs[i]->getFileID()) + " has " + i2s(acceptors[i].size()) + " acceptor(s) \n";
    }
    errorLogger.throwError("protonate::optimizePolarHs", eMessage, INFO);

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

    torsionParam* pTorsionParam;
    std::vector<std::vector<vector3d> > newCoordinates;

    // Calculate torsional energies
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      std::vector<vector3d> itsCoordinates;
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
        }
        newCoordinates.push_back(itsCoordinates);
      }
    }

#ifdef DEBUG
    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      for (int r = 0; r < segments; r++) {
        std::cout << torEnergies[i*segments+r] << " ";
      }
      std::cout << " " << std::endl;
    }

    for (unsigned int i = 0; i < newCoordinates.size(); i++) {
      for (unsigned int j = 0; j < newCoordinates[i].size(); j++) {
        std::cout << newCoordinates[i][j] << std::endl;
      }
    }
#endif

    double* eleEnergies;
    int* bestTors;

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

    eMessage = " Energy Minimization \n";
    int step = 1;

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
      eMessage += " Step:" + i2s(step) + " " + d2s(eTotal) + "\n";
      step++;
    }
    errorLogger.throwError("protonate::optimizePolarHs", eMessage, INFO);
    eMessage = "";

    for (unsigned int i = 0; i < polarHsTorsions.size(); i++) {
      if (polarHsTorsions[i].size() > 0) {
        double dist = polarHs[i]->getCoords()->dist(*(donors[i]->getCoords()));
        double ang = angle(*(polarHs[i]->getCoords()),
                           *(donors[i]->getCoords()), *(atoms13[i]->getCoords()));
        // Build atom coordinates at new dihedral angle
        buildCoord(*(polarHs[i]->getCoords()), *(donors[i]->getCoords()),
                     *(atoms13[i]->getCoords()), *(atoms14[i]->getCoords()),
                     dist, ang, static_cast<double>(segmentSize*bestTors[i]));
      }
    }
}

// ============================================================
// Function : buildMissingHeavyAtoms()
// ------------------------------------------------------------
// Attempt to build the missing heavy atoms of a fragment
// ============================================================
int protonate::buildMissingHeavyAtoms(submolecule* pSubMolecule)
{
    if (!pParam) {
      pParam = pCol->getParameters();
    }

    if (!pParam) {
      errorLogger.throwError("protonate::buildMissingHeavyAtoms", "Can't find required parameters", INFO);
    }

    std::string errorMessage = "\n Submolecule: " + pSubMolecule->getName() + ":" + i2s(pSubMolecule->getSubMolId());
    errorLogger.throwError("protonate::buildMissingHeavyAtoms", errorMessage, INFO);

    molecule* pLocalMol = pSubMolecule->getParent();
    collection* pLocalCol = pLocalMol->getParent();

    pStdFrag = pSubMolecule->getStdFrag();
    if (!pStdFrag) {
      std::cout << " Error in protonate::buildMissingHeavyAtoms " << std::endl;
      return 1;
    }

    int nHeavyAtoms = pSubMolecule->numHeavyAtoms();
    if (nHeavyAtoms < 3) {
      return 1;
    }
    int nStdHeavyAtoms = pStdFrag->numStdHeavyAtoms();
    errorMessage += "\n  Contains " + i2s(nHeavyAtoms) + " heavy atoms, it should have "+ i2s(nStdHeavyAtoms) + "\n";

    if (nHeavyAtoms == nStdHeavyAtoms) {
      std::cout << "  protonate::buildMissingHeavyAtoms: No atoms need to be built" << std::endl;
      return 0;
    }

    std::vector<stdAtom*> stdAtomList = pStdFrag->getStdAtomList();
    std::vector<stdAtom*> stdAtomsToAdd;
    std::vector<atom*> threeAtoms;
    std::vector<int> threeAtomsIndices;
    for (unsigned int x = 0; x < stdAtomList.size(); x++) {
      if (stdAtomList[x]->atSymbol != "H") {
        atom* pAtom = pSubMolecule->getAtom(stdAtomList[x]);
        if (!pAtom) {
          stdAtomsToAdd.push_back(stdAtomList[x]);
          errorMessage += "   Need to add " + stdAtomList[x]->identity + "\n";
        }
        else {
          if (threeAtoms.size() < 3) {
            threeAtoms.push_back(pAtom);
            threeAtomsIndices.push_back(stdAtomList[x]->index-1);
            errorMessage += "    Using " + stdAtomList[x]->identity + " as reference\n";
          }
        }
      }
    }
    unsigned int numToAdd = stdAtomsToAdd.size();

    // Create a vector of new coordinates
    std::vector<vector3d> newCoordinates;
    std::vector<int> stdAtomIndices;
    for (unsigned int i = 0; i < numToAdd; i++) {
      vector3d* itsCoord = new vector3d(0.0);
      newCoordinates.push_back(*(itsCoord));
      //std::cout << stdAtomsToAdd[i]->identity << " " << *(itsCoord) << std::endl;
      stdAtomIndices.push_back(stdAtomsToAdd[i]->index-1);
    }

    int err = pStdFrag->generateCoordinates();
    if (err) {
      std::cout << " Error in protonate::buildMissingHeavyAtoms " << std::endl;
      return err;
    }
    std::vector<vector3d*> stdFragCoords = pStdFrag->getCoordinates();

    for (unsigned int i = 0; i < numToAdd; i++) {
      stdAtom* pStdAtom = stdAtomsToAdd[i];
      double abDist = stdFragCoords[pStdAtom->index-1]->dist(*(stdFragCoords[threeAtomsIndices[0]]));
      //std::cout << abDist << std::endl;
      double abcAngle = angle(*(stdFragCoords[pStdAtom->index-1]),
                              *(stdFragCoords[threeAtomsIndices[0]]),
                              *(stdFragCoords[threeAtomsIndices[1]]));
      //std::cout << abcAngle << std::endl;

      double abcdTorsion = torsion(*(stdFragCoords[pStdAtom->index-1]),
                                   *(stdFragCoords[threeAtomsIndices[0]]),
                                   *(stdFragCoords[threeAtomsIndices[1]]),
                                   *(stdFragCoords[threeAtomsIndices[2]]));
      //std::cout << abcdTorsion << std::endl;

      buildCoord(newCoordinates[i], *(threeAtoms[0]->getCoords()),
                 *(threeAtoms[1]->getCoords()), *(threeAtoms[2]->getCoords()),
                 abDist, abcAngle, abcdTorsion);

      char temp[100];
      sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f\n",
              i+1,(pStdAtom->identity.c_str()),
              "XXX",
              "A",
              (1),
              "A",
              newCoordinates[i].getX(),newCoordinates[i].getY(),newCoordinates[i].getZ());
      errorMessage += temp;
    }

    // Carry out a dump check
    std::string errorMessage2 = "\n";
    std::vector<atom*> colAtoms = pLocalCol->getAtomList();
    for (unsigned int i = 0; i < numToAdd; i++) {
      for (unsigned int j = 0; j < colAtoms.size(); j++) {
        if (colAtoms[j]->getParent() == pSubMolecule) continue;
        if (colAtoms[j]->getCoords()->dist(newCoordinates[i]) < 1.0) {
          errorMessage2 += " New atom on top of atom " + i2s(colAtoms[j]->getFileID());
          errorLogger.throwError("protonate::buildMissingHeavyAtoms", errorMessage2, MTK_ERROR);
          return 1;
        }
      }
    }

    // Add heavy atoms to molecule
    for (unsigned int i = 0; i < numToAdd; i++) {
      atom* pAtom = pSubMolecule->addAtom();
      pAtom->setElement(pLocalCol->pElements->getElement(
             stdAtomsToAdd[i]->atSymbol));

      int ifileID = pLocalMol->getMaxFileID();
      pAtom->setFileID(ifileID+1);
      pLocalMol->setMaxFileID(ifileID+1);

      pAtom->setCoords(newCoordinates[i][0], newCoordinates[i][1],
             newCoordinates[i][2]);

      pAtom->setStdAtom(stdAtomsToAdd[i]);
      pAtom->setName(stdAtomsToAdd[i]->identity);

      atomType* pAtomType = pParam->getAtomType(stdAtomsToAdd[i]->type);
      if (!pAtomType) {
        std::string eM = "Can't find atomType for " + pAtom->getName() + " " + i2s(pAtom->getFileID());
        errorLogger.throwError("protonate::buildMissingHeavyAtoms", eM, MTK_ERROR);
        //exit(0);

        std::stringstream ss;
        ss << pAtom->getName() << std::endl;
        ss << stdAtomsToAdd[i]->type << std::endl;
        ss << " can't find atomType ...exiting..." << std::endl;
        std::cout << ss.str();
        throw MTKException(ss.str());
      }

      //std::cout << "   protonate::buildMissingHeavyAtoms::" << pSubMolecule->getName() << " "
      //          << stdAtomsToAdd[i]->identity << " Assigned " << std::endl;

      pAtom->setType(stdAtomsToAdd[i]->kind);

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

    errorLogger.throwError("protonate::buildMissingHeavyAtoms", errorMessage, INFO);

    // Update the connectivity
    connections* pConnections = new connections(pLocalCol);
    pConnections->bondByLibrary(pLocalMol, pSubMolecule, pStdFrag, 0, 0);
    //pConnections->assignBonds(pLocalMol);
    pConnections->assignAngles(pLocalMol);
    pConnections->assignTorsions(pLocalMol);
    pConnections->assignImpropers(pLocalMol);
    pConnections->assignStd(pLocalMol);
    delete pConnections;
    return 0;
}

} // MTKpp namespace
