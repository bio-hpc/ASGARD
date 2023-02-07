/*!
   \file zmatParser.cpp
   \brief Parses Z-Matrix files
   \author Martin Peters

   Reads and writes guassian files

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.7 $

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

#include <cmath>

#include "zmatParser.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"
#include "Molecule/bond.h"
#include "Molecule/angle.h"
#include "Molecule/torsion.h"
#include "Molecule/improper.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

#include "StringManip.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : zmatParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
zmatParser::zmatParser():baseParser() {}

// ============================================================
// Function : ~zmatParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
zmatParser::~zmatParser() {}

// ============================================================
// Function : Read
// ------------------------------------------------------------
// parsers a zmatrix files
// ------------------------------------------------------------
// Format:
// ============================================================
void zmatParser::Read(const std::string &zmatfile, collection* pCollection,
                      std::vector<std::vector<std::string> > &zmatrix,
                      std::map<std::string, double> &zmatData)
{
    std::ifstream izmat;
    izmat.open(zmatfile.c_str());

    if (!izmat) {
      std::cout << "\nUNABLE TO OPEN Z-MATRIX FILE"
                << "\nFILENAME = " << zmatfile << std::endl;
      return;
    }

    std::string fileline = "";
    std::vector<std::string> lineSplit;
    atom* pAtom1 = 0;
    atom* pAtom2 = 0;
    atom* pAtom3 = 0;
    atom* pAtom4 = 0;
    Bond* pBond  = 0;
    Angle* pAngle = 0;
    Torsion* pTorsion = 0;
    int nAtoms   = 1;

    molecule* pMolecule  = pCollection->addMolecule();
    pMolecule->setName("Mol");
    pMolecule->setMolId(pCollection->getNumberMolecules());

    submolecule* pSubMolecule = pMolecule->addSubMolecule();
    pSubMolecule->setName("sMol");
    pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

    // First line
    getline(izmat,fileline);
    splitString(fileline, " ", lineSplit, 0);
    if (lineSplit.size() == 1) {
      pAtom1 = pSubMolecule->addAtom();
      pAtom1->setElement(pCollection->pElements->getElement(lineSplit[0]));
      pAtom1->setFileID(nAtoms);
      zmatrix.push_back(lineSplit);
      nAtoms++;
    }
    else {
      std::cout << " Error reading z-matrix file : 1st line " << std::endl;
      return;
    }
    lineSplit.clear();

    // Second line
    getline(izmat,fileline);
    splitString(fileline, " ", lineSplit, 0);
    if (lineSplit.size() == 3) {
      pAtom1 = pSubMolecule->addAtom();
      pAtom1->setElement(pCollection->pElements->getElement(lineSplit[0]));
      pAtom1->setFileID(nAtoms);
      zmatrix.push_back(lineSplit);
      nAtoms++;

      pAtom2 = pMolecule->getAtom(atoi(lineSplit[1].c_str()), 1, 0);
      pBond = pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, 0.0);
      zmatData[lineSplit[2]] = 0.0;
    }
    else {
      std::cout << " Error reading z-matrix file : 2nd line " << std::endl;
      return;
    }
    lineSplit.clear();

    // Third line
    getline(izmat,fileline);
    splitString(fileline, " ", lineSplit, 0);
    if (lineSplit.size() == 5) {
      pAtom1 = pSubMolecule->addAtom();
      pAtom1->setElement(pCollection->pElements->getElement(lineSplit[0]));
      pAtom1->setFileID(nAtoms);
      zmatrix.push_back(lineSplit);
      nAtoms++;

      pAtom2 = pMolecule->getAtom(atoi(lineSplit[1].c_str()), 1, 0);
      pBond = pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, 0.0);
      zmatData[lineSplit[2]] = 0.0;

      pAtom3 = pMolecule->getAtom(atoi(lineSplit[3].c_str()), 1, 0);
      pAngle = pMolecule->addAngle(pAtom1, pAtom2, pAtom3, 0.0);
      zmatData[lineSplit[4]] = 0.0;
    }
    else {
      std::cout << " Error reading z-matrix file : 2nd line " << std::endl;
      return;
    }
    lineSplit.clear();

    // Fourth line to N atoms
    getline(izmat,fileline);
    splitString(fileline, " ", lineSplit, 0);
    while (lineSplit.size() == 7) {
      pAtom1 = pSubMolecule->addAtom();
      pAtom1->setElement(pCollection->pElements->getElement(lineSplit[0]));
      pAtom1->setFileID(nAtoms);
      zmatrix.push_back(lineSplit);
      nAtoms++;

      pAtom2 = pMolecule->getAtom(atoi(lineSplit[1].c_str()), 1, 0);
      pBond = pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, 0.0);
      zmatData[lineSplit[2]] = 0.0;

      pAtom3 = pMolecule->getAtom(atoi(lineSplit[3].c_str()), 1, 0);
      pAngle = pMolecule->addAngle(pAtom1, pAtom2, pAtom3, 0.0);
      zmatData[lineSplit[4]] = 0.0;

      pAtom4 = pMolecule->getAtom(atoi(lineSplit[5].c_str()), 1, 0);
      pTorsion = pMolecule->addTorsion(pAtom1, pAtom2, pAtom3, pAtom4, 0.0);
      zmatData[lineSplit[6]] = 0.0;

      lineSplit.clear();
      getline(izmat,fileline);
      splitString(fileline, " ", lineSplit, 0);
    }

    // Get bonds, angles, and torsions
    getline(izmat,fileline);
    splitString(fileline, " ", lineSplit, 0);
    while (lineSplit.size() == 2) {
      zmatData[lineSplit[0]] = strtod((char*)lineSplit[1].c_str(), 0);

      lineSplit.clear();
      getline(izmat,fileline);
      splitString(fileline, " ", lineSplit, 0);
    }

    // build coordinates
/*
\todo build coordinates.
\todo assign remaining bonds, angles, torsions, and impropers.
*/
    // Bond by distance
    //pConnections = new connections(pCollection);

    // Assign remaining angles, torsions, and impropers
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a zmatrix file.
// ------------------------------------------------------------
// Format:
// ============================================================
void zmatParser::Write(const std::string &zmatfile, molecule* pMolecule,
                       std::vector< vector3d > &coordinates)
{
    std::ofstream ogauss;
    ogauss.open(zmatfile.c_str());

    if (!ogauss or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN Z-MATRIX FILE"
                << "\nFILENAME = " << zmatfile << std::endl;
      return;
    }
}

// ============================================================
// Function : genZmatrix
// ------------------------------------------------------------
// Write a zmatrix file.
// ============================================================
int zmatParser::genZmatrix(molecule* pMolecule,
                           std::vector<std::vector<std::string> > &zmatrix,
                           std::map<std::string, double> &zmatData,
                           std::map<int, int> &atomMap)
{
#ifdef DEBUG
    std::cout << " zmatParser::genZmatrix " << std::endl;
#endif

    std::vector<std::string> line;
    std::vector<atom*> atomList = pMolecule->getAtomList();
    int nAtoms = atomList.size();
    std::vector<atom*> atomsDone;
    int nHeavyAtoms = 0;
    atom* pAtom1 = 0;
    atom* pAtom2 = 0;
    atom* pAtom3 = 0;

    // FIRST ATOM
    for (int i = 0; i < nAtoms; i++) {
      int currNHeavyAtoms = 0;
      std::vector<atom*> bondedAtoms = atomList[i]->getBondedAtoms();
      for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
        if (bondedAtoms[j]->getAtomicNum() != 1 ) {
          currNHeavyAtoms++;
        }
      }
      if (currNHeavyAtoms > nHeavyAtoms) {
        pAtom1 = atomList[i];
        nHeavyAtoms = currNHeavyAtoms;
      }
    }
    if ((pAtom1 == 0) and (nAtoms >= 1)) {
      pAtom1 = atomList[0];
    }
    line.push_back(pAtom1->getElement()->symbol);
    zmatrix.push_back(line);
    atomsDone.push_back(pAtom1);
    atomMap[pAtom1->getIndex()] = 1;
    line.clear();
#ifdef DEBUG
    std::cout << "  1st Atom " << pAtom1->getIndex() << std::endl;
#endif

    // SECOND ATOM
    std::vector<atom*> at1BondedAtoms = pAtom1->getBondedAtoms();
    for (unsigned int i = 0; i < at1BondedAtoms.size(); i++) {
      if (at1BondedAtoms[i]->getAtomicNum() != 1 ) {
        pAtom2 = at1BondedAtoms[i];
        break;
      }
    }
    if ((pAtom2 == 0) and (at1BondedAtoms.size() >= 1)) {
      pAtom2 = at1BondedAtoms[0];
    }
    line.push_back(pAtom2->getElement()->symbol);
    line.push_back("1");
    line.push_back("bd1");
    Bond* pBond = pMolecule->getBond(pAtom1, pAtom2);
    vector3d*      coord1;
    vector3d*      coord2;
    vector3d*      coord3;
    vector3d*      coord4;

    coord1 = pAtom1->getCoords();
    coord2 = pAtom2->getCoords();
    double bondLength = coord1->dist(*coord2);
    zmatData["bd1"] = bondLength;
    zmatrix.push_back(line);
    atomMap[pAtom2->getIndex()] = 2;
    atomsDone.push_back(pAtom2);
    line.clear();
#ifdef DEBUG
    std::cout << "  2nd Atom " << pAtom2->getIndex() << std::endl;
#endif

    // THIRD ATOM
    std::vector<atom*> at2BondedAtoms = pAtom2->getBondedAtoms();
    for (unsigned int i = 0; i < at2BondedAtoms.size(); i++) {
      if (at2BondedAtoms[i] != pAtom1) {
        if (at2BondedAtoms[i]->has13BondedAtom(pAtom1)) {
          pAtom3 = at2BondedAtoms[i];
          break;
        }
      }
    }
    if (pAtom3 == 0) {
      for (unsigned int i = 0; i < at1BondedAtoms.size(); i++) {
        if (at1BondedAtoms[i] != pAtom2) {
          if (at1BondedAtoms[i]->has13BondedAtom(pAtom2)) {
            pAtom3 = at1BondedAtoms[i];
            break;
          }
        }
      }
    }
    if (pAtom3 == 0) return 1;

    line.push_back(pAtom3->getElement()->symbol);
    pBond = pMolecule->getBond(pAtom2, pAtom3);
    double angleSize = 0.0;
    if (pBond) {
      line.push_back("2");
      line.push_back("bd2");
      coord1 = pAtom2->getCoords();
      coord2 = pAtom3->getCoords();
      bondLength = coord1->dist(*coord2);
      zmatData["bd2"] = bondLength;

      line.push_back("1");
      angleSize = angle(*(pAtom1->getCoords()), *(pAtom2->getCoords()),
                        *(pAtom3->getCoords()));
      line.push_back("ang1");
      zmatData["ang1"] = angleSize * RAD2DEG;
    }
    else {
      pBond = pMolecule->getBond(pAtom1, pAtom3);
      if (pBond) {
        line.push_back("1");
        line.push_back("bd2");
        coord1 = pAtom1->getCoords();
        coord2 = pAtom3->getCoords();
        bondLength = coord1->dist(*coord2);
        zmatData["bd2"] = bondLength;

        line.push_back("2");
        angleSize = angle(*(pAtom3->getCoords()), *(pAtom1->getCoords()),
                          *(pAtom2->getCoords()));
        line.push_back("ang1");
        zmatData["ang1"] = angleSize * RAD2DEG;
      }
    }
    zmatrix.push_back(line);
    atomsDone.push_back(pAtom3);
    atomMap[pAtom3->getIndex()] = 3;
    line.clear();
#ifdef DEBUG
    std::cout << "  3rd Atom " << pAtom3->getIndex() << std::endl;
#endif

    Torsion* pTorsion = 0;
    Angle* pAngle = 0;
    int nBonds = 3;
    int nAngles = 2;
    int nTorsions = 1;
    double torsionSize = 0.0;

    int allowedAttempts = 100;
    int nAttempts = 0;
    // FOURTH TO NTH ATOM
    std::vector<atom*>::iterator result;
    while (static_cast<int>(atomsDone.size()) < nAtoms) {
      nAttempts++;
      for (int j = 0; j < nAtoms; j++) {
        result = std::find(atomsDone.begin(), atomsDone.end(), atomList[j]);
        bool gotIt = false;
        if (result == atomsDone.end()) {
          for (unsigned int k = 0; k < atomsDone.size(); k++) {
            pBond = pMolecule->getBond(atomList[j], atomsDone[k]);
            if (pBond) {
              std::cout << "   Bond:" << atomList[j]->getIndex() << "-" << atomsDone[k]->getIndex() << std::endl;
              for (unsigned int l = 0; l < atomsDone.size(); l++) {
                pAngle = pMolecule->getAngle(atomList[j],
                         atomsDone[k], atomsDone[l]);
                if (pAngle) {
                 std::cout << "   Angle:" << atomList[j]->getIndex() << "-" << atomsDone[k]->getIndex() << "-"
                           << atomsDone[l]->getIndex() << std::endl;

                  for (unsigned int v = 0; v < atomsDone.size(); v++) {
                    pTorsion = pMolecule->getTorsion(atomList[j],
                               atomsDone[k], atomsDone[l], atomsDone[v]);
                    if (pTorsion) {
                      std::cout << "   Torsion:" << atomList[j]->getIndex() << "-" << atomsDone[k]->getIndex() << "-"
                                << atomsDone[l]->getIndex() << "-" << atomsDone[v]->getIndex() << std::endl;
                        gotIt = true;

                        coord1 = atomList[j]->getCoords();
                        coord2 = atomsDone[k]->getCoords();
                        coord3 = atomsDone[l]->getCoords();
                        coord4 = atomsDone[v]->getCoords();

                        bondLength = coord1->dist(*coord2);
                        angleSize = angle(*coord1, *coord2, *coord3);
                        torsionSize = torsion(*coord1, *coord2,
                                              *coord3, *coord4);
                        if (torsionSize != torsionSize) {
                          std::cout << " isnan " << std::endl;
                          //std::cout << *coord1 << " " << atomList[j]->getIndex() << " " << atomList[j]->getName() << std::endl;
                          //std::cout << *coord2 << " " << atomsDone[k]->getIndex() << " " << atomsDone[k]->getName() << std::endl;
                          //std::cout << *coord3 << " " << atomsDone[l]->getIndex() << " " << atomsDone[l]->getName() << std::endl;
                          //std::cout << *coord4 << " " << atomsDone[v]->getIndex() << " " << atomsDone[v]->getName() << std::endl;
                          //exit(0);
                          throw MTKException(" zmatParser::genZmatrix isnan error");
                          torsionSize = 0.0;
                        }

                        atomsDone.push_back(atomList[j]);
                        atomMap[atomList[j]->getIndex()] = static_cast<int>(atomsDone.size());
                        line.push_back(atomList[j]->getElement()->symbol);
                        std::stringstream sK;
                        sK << k+1;
                        line.push_back(sK.str().c_str());
                        std::stringstream sK2;
                        sK2 << nBonds;
                        std::string bndName = sK2.str().c_str();
                        bndName = "bd"+bndName;
                        line.push_back(bndName);
                        zmatData[bndName] = bondLength;

                        std::stringstream sL;
                        sL << l+1;
                        line.push_back(sL.str().c_str());
                        std::stringstream sL2;
                        sL2 << nAngles;
                        std::string angName = sL2.str().c_str();
                        angName = "ang"+angName;
                        line.push_back(angName);
                        zmatData[angName] = angleSize * RAD2DEG;

                        std::stringstream sV;
                        sV << v+1;
                        line.push_back(sV.str().c_str());
                        std::stringstream sV2;
                        sV2 << nTorsions;
                        std::string torName = sV2.str().c_str();
                        torName = "tor"+torName;
                        line.push_back(torName);
                        zmatData[torName] = torsionSize * RAD2DEG;
                        zmatrix.push_back(line);
                        line.clear();
                        nBonds++;
                        nAngles++;
                        nTorsions++;
#ifdef DEBUG
    std::cout << "  Added Atom " << atomList[j]->getIndex() << std::endl;
#endif
                        nAttempts = 0;
                      }
                    if (gotIt) break;
                  }
                }
                if (gotIt) break;
              }
            }
            if (gotIt) break;
          }
        }
        if (gotIt) break;
      }
      if (nAttempts == allowedAttempts) {
        std::cout << " Z-Matrix Build Failed .... returning 1 " << std::endl;
        return 1;
      }
    }
#ifdef DEBUG
    for (unsigned int i = 0; i < zmatrix.size(); i++) {
      for (unsigned int j = 0; j < zmatrix[i].size(); j++) {
        std::cout << zmatrix[i][j] << " ";
      }
      std::cout << " " << std::endl;
    }
    std::cout << " " << std::endl;


    typedef std::map<std::string, double>::iterator zmatDataIterator;
    for (zmatDataIterator z = zmatData.begin(); z != zmatData.end(); z++) {
      std::cout << z->first << " " << z->second << std::endl;
    }
#endif
    return 0;
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a zmatrix file.
// ============================================================
void zmatParser::Write(const std::string &zmatfile, molecule* pMolecule)
{
/*
    if (pMolecule != 0) {
      std::vector< vector3d > coordinates;
      pMolecule->getCoordinates(coordinates);
      this->Write(zmatfile, pMolecule, coordinates);
    }
*/
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a zmatrix file.
// ============================================================
void zmatParser::Write(const std::string &zmatfile,collection* pCollection, const int &molId)
{
    molecule* pMolecule = pCollection->getMolecule(molId);

    if (pMolecule != 0) {
      std::vector< vector3d > coordinates;
      pMolecule->getCoordinates(coordinates);
      this->Write(zmatfile, pMolecule, coordinates);
    }
}

} // MTKpp namespace

