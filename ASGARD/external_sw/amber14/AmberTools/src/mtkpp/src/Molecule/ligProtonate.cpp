/*!
   \file ligProtonate.cpp
   \brief Protonates small molecule
   \author Martin Peters

   $Date: 2010/03/29 20:44:26 $
   $Revision: 1.9 $

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

#include "ligProtonate.h"

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
#include "utility.h"
#include "Log/errorHandler.h"

#include "parameters.h"
#include "Utils/idObject.h"
#include "Utils/vector3d.h"

#include <math.h>

namespace MTKpp
{

// ============================================================
// Function : ligProtonate()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
ligProtonate::ligProtonate() {}

// ============================================================
// Function : ~ligProtonate()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
ligProtonate::~ligProtonate() {}

// ============================================================
// Function : addHydrogens()
// ------------------------------------------------------------
//
// ============================================================
void ligProtonate::addHydrogens(submolecule* pSubMolecule)
{
    std::string errorMessage = "\nAdding Hydrogens without libraries\n";
    pMol = pSubMolecule->getParent();
    pCol = pMol->getParent();

    atomList = pSubMolecule->getAtomList();
    pAtom1 = 0;
    pAtom2 = 0;
    pAtom3 = 0;
    pAtom4 = 0;

    double distance = 0.0;
    double dAngle = 0.0;
    std::vector<double> torsions;
    std::vector<atom*> bondedAtoms;

    errorMessage += " Add Hydrogens to polar atoms\n";
    std::vector<atom*> acceptors;
    std::vector<atom*> donors;
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom1 = atomList[i];
      std::string symbol = pAtom1->getElement()->symbol;
      if (symbol != "O" && symbol != "N" && symbol != "S" && symbol != "F") continue;
      if (pAtom1->getValence() == pAtom1->getElement()->filledShell) {
        acceptors.push_back(pAtom1);
      }
      else if (pAtom1->getValence() < pAtom1->getElement()->filledShell) {
        if (symbol == "N") {
          if (pAtom1->getHybridization() == 4) { // sp3
            continue;
          }
          else if (pAtom1->getHybridization() == 3 and pAtom1->getType() != 2) { // non-terminal sp2
            continue;
          }
        }

        // ring systems are dealt with below, ignore them for donors
        if ((pAtom1->getType() != 5) && (pAtom1->getType() != 6)) {
          donors.push_back(pAtom1);
        }
        if ((symbol == "O") || (symbol == "S")) {
          acceptors.push_back(pAtom1);
        }
      }
    }
    for (unsigned int i = 0; i < donors.size(); i++) {
      pAtom1 = donors[i];
      bondedAtoms = pAtom1->getBondedAtoms();
      for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
        if ((bondedAtoms[j]->getAtomicNum() != 1) and (bondedAtoms[j]->getType() != 2)) {
          pAtom2 = bondedAtoms[j];
        }
      }
      if (pAtom2 != 0) {
        bondedAtoms = pAtom2->getBondedAtoms();
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          if (bondedAtoms[j] != pAtom1) {
            pAtom3 = bondedAtoms[j];
          }
        }
      }

      if (pAtom3 != 0) {
        distance = this->getBestDistance(pAtom1);
        dAngle = this->getBestAngle(pAtom1, pAtom2);
        torsions = this->getHBDonorTorsions(pAtom1, pAtom2, pAtom3, distance, dAngle, acceptors);

        pSubMolecule = pAtom1->getParent();

        for (unsigned int k = 0; k < torsions.size(); k++) {
          pAtom4 = pSubMolecule->addAtom();
          pAtom4->setElement(pMol->getParent()->pElements->getElement("H"));
          pAtom4->setFileID(pMol->getNumAtoms());
          pAtom4->setCoords(0,0,0);

          buildCoord(*(pAtom4->getCoords()), *(pAtom1->getCoords()),
                     *(pAtom2->getCoords()), *(pAtom3->getCoords()),
                     distance, dAngle, torsions[k]);

          pBond = pMol->addBond(pAtom1, pAtom4, 1, 0, 2, 0.0);
          pAtom1->addBondedAtom(pAtom4);
          pAtom4->addBondedAtom(pAtom1);
          pAtom1->setValence(pAtom1->getValence()+1);
        }
      }
      pAtom2 = 0;
      pAtom3 = 0;
    }

    errorMessage += " Add Hydrogens to Ring systems\n";
    std::vector<ring*> molRings = pMol->getRings();
    ring* r;
    unsigned int k,l,s;
    for (unsigned int i = 0; i < molRings.size(); i++) {
      r = molRings[i];
      if (!r->aromatic and r->size > 8) continue;
      s = r->atoms.size();
      for (unsigned int j = 0; j < s; j++) {
        pAtom1 = r->atoms[j];
        if (pAtom1->getValence() == pAtom1->getElement()->filledShell) continue;
        if (j+1 == s-1) {
          k = j+1;
          l = 0;
        }
        else if (j == s-1) {
          k = 0;
          l = 1;
        }
        else {
          k = j+1;
          l = j+2;
        }

        pAtom2 = r->atoms[k];
        pAtom3 = r->atoms[l];

        distance = this->getBestDistance(pAtom1);
        dAngle = this->getBestAngle(pAtom1, pAtom2, s);
        torsions = this->getBestTorsions(pAtom1, pAtom2, pAtom3);
        pSubMolecule = pAtom1->getParent();

        for (unsigned int d = 0; d < torsions.size(); d++) {
          pAtom4 = pSubMolecule->addAtom();
          pAtom4->setElement(pCol->pElements->getElement("H"));
          pAtom4->setFileID(pMol->getNumAtoms());
          pAtom4->setCoords(0,0,0);

          buildCoord(*(pAtom4->getCoords()), *(pAtom1->getCoords()),
                     *(pAtom2->getCoords()), *(pAtom3->getCoords()),
                     distance, dAngle, torsions[d]);
          errorMessage += "  H-" + i2s(pAtom1->getIndex()) + "-" +
                          i2s(pAtom2->getIndex()) + "-" +
                          i2s(pAtom3->getIndex()) + " " +
                          d2s(distance) + " " +
                          d2s(dAngle * RAD2DEG) + " " +
                          d2s(torsions[d] * RAD2DEG) + "\n";
          pBond = pMol->addBond(pAtom1, pAtom4, 1, 0, 1, 0.0);
          pAtom1->addBondedAtom(pAtom4);
          pAtom4->addBondedAtom(pAtom1);
          pAtom4->setType(1);
          pAtom1->setValence(pAtom1->getValence()+1);
        }
      }
    }

    errorMessage += " Add Hydrogens to terminal atoms\n";
    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom1 = atomList[i];
      if (pAtom1->getType() == 2 and (pAtom1->getValence() != pAtom1->getElement()->filledShell)) {
        bondedAtoms = pAtom1->getBondedAtoms();
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          if (bondedAtoms[j]->getAtomicNum() != 1) {
            pAtom2 = bondedAtoms[j];
          }
        }
        if (pAtom2 != 0) {
          bondedAtoms = pAtom2->getBondedAtoms();
          for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
            if (bondedAtoms[j]->getAtomicNum() != 1 and bondedAtoms[j] != pAtom1) {
              pAtom3 = bondedAtoms[j];
            }
          }
        }
        if (pAtom3 != 0) {
          distance = this->getBestDistance(pAtom1);
          dAngle = this->getBestAngle(pAtom1, pAtom2);

          // The following hack ensures that the -NH2 group attached to pterin's is planar with the ring
          if ((pAtom1->getElementSymbol() == "N") and
              (pAtom1->getHybridization() == 4) and
              (pAtom1->getValence() == 6) and
              (pAtom2->getType() == 6)) {
            torsions.clear();
            torsions.push_back(0.0);
            torsions.push_back(PI);
            dAngle = 120.0 * DEG2RAD;
            pAtom1->setHybridization(3);
            errorMessage += "  Setting " + i2s(pAtom1->getIndex()) + " to sp2\n";
          }
          else {
            torsions = this->getBestTorsions(pAtom1, pAtom2, pAtom3);
          }
          pSubMolecule = pAtom1->getParent();

          for (unsigned int d = 0; d < torsions.size(); d++) {
            pAtom4 = pSubMolecule->addAtom();
            pAtom4->setElement(pMol->getParent()->pElements->getElement("H"));
            pAtom4->setFileID(pMol->getNumAtoms());
            pAtom4->setCoords(0,0,0);

            buildCoord(*(pAtom4->getCoords()), *(pAtom1->getCoords()),
                       *(pAtom2->getCoords()), *(pAtom3->getCoords()),
                       distance, dAngle, torsions[d]);

            pBond = pMol->addBond(pAtom1, pAtom4, 1, 0, 2, 0.0);
            pAtom1->addBondedAtom(pAtom4);
            pAtom4->addBondedAtom(pAtom1);
            pAtom4->setType(1);
            pAtom1->setValence(pAtom1->getValence()+1);
          }
        }
      }
      pAtom2 = 0;
      pAtom3 = 0;
    }

    errorMessage += " Add Hydrogens to the remaining atoms\n";

    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom1 = atomList[i];
      if (pAtom1->getValence() != pAtom1->getElement()->filledShell) {
        bondedAtoms = pAtom1->getBondedAtoms();
        if (pAtom1->getElementSymbol() == "N" and pAtom1->getHybridization() == 4) {
          double X_N_X = 0.0;
          int nX_N_X = 0;
          for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
            for (unsigned int j2 = j+1; j2 < bondedAtoms.size(); j2++) {
              X_N_X += angle(*(bondedAtoms[j]->getCoords()),
                             *(pAtom1->getCoords()),
                             *(bondedAtoms[j2]->getCoords())) * RAD2DEG;
              nX_N_X++;
            }
          }
          if (nX_N_X > 0) {
            if (nX_N_X > 1) X_N_X = X_N_X / double(nX_N_X);
            if (X_N_X > 115.0) {
              pAtom1->setHybridization(3);
              errorMessage += " Setting " + i2s(pAtom1->getIndex()) + " to sp2 \n";
            }
          }
        }

        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          if ((bondedAtoms[j]->getAtomicNum() != 1) and
              (bondedAtoms[j]->getType() != 2)) {
            pAtom2 = bondedAtoms[j];
          }
        }
        if (pAtom2 != 0) {
          bondedAtoms = pAtom2->getBondedAtoms();
          for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
            if (bondedAtoms[j] != pAtom1) {
              pAtom3 = bondedAtoms[j];
            }
          }
        }

        if (pAtom3 != 0) {
          distance = this->getBestDistance(pAtom1);
          dAngle = this->getBestAngle(pAtom1, pAtom2);
          torsions = this->getBestTorsions(pAtom1, pAtom2, pAtom3);
          pSubMolecule = pAtom1->getParent();

          for (unsigned int d = 0; d < torsions.size(); d++) {
            pAtom4 = pSubMolecule->addAtom();
            pAtom4->setElement(pCol->pElements->getElement("H"));
            pAtom4->setFileID(pMol->getNumAtoms());
            pAtom4->setCoords(0,0,0);

            buildCoord(*(pAtom4->getCoords()), *(pAtom1->getCoords()),
                       *(pAtom2->getCoords()), *(pAtom3->getCoords()),
                       distance, dAngle, torsions[d]);

            pBond = pMol->addBond(pAtom1, pAtom4, 1, 0, 2, 0.0);
            pAtom1->addBondedAtom(pAtom4);
            pAtom4->addBondedAtom(pAtom1);
            pAtom4->setType(1);
            pAtom1->setValence(pAtom1->getValence()+1);
          }
        }
      }
      pAtom2 = 0;
      pAtom3 = 0;
    }
    errorLogger.throwError("ligProtonate::addHydrogens", errorMessage, 2);

}

// ============================================================
// Function : getBestDistance()
// ------------------------------------------------------------
//
// ============================================================
double ligProtonate::getBestDistance(atom* a)
{
    std::string sym = a->getElement()->symbol;

    // group 4
    if (sym == "C") return 1.09;

    // group 5
    if (sym == "N") return 1.008;

    // group 6
    if (sym == "O") return 0.95;
    if (sym == "S") return 1.008;
    if (sym == "Se") return 1.10;
    return 1.05;
}

// ============================================================
// Function : getBestAngle()
// ------------------------------------------------------------
//
// ============================================================
double ligProtonate::getBestAngle(atom* a, atom* b, const int &ringSize)
{
#ifdef DEBUG
    std::cout << "   protonate::getBestAngle" << std::endl;
#endif
    pBond = 0;
    pBond = pMol->getBond(a, b);

    if (pBond) {
      if (a->getType() == 6 and b->getType() == 6) { // aromatic bond,  pBond->type == 4
        if (ringSize > 0) {
#ifdef DEBUG
          std::cout << "    Adding H to a Ring of Size " << ringSize << " with an Angle of "
                    << ((360 - ((ringSize-2)*180 )/ringSize)/2) << std::endl;
#endif
          return ((360 - ((ringSize-2)*180 )/ringSize)/2) * DEG2RAD;
        }
      }
      else if (pBond->type == 1) { // single bond
        return 109.47 * DEG2RAD;
      }
      else if (pBond->type == 2) { // double bond
        return 120.00 * DEG2RAD;
      }
      else if (pBond->type == 3) { // triple bond
        return 180.00 * DEG2RAD;
      }
    }
    else {
      int aHybrid = a->getHybridization();
      if (aHybrid == 2) { // sp
        return 180.00 * DEG2RAD;
      }
      else if (aHybrid == 3) { // sp2
        return 120.00 * DEG2RAD;
      }
      else if (aHybrid == 4) { // sp3
        return 109.47 * DEG2RAD;
      }
    }
    return 109.47 * DEG2RAD;
}

// ============================================================
// Function : getBestTorsion()
// ------------------------------------------------------------
//
// ============================================================
std::vector<double> ligProtonate::getBestTorsions(atom* a, atom* b, atom* c)
{
#ifdef DEBUG
    std::cout << " protonate::getBestTorsions" << std::endl;
#endif

    std::vector<atom*> aAtomList = a->getBondedAtoms();
    atom* a1;
    std::vector<double> newTorsions;
    std::vector<double> usedTorsions;
    std::vector<double>::iterator result;
    double dTorsion;
    Bond* pBond = 0;
    pBond = pMol->getBond(a, b);
    if (!pBond) return newTorsions;

    for (unsigned int i = 0; i < aAtomList.size(); i++) {
      a1 = aAtomList[i];
      if (a1 != b) {
        dTorsion = torsion(*(a1->getCoords()), *(a->getCoords()),
                           *( b->getCoords()), *(c->getCoords()));
        // make all torsion positive (from 0 to 360 degrees)
        if (dTorsion < 0) dTorsion = dTorsion + 2*PI;
        usedTorsions.push_back(dTorsion);
      }
    }

    int numToAdd = a->getElement()->filledShell - a->getValence();

#ifdef DEBUG
    std::cout << " a: " << a->getFileID()
              << " | type = " << a->getType()
              << " | valence = " << a->getValence()
              << " | hybrid = " << a->getHybridization()
              << " | number H to add = " << numToAdd
              << std::endl;

    std::cout << " a = " << a->getFileID()
              << " - b = " << b->getFileID()
              << " - c = " << c->getFileID()
              << " | a-b Type = " << pBond->type
              << " | a-b Top  = " << pBond->topology << std::endl;

    for (unsigned int x = 0; x < usedTorsions.size(); x++) {
      std::cout << " used torsion: " << usedTorsions[x] * RAD2DEG << std::endl;
    }
#endif

    int k = 0;
    int tempNumToAdd = numToAdd;
    if (pBond->topology == 2) { // Chain atom
      if (a->getHybridization() == 4) { // sp3 atom3
        if (a->getType() == 2 and usedTorsions.size() == 0) { // Terminal atom

          std::vector<atom*> lAtomList = pMol->getAtomList();
          vector3d* pHydrogen = new vector3d(0.0);
          buildCoord(*(pHydrogen), *(a->getCoords()),
                     *(b->getCoords()), *(c->getCoords()),
                     this->getBestDistance(a), getBestAngle(a,b,0), 180.0);
          bool bOK = true;

          for (unsigned int x = 0; x < lAtomList.size(); x++) {
            double covalentBondRadii = lAtomList[x]->getElement()->covalentRadius;
            //
            // 0.23 is the covalent radius of H from elements.xml
            // BONDTOLERANCE is defined as 0.4 in constants.h
            //
            double d = 0.23 + covalentBondRadii + BONDTOLERANCE;
            vector3d* atomCoord = lAtomList[x]->getCoords();
            double dist12 = pHydrogen->dist(*atomCoord);
            if (dist12 < d) {
              bOK = false;
            }
          }
          if (bOK) {
            usedTorsions.push_back(PI);
            newTorsions.push_back(PI);
          }
          else {
            usedTorsions.push_back(0.0);
            newTorsions.push_back(0.0);
          }
          tempNumToAdd--;
          delete pHydrogen;
        }
        bool used = false;
        while (k < tempNumToAdd) {
          for (unsigned int l = 0; l < usedTorsions.size(); l++) {
            double testTorsion = usedTorsions[l] + 2*(PI/3);
            if (testTorsion > 2*PI) testTorsion = usedTorsions[l] - 4*(PI/3);
            for (unsigned int r = 0; r < usedTorsions.size(); r++) {
              if ( (std::abs(usedTorsions[r] - testTorsion) < PI/5) or
                   (std::abs(usedTorsions[r] - (testTorsion + PIt2) ) < PI/5) or
                   (std::abs(usedTorsions[r] - (testTorsion - PIt2) ) < PI/5) ) {
                used = true;
              }
            }
            if (!used) {
              newTorsions.push_back(testTorsion);
              usedTorsions.push_back(testTorsion);
              k++;
            }
            used = false;
            if (newTorsions.size() == static_cast<unsigned int>(numToAdd)) break;
          }
        }
      }
      else if (a->getHybridization() == 3) { // sp2 atoms
        if (usedTorsions.size() == 0) {
          double p[2] = {PI, 0}; // 0 and 180
          bool used = false;
          for (unsigned int j = 0; j < 2; j++) {
            for (unsigned int i = 0; i < usedTorsions.size(); i++) {
              if (std::abs(usedTorsions[i] - p[j]) < PI/5) {
                used = true;
              }
            }
            if (!used) {
              newTorsions.push_back(p[j]);
            }
            used = false;
            if (newTorsions.size() == static_cast<unsigned int>(numToAdd)) break;
          }
        }
        else if (usedTorsions.size() == 1) {
          if (usedTorsions[0] > PI) {
            newTorsions.push_back(usedTorsions[0]-PI);
          }
          else {
            newTorsions.push_back(usedTorsions[0]+PI);
          }
        }
        return newTorsions;
      }
      else if (a->getHybridization() == 2) { // sp atoms
        for (unsigned int i = 0; i < usedTorsions.size(); i++) {
          if (std::abs(usedTorsions[i] - PI) > PI/5) {
            result = std::find(newTorsions.begin(),newTorsions.end(), PI);
            if (result == newTorsions.end()) {
              newTorsions.push_back(PI); // 180.00
            }
          }
        }
        return newTorsions;
      }
    }

    k = 0;
    if (pBond->topology == 1) { // ring systems
      if ( (pBond->type == 4 ) or (pBond->type == 6) or (pBond->type == 7) ) { // aromatic rings
        for (unsigned int i = 0; i < usedTorsions.size(); i++) {
          if (std::abs(usedTorsions[i] - PI) > PI/5) {
            result = std::find(newTorsions.begin(),newTorsions.end(), PI);
            if (result == newTorsions.end()) {
              newTorsions.push_back(PI); // 180.00
            }
          }
        }
        return newTorsions;
      }
      else { // non aromatic rings
        // Known issues with 3 membered rings ...
///////
        double idealTorsion = 120.0;
        double minScore = idealTorsion;
        double score = 0.0;
        double bestTorXXXX = 0.0;
//////
        if (a->getHybridization() == 4) { // sp3 atom
          //bool used = false;

          while (k < numToAdd) {
            for (int d = 0; d < 360; d++) {
              score = 0.0;
              //std::cout << d << " " << d+360 << " " << d-360 << " ";
              for (unsigned int l = 0; l < usedTorsions.size(); l++) {
                double torInDeg = usedTorsions[l] * RAD2DEG;
                double tp = std::min( std::min(  std::abs(torInDeg - d),
                                                 std::abs(torInDeg - (d+360))
                                              ), std::abs(torInDeg - (d-360))
                                    );
                //std::cout << tp << " " << " " << std::abs(tp - idealTorsion) << std::endl;
                score += std::abs(tp - idealTorsion);
              }
              score /= double(usedTorsions.size());
              if (score < minScore) {
                minScore = score;
                bestTorXXXX = double(d);
              }
            }
            newTorsions.push_back(bestTorXXXX/RAD2DEG);
            usedTorsions.push_back(bestTorXXXX/RAD2DEG);
            k++;

            //std::cout << "   " << minScore << "   new ::: " << bestTorXXXX << std::endl;
            minScore = idealTorsion;

/*
            for (unsigned int l = 0; l < usedTorsions.size(); l++) {
              double testTorsion = usedTorsions[l] + 2*(PI/3); // 120.0
              if (testTorsion > 2*PI) testTorsion = usedTorsions[l] - 4*(PI/3);
              for (unsigned int r = 0; r < usedTorsions.size(); r++) {
                if ( (std::abs(usedTorsions[r] - testTorsion) < PI/5) or
                     (std::abs(usedTorsions[r] - (testTorsion + PIt2) ) < PI/5) or
                     (std::abs(usedTorsions[r] - (testTorsion - PIt2) ) < PI/5) ) {
                  used = true;
                }
              }

              if (!used) {
                newTorsions.push_back(testTorsion);
                usedTorsions.push_back(testTorsion);
                k++;
                break;
              }
              used = false;
            }
*/
          }
        }
        else if (a->getHybridization() == 3) { // sp2 atoms
          if (usedTorsions.size() == 0) {
            double p[2] = {0, PI}; // 0 and 180
            bool used = false;
            for (unsigned int j = 0; j < 2; j++) {
              for (unsigned int i = 0; i < usedTorsions.size(); i++) {
                if (std::abs(usedTorsions[i] - p[j]) < PI/5) {
                  used = true;
                }
              }
              if (!used) {
                newTorsions.push_back(p[j]);
              }
              used = false;
            }
          }
          else if (usedTorsions.size() == 1) {
            if (usedTorsions[0] > PI) {
              newTorsions.push_back(usedTorsions[0]-PI);
            }
            else {
              newTorsions.push_back(usedTorsions[0]+PI);
            }
          }
          return newTorsions;
        }
      }
    }
#ifdef DEBUG
    for (unsigned int b = 0; b < newTorsions.size(); b++) {
      std::cout << " torsion to be added: " << newTorsions[b] * RAD2DEG << std::endl;
    }
#endif
    return newTorsions;
}

// ============================================================
// Function : getHBDonorTorsion()
// ------------------------------------------------------------
//
// ============================================================
std::vector<double> ligProtonate::getHBDonorTorsions(atom* pAtom1, atom* pAtom2,
                    atom* pAtom3, double HDdist, double ang,
                    std::vector<atom*> acceptors)
{
    Bond* pBond = 0;
    pBond = pMol->getBond(pAtom1, pAtom2);
    if (!pBond) {
      std::cout << " Error in ligProtonate::getHBDonorTorsions ... exiting " << std::endl;
    }

#ifdef DEBUG
    std::cout << " protonate::getHBDonorTorsions" << std::endl;
    std::cout << " a: " << pAtom1->getFileID()
              << " | type = " << pAtom1->getType()
              << " | valence = " << pAtom1->getValence()
              << " | hybrid = " << pAtom1->getHybridization()
              << std::endl;

    std::cout << " a = " << pAtom1->getFileID()
              << " - b = " << pAtom2->getFileID()
              << " - c = " << pAtom3->getFileID()
              << " | a-b Type = " << pBond->type
              << " | a-b Top  = " << pBond->topology << std::endl;
#endif

    coord1 = pAtom1->getCoords();
    std::vector<double> newTorsions;
    std::vector<atom*>  aa;
    double bestHBEnergy = 0.0;
    double bestTor = 180.0;
    bool torsionChanged = false;
    double cutoffDist = 3.5;

    for (unsigned int i = 0; i < acceptors.size(); i++) {
      if (acceptors[i]->getIndex() == pAtom1->getIndex() or
          acceptors[i]->getIndex() == pAtom2->getIndex() or
          acceptors[i]->getIndex() == pAtom3->getIndex()) {
        continue;
      }
      else {
        coord2 = acceptors[i]->getCoords();
        double distance = coord1->dist(*coord2);

        if (distance <= cutoffDist) {
          aa = acceptors[i]->getBondedAtoms(); //aa- acceptor antecedants
          for (unsigned int j = 0; j < aa.size(); j++) {
            double CheckBadHBangle = angle(*(pAtom1->getCoords()),
                   *(acceptors[i]->getCoords()),  *(aa[j]->getCoords()));
            if (CheckBadHBangle < 0.785) { // 90 degrees
              continue;
            }
            vector3d* pHydrogen = new vector3d(0.0);
            for (unsigned int k = 0; k < 360; k += 6) {
              buildCoord(*(pHydrogen), *(pAtom1->getCoords()),
                         *(pAtom2->getCoords()), *(pAtom3->getCoords()),
                         HDdist, ang, static_cast<double>(k));
              double angle3 = angle(*(pAtom1->getCoords()), *(pHydrogen),
                                    *(acceptors[i]->getCoords()));
              double HAdist = pHydrogen->dist(*acceptors[i]->getCoords());
              double cos2 = cos(angle3);
              double hbe = -1 * (cos2*cos2) * exp( -1.0 * (HAdist - 2.0) *
                             (HAdist - 2.0));
              if (hbe < bestHBEnergy) {
                bestHBEnergy = hbe;
                bestTor = static_cast<double>(k);
                torsionChanged = true;
              }
            }
          }
        }
      }
    }
    if (!torsionChanged) {
#ifdef DEBUG
      std::cout << "normal torsion: " << bestTor << "  at: " <<
                   pAtom1->getElement()->symbol << " parent: " <<
                   pAtom1->getParent()->getName() << std::endl;
#endif
      newTorsions = this->getBestTorsions(pAtom1, pAtom2, pAtom3);
    }
    else {
#ifdef DEBUG
      std::cout << "adjust torsion: " << bestTor << "  at: " <<
                   pAtom1->getElement()->symbol << " parent: " <<
                   pAtom1->getParent()->getName() << std::endl;
#endif
      newTorsions.push_back(bestTor);
    }
    return newTorsions;
}

} // MTKpp namespace
