/*!
   \file metalCenter.cpp
   \brief Container for metaloprotein information
   \author Martin Peters

   Container for metalloprotein information

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.10 $

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

#include "metalCenter.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "bond.h"
#include "angle.h"
#include "torsion.h"
#include "element.h"
#include "Utils/vector3d.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "parameters.h"
#include "utility.h"

#include "Utils/index.h"

namespace MTKpp
{

// ============================================================
// Function : metalCenter()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
metalCenter::metalCenter(atom* metal, std::map<std::string, double> m_d)
{
    metalAtom = metal;
    metalCoords = metalAtom->getCoords();
    metal_donors = m_d;
    assigned.assign(15, false);
    coordinationType.assign(15,"n"); // not bonded
    pMolecule = metalAtom->getParent()->getParent();
    pCollection = pMolecule->getParent();
    pStdLibrary = pCollection->getStdLibrary();
    primaryShell = "";
    secondaryShell = "";
    geometry1 = "unk";
    geometry1RMS = 0.0;
    geometry2 = "unk";
    geometry2RMS = 0.0;
    pdbFile = "";
    error = 0;
    pStdGroup = 0;
}

// ============================================================
// Function : ~metalCenter()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
metalCenter::~metalCenter() {}

// ============================================================
// Function : getMetalAtom()
// ------------------------------------------------------------
//
// ============================================================
atom* metalCenter::getMetalAtom()
{
    return this->metalAtom;
}

// ============================================================
// Function : setPdb()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::setPdb(std::string f)
{
    this->pdbFile = f;
}

// ============================================================
// Function : setError()
// ------------------------------------------------------------
//  sets the value of error
// ============================================================
void metalCenter::setError(int e)
{
    this->error = e;
}

// ============================================================
// Function : getError()
// ------------------------------------------------------------
//  Returns the value of error
// ============================================================
int metalCenter::getError()
{
    return this->error;
}

// ============================================================
// Function : addAtom()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::addAtom(atom* a)
{
    this->itsBondedAtoms.push_back(a);
    this->itsResNums.push_back(a->getParent()->getSubMolId());
    this->itsResNames.push_back(a->getParent()->getName());
}

// ============================================================
// Function : addLabel()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::addLabel(std::string s)
{
    this->itsLabels.push_back(s);
}

// ============================================================
// Function : getPrimaryShell()
// ------------------------------------------------------------
//
// ============================================================
std::string metalCenter::getPrimaryShell()
{
    return this->primaryShell;
}

// ============================================================
// Function : getPrimaryShellAtoms()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::getPrimaryShellAtoms(std::vector<atom*>& pSAtoms)
{
    pSAtoms.clear();
    for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
      if (coordinationType[i] == "p") {
        pSAtoms.push_back(this->itsBondedAtoms[i]);
      }
    }
}

// ============================================================
// Function : getSecondaryShell()
// ------------------------------------------------------------
//
// ============================================================
std::string metalCenter::getSecondaryShell()
{
    return this->secondaryShell;
}

// ============================================================
// Function : getSecondaryShellAtoms()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::getSecondaryShellAtoms(std::vector<atom*>& pSAtoms)
{
    pSAtoms.clear();
    for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
      if (coordinationType[i] == "s") {
        pSAtoms.push_back(this->itsBondedAtoms[i]);
      }
    }
}

// ============================================================
// Function : getGeometry()
// ------------------------------------------------------------
//
// ============================================================
std::string metalCenter::getGeometry()
{
    return this->geometry1;
}

// ============================================================
// Function : getGeometryRMS()
// ------------------------------------------------------------
//
// ============================================================
double metalCenter::getGeometryRMS()
{
    return this->geometry1RMS;
}

// ============================================================
// Function : checkBidentate()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::checkBidentate(int at1, int at2)
{
    atom* atom1 = this->itsBondedAtoms[at1];
    atom* atom2 = this->itsBondedAtoms[at2];

#ifdef DEBUG
    std::cout << " Check bidentate " << atom1->getName() << " "
              << atom2->getName() << std::endl;
#endif

    double myDist1 = metalCoords->dist( (*atom1->getCoords()) );
    double myDist2 = metalCoords->dist( (*atom2->getCoords()) );

    double t1 = std::abs(myDist1 - 2.0);
    double t2 = std::abs(0.04/(myDist2 - 2.0));

#ifdef DEBUG
    std::cout << "    dist1 = " << myDist1 << " dist2 = " << myDist2 << std::endl;
    std::cout << "    t1 = " << t1 << " t2 = " << t2 << std::endl;
#endif

    assigned[at1] = true;
    assigned[at2] = true;
    // Both are primary bonded
    if (std::abs(myDist1 - myDist2) < 0.2) {
      coordinationType[at1] = "p"; // primary
      coordinationType[at2] = "p";
    }
    else if (t1 - t2 < 0.15) {
      if (myDist1 < myDist2) {
        coordinationType[at1] = "p";
        coordinationType[at2] = "s"; // secondary
      }
      else {
        coordinationType[at1] = "s";
        coordinationType[at2] = "p";
      }
    }
    //else {
    //}

#ifdef DEBUG
    std::cout << " Check bidentate: Assigned " << coordinationType[at1]
              << " " << coordinationType[at2] << std::endl;
#endif

}

// ============================================================
// Function : assignBondType()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::assignBondType()
{
    for (unsigned int a = 0; a < this->itsBondedAtoms.size(); a++) {

#ifdef DEBUG
        std::cout << " metalCenter::assignBondType "
                  << this->itsBondedAtoms[a]->getName() << std::endl;
#endif

      if (!assigned[a]) {
        for (unsigned int a1 = a+1; a1 < this->itsBondedAtoms.size(); a1++) {
          if (this->itsResNames[a] == this->itsResNames[a1] and
              this->itsResNums[a] == this->itsResNums[a1] and
              itsLabels[a] == itsLabels[a1]) {
            this->checkBidentate(a, a1);

            if (this->itsResNames[a] == "HIS") {
              std::cout << " metalCenter::assignBondType "
                        << this->itsResNames[a] << "\n";
              this->error = 1;
            }
          }
        }

        if (!assigned[a]) {
          double myDist = metalCoords->dist( (*itsBondedAtoms[a]->getCoords()) );
          if (myDist < metal_donors[itsLabels[a]] + 0.5) {
            coordinationType[a] = "p";
            assigned[a] = true;
          }
          else if (myDist < metal_donors[itsLabels[a]] + 1.0) {
            coordinationType[a] = "s";
            assigned[a] = true;
          }
        }
      }

#ifdef DEBUG
        std::cout << " metalCenter::assignBondType assigned "
                  << coordinationType[a] << std::endl;
#endif

    }
    for (unsigned int a = 0; a < this->itsBondedAtoms.size(); a++) {
      if (coordinationType[a] == "p") {
        primaryShell += pStdLibrary->getL(this->itsResNames[a]);
      }
      else if (coordinationType[a] == "s") {
        secondaryShell += pStdLibrary->getL(this->itsResNames[a]);
      }
    }
}

// ============================================================
// Function : writeEnvironment()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::writeEnvironment(std::string file_name)
{
    std::vector<submolecule*> shellSubMols;
    std::vector<submolecule*>::iterator subMolIterator;

    for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
      if (coordinationType[i] == "p" or coordinationType[i] == "s") {
        submolecule* pSM = this->itsBondedAtoms[i]->getParent();
        subMolIterator = std::find(shellSubMols.begin(), shellSubMols.end(), pSM);
        if (subMolIterator == shellSubMols.end()) {
          shellSubMols.push_back(pSM);
        }
      }
    }

    // Add metal ion
    shellSubMols.push_back(metalAtom->getParent());
/*
    pdbParser* pPdbP = new pdbParser();
    molecule* pM = new molecule(pCollection);
    for (unsigned int j = 0; j < shellSubMols.size(); j++) {
      submolecule* pSM = pM->addSubMolecule();
      pSM->setSubMolId(j+1);
      pSM->setName(shellSubMols[j]->getName());
      std::vector<atom*> smAtoms = shellSubMols[j]->getAtomList();
      for (unsigned int k = 0; k < smAtoms.size(); k++) {
        atom* pA = pSM->addAtom(smAtoms[k]);
        if (!pA) {
          std::cout << " Error in writeEnvironment " << std::endl;
        }
      }
    }
    pPdbP->Write(file_name, pM);
    delete pPdbP;
*/
}

// ============================================================
// Function : assignGeometry()
// ------------------------------------------------------------
//  Assign coordination state
// ============================================================
void metalCenter::assignGeometry()
{
    /*!
        Coordination | Primary  | Secondary
        oct          |  6       |   0
    */
    if (primaryShell.size() == 6) {
      geometry1 = "oct";
      std::vector<atom*> primaryShellAtoms;
      for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
        if (coordinationType[i] == "p") {
          primaryShellAtoms.push_back(this->itsBondedAtoms[i]);
        }
      }
      this->test6Coord(primaryShellAtoms);
    }
    /*!
          Coordination | Primary  | Secondary
          tbp/ttp      |  5       |   0
          oct          |  5       |   1
    */
    else if (primaryShell.size() == 5) {
      geometry1 = "tbp";
      geometry2 = "ttp";
      std::vector<atom*> shellAtoms;
      for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
        if (coordinationType[i] == "p") {
          shellAtoms.push_back(this->itsBondedAtoms[i]);
        }
      }
      if (secondaryShell.size() == 1) {
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test6Coord(shellAtoms);
      }
      else {
        this->test5Coord(shellAtoms);
      }
    }
    /*!
        Coordination | Primary  | Secondary
        tet          |  4       |   0
        sqp          |  4       |   0
        oct          |  4       |   2
    */
    else if (primaryShell.size() == 4) {
      std::vector<atom*> shellAtoms;
      for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
        if (coordinationType[i] == "p") {
          shellAtoms.push_back(this->itsBondedAtoms[i]);
        }
      }
      if (secondaryShell.size() == 2) {
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test6Coord(shellAtoms);
      }
      else {
        this->test4Coord(shellAtoms);
      }
    }
    /*!
        Coordination | Primary  | Secondary
        trp          |  3       |   0
        tet          |  3       |   1
        oct          |  3       |   3
    */
    else if (primaryShell.size() == 3) {
      std::vector<atom*> shellAtoms;
      for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
        if (coordinationType[i] == "p") {
          shellAtoms.push_back(this->itsBondedAtoms[i]);
        }
      }

      if (secondaryShell.size() == 1) {
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test4Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 3) {
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test6Coord(shellAtoms);
      }
      else {
        geometry1 = "trp";
      }
    }
    /*!
        Coordination | Primary  | Secondary
        tet          |  2       |   2
        tbp/ttp     |  2       |   3
        oct          |  2       |   4
    */
    else if (primaryShell.size() == 2) {
      if (secondaryShell.size() == 2) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "p" or coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test4Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 3) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "p" or coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test5Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 4) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "p" or coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test6Coord(shellAtoms);
      }
      else {
        geometry1 = "unk";
      }
    }
    /*!
        Coordination | Primary  | Secondary
        tet          |  1       |   3
        ttp          |  1       |   4
        oct          |  1       |   5
    */
    else if (primaryShell.size() == 1) {
      if (secondaryShell.size() == 3) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "p" or coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test4Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 4) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "p" or coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test5Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 5) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "p" or coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test6Coord(shellAtoms);
      }
      else {
        geometry1 = "unk";
      }
    }
    /*!
        Coordination | Primary  | Secondary
        tet          |  0       |   4
        ttp/tbp      |  0       |   5
        oct          |  0       |   6
    */
    else if (primaryShell.size() == 0) {
      if (secondaryShell.size() == 4) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test4Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 5) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test5Coord(shellAtoms);
      }
      else if (secondaryShell.size() == 6) {
        std::vector<atom*> shellAtoms;
        for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
          if (coordinationType[i] == "s") {
            shellAtoms.push_back(this->itsBondedAtoms[i]);
          }
        }
        this->test6Coord(shellAtoms);
      }
      else {
        geometry1 = "unk";
      }
    }
    else {
      geometry1 = "unk";
    }
}

// ============================================================
// Function : print()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::print(std::ostream& os)
{
    std::string priShell = this->getPrimaryShell();
    std::sort(priShell.begin(), priShell.end());
    std::string secShell = this->getSecondaryShell();
    std::sort(secShell.begin(), secShell.end());

    for (unsigned int a = 0; a < this->itsBondedAtoms.size(); a++) {
      double myDistance = metalCoords->dist( (*itsBondedAtoms[a]->getCoords()) );
      if (coordinationType[a] != "n") {
        os << pdbFile << ","
           << metalAtom->getParent()->getName() << ","
           << metalAtom->getParent()->getSubMolId() << ","
           << metalAtom->getName() << ","
           << metalAtom->getFileID() << ","
           << itsBondedAtoms[a]->getParent()->getName() << ","
           << itsBondedAtoms[a]->getParent()->getSubMolId() << ","
           << itsBondedAtoms[a]->getName() << ","
           << itsBondedAtoms[a]->getFileID() << ","
           << myDistance << ","
           << coordinationType[a] << ","
           << itsBondedAtoms[a]->getTempFactor() << ","
           << priShell << ","
           << secShell << ","
           << geometry1 << ","
           << geometry1RMS << ","
           << geometry2 << ","
           << geometry2RMS
           << std::endl;
      }
    }
}

// ============================================================
// Function : getInfo()
// ------------------------------------------------------------
//
// ============================================================
std::string metalCenter::getInfo()
{
    std::string priShell = this->getPrimaryShell();
    std::sort(priShell.begin(), priShell.end());
    std::string secShell = this->getSecondaryShell();
    std::sort(secShell.begin(), secShell.end());

    std::string info = "";

    info += pdbFile + "," +
            priShell + "," +
            secShell + "," +
            geometry1 + "," +
            d2s(geometry1RMS) + "," +
            geometry2 + "," +
            d2s(geometry2RMS) + "\n";

    for (unsigned int a = 0; a < this->itsBondedAtoms.size(); a++) {
      double myDistance = metalCoords->dist( (*itsBondedAtoms[a]->getCoords()) );
      if (coordinationType[a] != "n") {
        info += pdbFile + "," +
            metalAtom->getParent()->getName() + "," +
                i2s(metalAtom->getParent()->getSubMolId()) + "," +
                metalAtom->getName() + "," +
                i2s(metalAtom->getFileID()) + "," +
                itsBondedAtoms[a]->getParent()->getName() + "," +
                i2s(itsBondedAtoms[a]->getParent()->getSubMolId()) + "," +
                itsBondedAtoms[a]->getName() + "," +
                i2s(itsBondedAtoms[a]->getFileID()) + "," +
                d2s(myDistance) + "," +
                coordinationType[a] + "," +
                d2s(itsBondedAtoms[a]->getTempFactor()) + "\n";
      }
    }
    return info;
}
// ============================================================
// Function : test4Coord()
// ------------------------------------------------------------
//  Test if 4/5 coordinate is tet, sqp, or tnb
// ============================================================
void metalCenter::test4Coord(std::vector<atom*> atoms)
{
    double deltaTet = 0.0;
    double deltaSqp = 0.0;
    double LargestAngle = 0.0;
    int largestAngleI = 0;
    int largestAngleJ = 0;
    double angles[16];
    int sqpAngles[16];
    int fourCoordIndex = 0;

    //std::cout << "\n\n metalCenter::test4Coord \n";
    for (unsigned int i = 0; i < 4; i++) {
      for (unsigned int j = i+1; j < 4; j++) {
        fourCoordIndex = i * 4 + j;

        double myAngle = RAD2DEG * angle(*(atoms[i]->getCoords()),
                                 *(metalCoords),
                                 *(atoms[j]->getCoords()));
        angles[fourCoordIndex] = myAngle;
        sqpAngles[fourCoordIndex] = 0;

        deltaTet += pow((myAngle - 109.5),2);
        //std::cout << i << " " << j << " " << fourCoordIndex << " " << myAngle << "\n";
        if (myAngle > LargestAngle) {
          largestAngleI = i;
          largestAngleJ = j;
          LargestAngle = myAngle;
        }
      }
    }
    deltaTet = sqrt(deltaTet / 6.0);
    //std::cout << " largestAngleI " << largestAngleI << " largestAngleJ " << largestAngleJ << "\n";
    sqpAngles[largestAngleI * 4 + largestAngleJ] = 1;
    std::vector<atom*> sqpAtoms;
    sqpAtoms.push_back(atoms[largestAngleI]);
    sqpAtoms.push_back(atoms[largestAngleJ]);

    int sqp[2];
    int sqpIndex = 0;
    for (unsigned int i = 0; i < 4; i++) {
      if (atoms[i] != sqpAtoms[0] and atoms[i] != sqpAtoms[1]) {
        sqp[sqpIndex] = i;
        sqpIndex++;
      }
    }
    sqpAngles[sqp[0] * 4 + sqp[1]] = 1;

    for (unsigned int i = 0; i < 4; i++) {
      for (unsigned int j = i+1; j < 4; j++) {
        fourCoordIndex = i * 4 + j;
        if (sqpAngles[fourCoordIndex]) {
          deltaSqp += pow((angles[fourCoordIndex] - 180.0),2);
          //std::cout << angles[fourCoordIndex] << "(180.0) " << deltaSqp << std::endl;
        }
        else {
          deltaSqp += pow((angles[fourCoordIndex] - 90.0),2);
          //std::cout << angles[fourCoordIndex] << "(90.0) " << deltaSqp << std::endl;
        }
      }
    }

    deltaSqp = sqrt(deltaSqp / 6.0);
    if (deltaTet < deltaSqp) {
      geometry1 = "tet";
      geometry2 = "sqp";
      geometry1RMS = deltaTet;
      geometry2RMS = deltaSqp;

      if (primaryShell.size() == 4 and secondaryShell.size() == 1) {
        geometry1 = "tnb";
        geometry2 = "";
      }
    }
    else {
      geometry1 = "sqp";
      geometry2 = "tet";
      geometry1RMS = deltaSqp;
      geometry2RMS = deltaTet;
    }
    //std::cout << geometry1 << " " << geometry1RMS << " " << geometry2 << " " << geometry2RMS << std::endl;
}

// ============================================================
// Function : test5Coord()
// ------------------------------------------------------------
//  Determine the tbp_rms and ttp_rms
// ============================================================
void metalCenter::test5Coord(std::vector<atom*> atoms)
{
    double deltaTbp = 0.0;
    double deltaTtp = 0.0;
    double largestAngle = 0.0;
    int largestAngleIndex = 0;
    double largestAngle2 = 0.0;
    int largestAngleIndex2 = 0;

    double angles[25];
    for (unsigned int i = 0; i < 25; i++) {
      angles[i] = 0.0;
    }

    int largeAngles[4];
    for (unsigned int i = 0; i < 4; i++) {
      largeAngles[i] = 0;
    }
    int fiveCoordIndex = 0;

    for (unsigned int i = 0; i < 5; i++) {
      for (unsigned int j = i+1; j < 5; j++) {
        fiveCoordIndex = i * 5 + j;
        double myAngle = RAD2DEG * angle(*(atoms[i]->getCoords()),
                                 * (metalCoords),
                                 * (atoms[j]->getCoords()));
        angles[fiveCoordIndex] = myAngle;
        angles[j*5+i] = myAngle;

        double smallestAngle = 180.0;
        int indexSmallestAngle = 0;
        for (unsigned int n = 0; n < 4; n++) {
          if (angles[largeAngles[n]] < smallestAngle) {
            smallestAngle = angles[largeAngles[n]];
            indexSmallestAngle = n;
          }
        }

        if (myAngle > smallestAngle) {
          largeAngles[indexSmallestAngle] = fiveCoordIndex;
        }
        if (myAngle > largestAngle) {
          largestAngle = myAngle;
          largestAngleIndex = fiveCoordIndex;
        }
        if (myAngle > largestAngle2 and myAngle < largestAngle) {
          largestAngle2 = myAngle;
          largestAngleIndex2 = fiveCoordIndex;
        }
      }
    }

    for (unsigned int i = 0; i < 5; i++) {
      for (unsigned int j = i+1; j < 5; j++) {
        fiveCoordIndex = i * 5 + j;
        bool bLarge = false;
        for (unsigned int n = 0; n < 4; n++) {
          if (largeAngles[n] == fiveCoordIndex) {
            if (largeAngles[n] == largestAngleIndex) {
              deltaTbp += pow((angles[fiveCoordIndex] - 180.0),2);
            }
            else {
              deltaTbp += pow((angles[fiveCoordIndex] - 120.0),2);
            }
            bLarge = true;
          }
        }
        if (!bLarge) {
          deltaTbp += pow((angles[fiveCoordIndex] - 90.0),2);
        }
      }
    }
    deltaTbp = sqrt(deltaTbp / 10.0);

    unsigned int apex = 0;
    std::vector<unsigned int> used;
    for (unsigned int i = 0; i < 5; i++) {
      for (unsigned int j = i+1; j < 5; j++) {
        fiveCoordIndex = i * 5 + j;
        if (fiveCoordIndex == largestAngleIndex or fiveCoordIndex == largestAngleIndex2) {
          intVectorIterator vIt_i = std::find(used.begin(), used.end(), i);
          intVectorIterator vIt_j = std::find(used.begin(), used.end(), j);
          if (vIt_i == used.end() and vIt_j == used.end()) {
            used.push_back(i);
            used.push_back(j);
          }
        }
      }
    }

    for (unsigned int i = 0; i < 5; i++) {
      intVectorIterator vIt = std::find(used.begin(), used.end(), i);
      if (vIt == used.end()) {
        apex = i;
      }
    }

    double b_m = 0.0;
    for (unsigned int i = 0; i < 5; i++) {
      if (i == apex) continue;
      fiveCoordIndex = i * 5 + apex;
      b_m += angles[fiveCoordIndex];
    }
    double idealValue = b_m / 4;

    for (unsigned int i = 0; i < 5; i++) {
      if (i == apex) continue;
      fiveCoordIndex = i * 5 + apex;
      deltaTtp += pow((angles[fiveCoordIndex] - idealValue),2);
    }

    idealValue = 360.0 - (2*idealValue);
    for (unsigned int i = 0; i < 5; i++) {
      for (unsigned int j = i+1; j < 5; j++) {
        fiveCoordIndex = i * 5 + j;
        if (fiveCoordIndex == largestAngleIndex or fiveCoordIndex == largestAngleIndex2) {
          deltaTtp += pow((angles[fiveCoordIndex] - idealValue),2);
        }
      }
    }
    idealValue = b_m / 4;
    idealValue = (2.0*asin(pow(2.0,-0.5) * sin((180.0 - idealValue) * DEG2RAD))) * RAD2DEG;

    for (unsigned int i = 0; i < 5; i++) {
      if (i == apex) continue;
      for (unsigned int j = i+1; j < 5; j++) {
        if (j == apex) continue;
        fiveCoordIndex = i * 5 + j;
        if (fiveCoordIndex == largestAngleIndex or fiveCoordIndex == largestAngleIndex2) continue;
        deltaTtp += pow((angles[fiveCoordIndex] - idealValue),2);
      }
    }

    if (deltaTbp < deltaTtp) {
      geometry1 = "tbp";
      geometry2 = "ttp";
      geometry1RMS = deltaTbp;
      geometry2RMS = deltaTtp;
    }
    else {
      geometry1 = "ttp";
      geometry2 = "tbp";
      geometry1RMS = deltaTtp;
      geometry2RMS = deltaTbp;
    }
}

// ============================================================
// Function : test6Coord()
// ------------------------------------------------------------
//  Determine the oct_rms
// ============================================================
void metalCenter::test6Coord(std::vector<atom*> atoms)
{
    double deltaOct = 0.0;
    double angles[36];
    for (unsigned int i = 0; i < 36; i++) {
      angles[i] = 0.0;
    }
    int angles180[3];
    for (unsigned int i = 0; i < 3; i++) {
      angles180[i] = 0;
    }
    int sixCoordIndex = 0;

    for (unsigned int i = 0; i < 6; i++) {
      for (unsigned int j = i+1; j < 6; j++) {
        sixCoordIndex = i * 6 + j;
        double myAngle = RAD2DEG * angle(*(atoms[i]->getCoords()),
                                 * (metalCoords),
                                 * (atoms[j]->getCoords()));
        angles[sixCoordIndex] = myAngle;

        double smallestAngle = 180.0;
        int indexSmallestAngle = 0;
        for (unsigned int n = 0; n < 3; n++) {
          if (angles[angles180[n]] < smallestAngle) {
            smallestAngle = angles[angles180[n]];
            indexSmallestAngle = n;
          }
        }

        if (myAngle > smallestAngle) {
          angles180[indexSmallestAngle] = sixCoordIndex;
        }
      }
    }

    for (unsigned int i = 0; i < 6; i++) {
      for (unsigned int j = i+1; j < 6; j++) {
        sixCoordIndex = i * 6 + j;
        bool b180 = false;
        for (unsigned int n = 0; n < 3; n++) {
          if (angles180[n] == sixCoordIndex) {
            deltaOct += pow((angles[sixCoordIndex] - 180.0),2);
            b180 = true;
          }
        }
        if (!b180) {
          deltaOct += pow((angles[sixCoordIndex] - 90.0),2);
        }
      }
    }
    deltaOct = sqrt(deltaOct / 15.0);
    geometry1RMS = deltaOct;
}

// ============================================================
// Function : setStdGroup()
// ------------------------------------------------------------
//
// ============================================================
void metalCenter::setStdGroup(stdGroup* pStdG)
{
    this->pStdGroup = pStdG;
}

// ============================================================
// Function : getStdGroup()
// ------------------------------------------------------------
//
// ============================================================
stdGroup* metalCenter::getStdGroup()
{
    return this->pStdGroup;
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// BONDS //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// ============================================================
// Function : addBond()
// ------------------------------------------------------------
//
// ============================================================
Bond* metalCenter::addBond(atom* at1, atom* at2, const int &type,
                    const int &stereo, const int &topology, const double &b)
{
    if (at1 and at2) {
      if (!this->hasBond(at1, at2)) {
        pBond = new Bond();
        pBond->atom1 = at1;
        pBond->atom2 = at2;
        pBond->type  = type;
        pBond->stereo = stereo;
        pBond->topology = topology;
        if (b == 0.0) {
          vector3d* coords1 = at1->getCoords();
          double d = coords1->dist((*at2->getCoords()));
          pBond->size  = d;
        }
        pBond->rotatable = 1;

        int bondIndex = indexAB(at1->getColIndex(), at2->getColIndex(), MAXATOMS);
        this->itsBondMap[bondIndex] = pBond;
        return pBond;
      }
    }
    return 0;
}

// ============================================================
// Function : addBond()
// ------------------------------------------------------------
//
// ============================================================
bool metalCenter::hasBond(atom* at1, atom* at2)
{
    int bondIndex = indexAB(at1->getColIndex(), at2->getColIndex(), MAXATOMS);
    BondMapIterator b = this->itsBondMap.find(bondIndex);

    if (b != this->itsBondMap.end()) {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getBondMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<int, Bond*> metalCenter::getBondMap()
{
    return this->itsBondMap;
}

// ============================================================
// Function : numBonds()
// ------------------------------------------------------------
//
// ============================================================
int metalCenter::numBonds()
{
    if (!this->itsBondMap.empty()) {
      return this->itsBondMap.size();
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ANGLES /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// ============================================================
// Function : addAngle()
// ------------------------------------------------------------
//
// ============================================================
Angle* metalCenter::addAngle(atom* at1, atom* at2, atom* at3,
                             const double& /*ang*/)
{
    pAngle = new Angle();
    pAngle->atom1 = at1;
    pAngle->atom2 = at2;
    pAngle->atom3 = at3;

    double dAngle = angle(*(at1->getCoords()), *(at2->getCoords()), *(at3->getCoords()));
    pAngle->size  = dAngle;
/*
    if (ang != dAngle) {
      std::cout << " metalCenter::addAngle warning " << std::endl;
    }
*/
    try {
      ULONG_KIND angleIndex = indexABC_ULL(at1->getColIndex(), at2->getColIndex(),
                                at3->getColIndex(), MAXATOMS, MAXATOMS);

      this->itsAngleMap[angleIndex] = pAngle;
    }
    catch (MTKException& e) {
      std::cout << "Error in addAngle: " << e.message << std::endl;
    }

    return pAngle;
}

// ============================================================
// Function : hasAngle()
// ------------------------------------------------------------
//
// ============================================================
bool metalCenter::hasAngle(atom* at1, atom* at2, atom* at3)
{
    try {
      ULONG_KIND angleIndex = indexABC_ULL(at1->getColIndex(), at2->getColIndex(),
                                at3->getColIndex(), MAXATOMS, MAXATOMS);

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
// Function : getAngleMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<ULONG_KIND, Angle*> metalCenter::getAngleMap()
{
    return this->itsAngleMap;
}

// ============================================================
// Function : numAngles()
// ------------------------------------------------------------
//
// ============================================================
int metalCenter::numAngles()
{
    if (!this->itsAngleMap.empty()) {
      return this->itsAngleMap.size();
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// TORSIONS ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// ============================================================
// Function : addTorsion()
// ------------------------------------------------------------
//
// ============================================================
Torsion* metalCenter::addTorsion(atom* at1, atom* at2, atom* at3, atom* at4,
                                 const double &tor)
{
    pTorsion = new Torsion();
    pTorsion->atom1 = at1;
    pTorsion->atom2 = at2;
    pTorsion->atom3 = at3;
    pTorsion->atom4 = at4;
    pTorsion->size = tor;

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
// Function : hasTorsion()
// ------------------------------------------------------------
//
// ============================================================
bool metalCenter::hasTorsion(atom* at1, atom* at2, atom* at3, atom* at4)
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
// Function : getTorsionMap()
// ------------------------------------------------------------
//
// ============================================================
std::map<ULONG_KIND, Torsion*> metalCenter::getTorsionMap()
{
    return this->itsTorsionMap;
}

// ============================================================
// Function : numTorsions()
// ------------------------------------------------------------
//
// ============================================================
int metalCenter::numTorsions()
{
    return this->itsTorsionMap.size();
}

// ============================================================
// Function : metalGroup()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
metalGroup::metalGroup()
{

}

// ============================================================
// Function : ~metalGroup()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
metalGroup::~metalGroup() {}

// ============================================================
// Function : addMetalCenter()
// ------------------------------------------------------------
//
// ============================================================
void metalGroup::addMetalCenter(metalCenter* m)
{
    this->metalCenters.push_back(m);
}

// ============================================================
// Function : getMetalCenters()
// ------------------------------------------------------------
//
// ============================================================
std::vector<metalCenter*> metalGroup::getMetalCenters()
{
    return this->metalCenters;
}

} // MTKpp namespace
