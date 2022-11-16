/*!
   \file conformer.cpp
   \brief Rotatable bond based conformer class
   \author Martin Peters

    - For each torsion we will only rotate on the rhs, X1, 4, X2, ... etc.
    \code
                        X1
                      /
                2 -- 3     X2
              /       \  /
            1          4
    \endcode
    - We need to make a list of these atoms for each torsion

    All torsions are stored in radians:
     360 deg = 2 pi radian
     180 deg = pi radian
      90 deg = pi/2 radian
       1 deg = pi/180 radian = 0.0174533 radian
      10 deg = 10 * 0.0174533 = 0.174533
    1 radian = 180/pi deg = 57.29578 deg

   $Date: 2010/03/29 20:42:28 $
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

#include "conformer.h"
#include "molecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"
#include "torsion.h"

#include "Utils/constants.h"

#include "Utils/vector3d.h"

namespace MTKpp
{

// ============================================================
// Function : conformers(molecule*)
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
conformers::conformers(molecule *parent):pParent(parent)
{
    pParent->getCoordinates(this->itsCoords);
}

// ============================================================
// Function : ~conformers()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
conformers::~conformers() {}

// ============================================================
// Function : setRotatableTorsions()
// ------------------------------------------------------------
//
// ============================================================
void conformers::addRotatableTorsion(int a, int b, int c, int d)
{
    std::vector<int> r;
    r.push_back(a);
    r.push_back(b);
    r.push_back(c);
    r.push_back(d);
    this->itsRotatableTorsions.push_back(r);
}

// ============================================================
// Function : determineAtomList()
// ------------------------------------------------------------
//
// ============================================================
void conformers::determineAtomList()
{
    for (unsigned int i = 0; i < this->itsRotatableTorsions.size(); i++) {
      std::vector<int> iTor = this->itsRotatableTorsions[i];

      std::vector<unsigned int> rhsAtoms;
      std::vector<unsigned int>::iterator unSignedIntIterator;

      std::vector<atom*> atomList;
      atomList = pParent->getAtomList();

      atom* torAtom2 = pParent->getAtom(iTor[1], 1, 0);
      atom* torAtom3 = pParent->getAtom(iTor[2], 1, 0);
      std::vector<atom*> at3_bondedAtoms = torAtom3->getBondedAtoms();

      /*
                        Y
                      /
                2 -- 3     X1
              /       \  /
            1          4
                       |
                       X2

         Here Y and 4 are added to the RHS atoms list for the torsion 1-2-3-4
      */
      for (unsigned int j = 0; j < at3_bondedAtoms.size(); j++) {
        if (at3_bondedAtoms[j] != torAtom2) {
          rhsAtoms.push_back(static_cast<unsigned int>(at3_bondedAtoms[j]->getIndex()));
        }
      }

      bool done = false;
      if (!atomList.empty()) {
        while (!done) {
          done = true;
          for (unsigned int j = 0; j < atomList.size(); j++) {
            pAtom1 = atomList[j];
            if (pAtom1 == torAtom3) continue;
            for (unsigned int k = 0; k < rhsAtoms.size(); k++) {
              pAtom2 = pParent->getAtom(rhsAtoms[k], 1, 0);
              if (pAtom1->hasBondedAtom(pAtom2)) {
                unSignedIntIterator = std::find(rhsAtoms.begin(), rhsAtoms.end(),
                               static_cast<unsigned int>(pAtom1->getIndex()));
                if (unSignedIntIterator == rhsAtoms.end()) {
                  rhsAtoms.push_back(static_cast<unsigned int>(pAtom1->getIndex()));
                  done = false;
                }
              }
            }
          }
        }
      }
      rhsAtomList.push_back(rhsAtoms);
    }
/*
#ifdef DEBUG
    for (unsigned int x = 0; x < rhsAtomList.size(); x++) {
      for (unsigned int x1 = 0; x1 < rhsAtomList[x].size(); x1++) {
        std::cout << rhsAtomList[x][x1] << " ";
      }
      std::cout << " \n" << std::endl;
    }
#endif
*/
}

// ============================================================
// Function : getConformerCoordinates()
// ------------------------------------------------------------
// Generates the atomic coordinates of the conformer
// ============================================================
void conformers::getConformerCoordinates(conformer* c, std::vector< vector3d > &newCoords)
{
    newCoords = this->itsCoords;

    if (!this->rhsAtomList.empty() and !this->itsRotatableTorsions.empty()) {
      for (unsigned int i = 0; i < this->itsRotatableTorsions.size(); i++) {

        std::vector<int> rotTor = this->itsRotatableTorsions[i];
        double rotMat[9];

        int success = formRotMat(newCoords[rotTor[0]-1], newCoords[rotTor[1]-1],
                      newCoords[rotTor[2]-1], newCoords[rotTor[3]-1],
                      c->torsions[i], rotMat);

        if (success) {
          std::vector<unsigned int> rhsAtoms = rhsAtomList[i];
          for (unsigned int k = 0; k < rhsAtoms.size(); k++) {
            newCoords[rhsAtoms[k]-1].rotateXYZ(rotMat, newCoords[rotTor[1]-1]);
          }
        }
      }
    }
}

// ============================================================
// Function : getConformerCoordinates()
// ------------------------------------------------------------
// Generates the atomic coordinates of the conformer
// ============================================================
void conformers::getConformerCoordinates(conformer* c, double coords[][3])
{
    std::vector< vector3d > newCoords;
    newCoords = this->itsCoords;

    if (!this->rhsAtomList.empty() and !this->itsRotatableTorsions.empty()) {
      for (unsigned int i = 0; i < this->itsRotatableTorsions.size(); i++) {

        std::vector<int> rotTor = this->itsRotatableTorsions[i];
        double rotMat[9];

        int success = formRotMat(newCoords[rotTor[0]-1], newCoords[rotTor[1]-1],
                      newCoords[rotTor[2]-1], newCoords[rotTor[3]-1],
                      c->torsions[i], rotMat);

        if (success) {
          std::vector<unsigned int> rhsAtoms = rhsAtomList[i];
          for (unsigned int k = 0; k < rhsAtoms.size(); k++) {
            newCoords[rhsAtoms[k]-1].rotateXYZ(rotMat, newCoords[rotTor[1]-1]);
          }
        }
      }
    }
    for (unsigned int i = 0; i < newCoords.size(); i++) {
      coords[i][0] = newCoords[i].getX();
      coords[i][1] = newCoords[i].getY();
      coords[i][2] = newCoords[i].getZ();
    }
}

// ============================================================
// Function : getConformerCenterOfMass()
// ------------------------------------------------------------
//
// ============================================================
void conformers::centerOfMass(std::vector< vector3d > coords, double c[3])
{
    vector3d  Center;

    for (unsigned int i = 0; i < coords.size(); i++) {
      Center = Center+coords[i];
    }

    Center = Center / coords.size();
    c[0] =  Center[0];
    c[1] =  Center[1];
    c[2] =  Center[2];
}

} // MTKpp namespace

