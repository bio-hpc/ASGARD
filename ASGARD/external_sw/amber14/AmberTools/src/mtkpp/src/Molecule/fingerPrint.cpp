/*!
   \file fingerPrint.cpp
   \brief Generates very simple molecular fingerprints
   \author Martin Peters

   $Date: 2010/04/29 18:59:18 $
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

#include "fingerPrint.h"

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
#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "parameters.h"

#include "Utils/vector3d.h"
#include "Utils/constants.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : fingerPrint()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
fingerPrint::fingerPrint() {}

// ============================================================
// Function : ~fingerPrint()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
fingerPrint::~fingerPrint() {}

// ============================================================
// Function : generateSimpleFP()
// ------------------------------------------------------------
// Generate simple fingerprint:
// [atom info, bond type, # of rings[ring info]]
// [H through to I (0-52),
//  B-H (53), C-H (54), N-H (55), O-H, (56) S-H (57), B-C (58),
//  B=C (59), B-O (60), B-N (61), B-F (62), B-S (63), B-Cl (64),
//  B-Br (65), B-I (66), C-C (67), C=C (68), C%C (69), N-N (70),
//  N=N (71), C-N (72), C=N (73), C%N (74), N-O (75), N=0 (76),
//  N-P (77), N-Se (78), N=Se (79), O-O (80), C-O (81), C=O (82),
//  O-Si (83), O-S (84), O=S (85), O-Se (86), O=Se (87), C-F (88),
//  C-Cl (89), C-Br (90), C-I (91), C-S (92), C=S (93), C-P (94),
//  C-Se (95), C=Se (96), S-S (97), S-N (98), S-P (99), P-P (100),
//  P-O (101), P=O (102), P-Se (103), Se-Se (104),
//  # rings (105)*[size, planarity, aromaticity, heterocyclic,
//  nNitrogen, nOxygen, nSulfur]]
// ============================================================
void fingerPrint::generateSimpleFP(molecule* pMolecule, std::vector<unsigned int> &fp)
{
    std::string errMessage = " Molecule: " + pMolecule->getName();
    errorLogger.throwError(" fingerPrint::generateSimpleFP ", errMessage, 4);
    fp.clear();

    std::vector<atom*>   atomList = pMolecule->getAtomList();
    std::map<int, Bond*> bondMap  = pMolecule->getBondMap();
    std::vector<ring*>   ringList = pMolecule->getRings();

    unsigned int nAtoms = atomList.size();
    unsigned int nRings = ringList.size();

    // Initialize vector
    for (unsigned int i = 0; i < 106; i++) {
      fp.push_back(0);
    }

    // Atoms
    int el = 0;
    for (unsigned int i = 0; i < nAtoms; i++) {
      el = atomList[i]->getElement()->number;
      if (el < 54) fp[el - 1]++;
    }

    // Bonds
    std::string at1 = "";
    std::string at2 = "";
    std::string at3 = "";
    for (std::map<int, Bond*>::iterator iter = bondMap.begin(); iter != bondMap.end(); iter++ ) {
      pBond = (*iter).second;
      at1 = pBond->atom1->getElement()->symbol;
      at2 = pBond->atom2->getElement()->symbol;

      if (at1 == "H" or at2 == "H") {
        if (at1 == "B" or at2 == "B" ) fp[53]++;
        if (at1 == "C" or at2 == "C" ) fp[54]++;
        if (at1 == "N" or at2 == "N" ) fp[55]++;
        if (at1 == "O" or at2 == "O" ) fp[56]++;
        if (at1 == "S" or at2 == "S" ) fp[57]++;
        continue;
      }

      if (at1 == "B" or at2 == "B") {
        if (at1 == "B") at3 = at2;
        if (at2 == "B") at3 = at1;
        if (at3 == "C") {
          if (pBond->type == 1) fp[58]++;
          if (pBond->type == 2) fp[59]++;
        }
        if (at3 == "O")  fp[60]++;
        if (at3 == "N")  fp[61]++;
        if (at3 == "F")  fp[62]++;
        if (at3 == "S")  fp[63]++;
        if (at3 == "Cl") fp[64]++;
        if (at3 == "Br") fp[65]++;
        if (at3 == "I")  fp[66]++;
      }

      if (at1 == "C" and at2 == "C") {
        if (pBond->type == 1) fp[67]++;
        if (pBond->type == 2) fp[68]++;
        if (pBond->type == 3) fp[69]++;
      }

      if (at1 == "N" and at2 == "N") {
        if (pBond->type == 1) fp[70]++;
        if (pBond->type == 2) fp[71]++;
      }

      if (at1 == "C" or at2 == "C") {
        if (at1 == "C") at3 = at2;
        if (at2 == "C") at3 = at1;

        if (at3 == "N") {
          if (pBond->type == 1) fp[72]++;
          if (pBond->type == 2) fp[73]++;
          if (pBond->type == 3) fp[74]++;
        }
        if (at3 == "F")  fp[88]++;
        if (at3 == "Cl") fp[89]++;
        if (at3 == "Br") fp[90]++;
        if (at3 == "I")  fp[91]++;
        if (at3 == "S") {
          if (pBond->type == 1) fp[92]++;
          if (pBond->type == 2) fp[93]++;
        }
        if (at3 == "P") fp[94]++;
        if (at3 == "Se") {
          if (pBond->type == 1) fp[95]++;
          if (pBond->type == 2) fp[96]++;
        }
      }

      if (at1 == "N" or at2 == "N") {
        if (at1 == "N") at3 = at2;
        if (at2 == "N") at3 = at1;
        if (at3 == "O") {
          if (pBond->type == 1) fp[75]++;
          if (pBond->type == 2) fp[76]++;
        }
        if (at3 == "P") fp[77]++;

        if (at3 == "Se") {
          if (pBond->type == 1) fp[78]++;
          if (pBond->type == 2) fp[79]++;
        }
      }

      if (at1 == "O" and at2 == "O") fp[80]++;

      if (at1 == "O" or at2 == "O") {
        if (at1 == "O") at3 = at2;
        if (at2 == "O") at3 = at1;
        if (at3 == "C") {
          if (pBond->type == 1) fp[81]++;
          if (pBond->type == 2) fp[82]++;
        }

        if (at3 == "Si") fp[83]++;

        if (at3 == "S") {
          if (pBond->type == 1) fp[84]++;
          if (pBond->type == 2) fp[85]++;
        }

        if (at3 == "Se") {
          if (pBond->type == 1) fp[86]++;
          if (pBond->type == 2) fp[87]++;
        }
      }

      if (at1 == "S" and at2 == "S") fp[97]++;

      if (at1 == "S" or at2 == "S") {
        if (at1 == "S") at3 = at2;
        if (at2 == "S") at3 = at1;
        if (at3 == "N") fp[98]++;
        if (at3 == "P") fp[99]++;
      }

      if (at1 == "P" and at2 == "P") fp[100]++;

      if (at1 == "P" or at2 == "P") {
        if (at1 == "P") at3 = at2;
        if (at2 == "P") at3 = at1;

        if (at3 == "O") {
          if (pBond->type == 1) fp[101]++;
          if (pBond->type == 2) fp[102]++;
        }
        if (at3 == "Se") fp[103]++;
      }

      if (at1 == "Se" and at2 == "Se") fp[104]++;
    }

    // Rings
    if (nRings > 0) {
      fp[105] = nRings;
      for (unsigned int j = 0; j < nRings; j++) {
        fp.push_back(ringList[j]->size);
        fp.push_back(ringList[j]->planar);
        fp.push_back(ringList[j]->aromatic);
        fp.push_back(ringList[j]->hetero);
        fp.push_back(ringList[j]->nNitrogen);
        fp.push_back(ringList[j]->nOxygen);
        fp.push_back(ringList[j]->nSulfur);
      }
    }
}

// ============================================================
// Function : compareSimpleFP()
// ------------------------------------------------------------
// Compare fragment fingerprints
// ============================================================
int fingerPrint::compareSimpleFP(std::vector<unsigned int> &fp1, std::vector<unsigned int> &fp2)
{
/*
#ifdef DEBUG
    std::cout << " fingerPrint::compareSimpleFP \n molecule fp " << std::endl;
    for (unsigned int j = 0; j < fp1.size(); j++) {
      std::cout << fp1[j] << " ";
    }
    std::cout << " \n fragment fp" << std::endl;
    for (unsigned int j = 0; j < fp2.size(); j++) {
      std::cout << fp2[j] << " ";
    }
    std::cout << " " << std::endl;
#endif
*/
    if (fp1.size() < 105) {
      return -1;
    }

    if (fp2.size() < 105) {
      return -2;
    }

    // Loop over all atoms, bonds
    for (unsigned int i = 0; i < 104; i++) {
      if (fp2[i] > fp1[i]) {
        //std::cout << "  atom/bond is different " << i << std::endl;
        return -3;
      }
    }

    // Rings
    if (fp2[105] > fp1[105]) {
      return -4;
    }

    unsigned int index1 = 105;
    unsigned int index2 = 105;
    int ringMatch = 0;
    for (unsigned int i = 0; i < fp2[105]; i++) {
      for (unsigned int j = 0; j < fp1[105]; j++) {
        if (fp2[index2+1] == fp1[index1+1] &&
            fp2[index2+2] == fp1[index1+2] &&
            fp2[index2+3] == fp1[index1+3] &&
            fp2[index2+4] == fp1[index1+4] &&
            fp2[index2+5] == fp1[index1+5] &&
            fp2[index2+6] == fp1[index1+6] &&
            fp2[index2+7] == fp1[index1+7] ) {
          ringMatch = 1;
        }
        index1 += 7;
      }
      if (ringMatch == 0) {
        //std::cout << " ring is different " << index2 << std::endl;
        return -5;
      }
      index2 += 7;
      index1 = 105;
      ringMatch = 0;
    }
    return 1;
}

// ============================================================
// Function : generateFragmentFP()
// ------------------------------------------------------------
// Generate fragment fingerprint
// ============================================================
void fingerPrint::generateFragmentFP(molecule* pMolecule, std::vector<unsigned int> fp)
{
    errorLogger.throwError(" fingerPrint::generateFragmentFP", " Begin ", 4);
}

// ============================================================
// Function : generate1DFP()
// ------------------------------------------------------------
// Generate 1D fingerprint
// ============================================================
void fingerPrint::generate1DFP(molecule* pMolecule, std::vector<unsigned int> fp)
{
    errorLogger.throwError(" fingerPrint::generate1DFP ", " Begin ", 4);
}

} // MTKpp namespace

