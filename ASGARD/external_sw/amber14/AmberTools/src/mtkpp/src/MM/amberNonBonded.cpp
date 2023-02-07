/*!
   \file amberNonBonded.cpp
   \brief AMBER non-bonded energy and gradient
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
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

//=============================================================
// description:
//               ___  [    A        C      q1 q2     ]
// E          =  \    [   ---   -  ---  +  -----     ]
//  nb pairs     /__  [   r^12     r^6      D r      ]
//             nonbonded
//              pairs
//
// In the AMBER potential 1-4 interactions are included in the non-bonded term,
// there are scaled by 1/1.2.  No 1-2 or 1-3 terms are included.
// All other interactions are scaled by 1.

// Lennard-Jones Potential:
//
// 1/r^6  --> attractive part (-ve)
//        --> dominates at large distances
//
// 1/r^12 --> repulsive part (+ve)
//        --> dominates at short distances
//
// Coulomb Potential:
//
// if q1 and q2 are the same sign it leads to a repulsive interaction
// however if they have different signs an attractive interaction results
//=============================================================

#include "amberNonBonded.h"
#include "amber.h"
#include "math.h"
#include <iostream>

#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : amberNonBonded()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberNonBonded::amberNonBonded() {}

// ============================================================
// Function : amberNonBonded()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberNonBonded::amberNonBonded(amber* p)
{
    pAmber = p;
}

// ============================================================
// Function : ~amberNonBonded()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
amberNonBonded::~amberNonBonded() {}

// ============================================================
// Function : calculateE()
// ------------------------------------------------------------
//
// ============================================================
double amberNonBonded::calculateE()
{
    int nAtoms = pAmber->getNumAtoms();
    int nUniqueTypes = pAmber->getNumTypes();
    int *types = pAmber->getIntTypes();
    int *nonBondedParameterIndex = pAmber->getNonBondedParameterIndex();

    double *xyz = pAmber->getCoords();
    double *charges = pAmber->getCharges();
    double *r6Params = pAmber->getR6Params();
    double *r12Params = pAmber->getR12Params();

    int *numExcluded = pAmber->getNumExcluded();
    int *excluded = pAmber->getExcluded();

    int *numExcluded14 = pAmber->getNumExcluded14();
    int *excluded14 = pAmber->getExcluded14();

    int *atomFlags = pAmber->getAtomFlags();

/*
std::cout << " amberNonBonded::calculateE \n"
          << " nAtoms = " << nAtoms
          << " nUniqueTypes = " << nUniqueTypes << std::endl;
*/

    double nonBondedEnergy = 0;
    double atom1Charge = 0.0;
    double atom1_14Charge = 0.0;
    double atom2Charge = 0.0;
    double r = 0.0;
    double r2 = 0.0;
    double r6 = 0.0;
    double r12 = 0.0;
    double vdw_energy = 0.0;
    double ele_energy = 0.0;

    double R6 = 0.0;
    double R12 = 0.0;
    double vdW14Energy = 0.0;
    double ele14Energy= 0.0;
    double vdWEnergy = 0.0;
    double eleEnergy = 0.0;

    double r6Param = 0.0;
    double r12Param = 0.0;

    int excludedAtomCounter = 0;
    int excluded14AtomCounter = 0;
    for (int i = 0; i < nAtoms; i++) {
      int atomIType = types[i];
      setToOne(atomFlags, nAtoms);

      int nExcludedAtoms = numExcluded[i];
      for (int e = 0; e < nExcludedAtoms; e++) {
        atomFlags[excluded[excludedAtomCounter]] = 0;
        excludedAtomCounter++;
      }

      int nExcludedAtoms14 = numExcluded14[i];
      for (int e = 0; e < nExcludedAtoms14; e++) {
        atomFlags[excluded14[excluded14AtomCounter]] = -1;
        excluded14AtomCounter++;
      }

      atom1Charge = charges[i];
      atom1_14Charge = charges[i] * 0.83333; // scale factor, ELE == 1.2

      for (int j = i+1; j < nAtoms; j++) {
        if (atomFlags[j] == 1) {
          int atomJType = types[j];
          int paramIndex = nonBondedParameterIndex[atomIType * nUniqueTypes + atomJType];

          r6Param = r6Params[paramIndex];   // 2 * epsilon * sigma
          r12Param = r12Params[paramIndex]; // epsilon *sigma^2
          r = sqrt ( pow((xyz[i*3  ] - xyz[j*3  ]),2) +
                     pow((xyz[i*3+1] - xyz[j*3+1]),2) +
                     pow((xyz[i*3+2] - xyz[j*3+2]),2) );
          r2 = r*r;

          // van der Waals
          r6 = r6Param / pow(r2, 3);
          r12 = r12Param / pow(r2, 6);
          vdw_energy = r12 - r6;

          // Electrostatic
          atom2Charge = charges[j];
          ele_energy = atom2Charge / r;

          // van der Waals
          R6  -= r6;
          R12 += r12;
          vdWEnergy += vdw_energy;

          // Electrostatic
          eleEnergy += atom1Charge * ele_energy;
        }
        else if (atomFlags[j] == -1) { // 1-4 Interaction
          int atomJType = types[j];
          int paramIndex = nonBondedParameterIndex[atomIType * nUniqueTypes + atomJType];

          r6Param = r6Params[paramIndex];   // 2 * epsilon * sigma
          r12Param = r12Params[paramIndex]; // epsilon *sigma^2
          r = sqrt ( pow((xyz[i*3  ] - xyz[j*3  ]),2) +
                     pow((xyz[i*3+1] - xyz[j*3+1]),2) +
                     pow((xyz[i*3+2] - xyz[j*3+2]),2) );
          r2 = r*r;

          // van der Waals
          r6 = r6Param / pow(r2, 3);
          r12 = r12Param / pow(r2, 6);
          vdW14Energy += 0.5 * (r12 - r6);

          R6  -= 0.5 * r6;
          R12 += 0.5 * r12;

          // Electrostatic
          atom2Charge = charges[j];
          ele14Energy += atom1_14Charge * (atom2Charge / r);
        }
      }
    }
    nonBondedEnergy = vdWEnergy + eleEnergy + vdW14Energy + ele14Energy;

    pAmber->setR6(R6);
    pAmber->setR12(R12);
    pAmber->setVDW(vdWEnergy);
    pAmber->setEle(eleEnergy);

    pAmber->setVDW14(vdW14Energy);
    pAmber->setEle14(ele14Energy);

    return nonBondedEnergy;
}

// ============================================================
// Function : calculateE()
// ------------------------------------------------------------
//
// ============================================================
/*
double amberNonBonded::calculateE()
{
    int nAtoms = pAmber->getNumAtoms();
    double *xyz = pAmber->getCoords();
    double *charges = pAmber->getCharges();

    int* nonBondedPtrs = pAmber->getNonBondedPtrs();
    int* nonBonded14Ptrs = pAmber->getNonBonded14Ptrs();

    int *nonBonded = pAmber->getNonBonded();
    double *nonBondedParams = pAmber->getNonBondedParams();

    double nonBondedEnergy = 0;
    double atom1Charge = 0.0;
    double atom1_14Charge = 0.0;
    double atom2Charge = 0.0;
    double r = 0.0;
    double r2 = 0.0;
    double P6 = 0.0;
    double P12 = 0.0;
    double r6 = 0.0;
    double r12 = 0.0;
    double vdw_energy = 0.0;
    double ele_energy = 0.0;

    double R6 = 0.0;
    double R12 = 0.0;
    double vdW14Energy = 0.0;
    double ele14Energy= 0.0;
    double vdWEnergy = 0.0;
    double eleEnergy = 0.0;

    int atJ = 0;

    int c = 0;
    int cP = 0;

    for (int i = 0; i < nAtoms; i++) {
      atom1Charge = charges[i];
      atom1_14Charge = charges[i] * 0.83333; // scale factor, ELE == 1.2

      for (int j = 0; j < nonBondedPtrs[i]; j++) {
        atJ = nonBonded[c];
        if (atJ > -1) {
          r = sqrt ( pow((xyz[i*3]   - xyz[atJ*3  ]),2) +
                     pow((xyz[i*3+1] - xyz[atJ*3+1]),2) +
                     pow((xyz[i*3+2] - xyz[atJ*3+2]),2) );
          r2 = r*r;

          // van der Waals
          P6 = nonBondedParams[cP] / pow(r2, 3);   // sigma/r^6
          P12 = P6 * P6;
          r6 = nonBondedParams[cP+1] * (2.0 * P6); // epsilon * (2.0 sigma/r^6)
          r12 = nonBondedParams[cP+1] * P12;       // epsilon * (sigma^2/r^12)
          vdw_energy = r12 - r6;

          // Electrostatic
          atom2Charge = charges[atJ];
          ele_energy = atom2Charge / r;

          if (nonBonded14Ptrs[c]) {
            // van der Waals
            R6  -= 0.5 * r6; // scale factor, VDW == 2
            R12 += 0.5 * r12;
            vdW14Energy += 0.5 * vdw_energy;

            // Electrostatic
            ele14Energy += atom1_14Charge * ele_energy;
          }
          else {
            // van der Waals
            R6  -= r6;
            R12 += r12;
            vdWEnergy += vdw_energy;

            // Electrostatic
            eleEnergy += atom1Charge * ele_energy;
          }
        }
        c++;
        cP+=2;
      }
    }
    nonBondedEnergy = vdWEnergy + vdW14Energy + eleEnergy + ele14Energy;

    pAmber->setR6(R6);
    pAmber->setR12(R12);
    pAmber->setVDW14(vdW14Energy);
    pAmber->setVDW(vdWEnergy);
    pAmber->setEle14(ele14Energy);
    pAmber->setEle(eleEnergy);
    return nonBondedEnergy;
}
*/

// ============================================================
// Function : calculateG()
// ------------------------------------------------------------
//
// ============================================================
double amberNonBonded::calculateG()
{
    /*
    int nAtoms = pAmber->getNumAtoms();
    double *xyz = pAmber->getCoords();
    double *grad = pAmber->getGradients();
    double *charges = pAmber->getCharges();

    int* nonBondedPtrs = pAmber->getNonBondedPtrs();
    int* nonBonded14Ptrs = pAmber->getNonBonded14Ptrs();

    int *nonBonded = pAmber->getNonBonded();
    double *nonBondedParams = pAmber->getNonBondedParams();

    double nonBondedEnergy = 0;
    double atom1Charge = 0.0;
    double atom1_14Charge = 0.0;
    //double atom2Charge = 0.0;
    double r = 0.0;
    double r2 = 0.0;
    double P6 = 0.0;
    double P12 = 0.0;
    double r6 = 0.0;
    double r12 = 0.0;
    double vdw_energy = 0.0;
    double ele_energy = 0.0;

    double R6 = 0.0;
    double R12 = 0.0;
    double vdW14Energy = 0.0;
    double ele14Energy= 0.0;
    double vdWEnergy = 0.0;
    double eleEnergy = 0.0;

    double x12 = 0.0;
    double y12 = 0.0;
    double z12 = 0.0;

    double dX = 0.0;
    double dY = 0.0;
    double dZ = 0.0;

    double a1Charge = 0.0;
    double qiqj = 0.0;

    double eps = 0.0;
    double dEvdwdr = 0.0;
    double dEeledr = 0.0;

    int atJ = 0;

    int c = 0;
    int cP = 0;

    double dielectricConstant = 1.0;

    double ele14ScaleFactor = 1.2;
    double vdw14ScaleFactor = 2.0;
    double oneOverEle14ScaleFactor = 1.0/ele14ScaleFactor;
    double oneOverVdw14ScaleFactor = 1.0/vdw14ScaleFactor;

    for (int i = 0; i < nAtoms; i++) {
      atom1Charge = charges[i];
      atom1_14Charge = charges[i] * oneOverEle14ScaleFactor;
      for (int j = 0; j < nonBondedPtrs[i]; j++) {
        atJ = nonBonded[c];
        if (atJ > -1) {
          eps = 0.0;
          if (nonBonded14Ptrs[c]) {
            eps = oneOverVdw14ScaleFactor * nonBondedParams[cP+1];
            a1Charge = atom1_14Charge;
          }
          else {
            eps = nonBondedParams[cP+1];
            a1Charge = atom1Charge;
          }

          x12 = xyz[i*3  ] - xyz[atJ*3  ];
          y12 = xyz[i*3+1] - xyz[atJ*3+1];
          z12 = xyz[i*3+2] - xyz[atJ*3+2];

          r = sqrt(x12*x12 + y12*y12 + z12*z12);
          r2 = r*r;

          // van der Waals
          P6 = nonBondedParams[cP] / pow(r2, 3);
          P12 = P6 * P6;
          r6 = eps * (2.0 * P6);
          r12 = eps * P12;
          R6  -= r6;
          R12 += r12;

          vdw_energy = r12 - r6;

          dEvdwdr = (-12 * eps * (P12 - P6))/r;

          // Electrostatic
          qiqj = a1Charge * charges[atJ];
          ele_energy = qiqj / (dielectricConstant * r);
          dEeledr = - qiqj / (dielectricConstant * r2);

          dX = (dEvdwdr + dEeledr) * x12;
          dY = (dEvdwdr + dEeledr) * y12;
          dZ = (dEvdwdr + dEeledr) * z12;

          grad[i*3  ] += dX;
          grad[i*3+1] += dY;
          grad[i*3+2] += dZ;

          grad[atJ*3  ] -= dX;
          grad[atJ*3+1] -= dY;
          grad[atJ*3+2] -= dZ;

          if (nonBonded14Ptrs[c]) {
            vdW14Energy += vdw_energy; // van der Waals
            ele14Energy += ele_energy; // Electrostatic
          }
          else {
            vdWEnergy += vdw_energy; // van der Waals
            eleEnergy += ele_energy; // Electrostatic
          }
        }
        c++;
        cP+=2;
      }
    }
    nonBondedEnergy = vdWEnergy + vdW14Energy + eleEnergy + ele14Energy;

    pAmber->setR6(R6);
    pAmber->setR12(R12);
    pAmber->setVDW14(vdW14Energy);
    pAmber->setVDW(vdWEnergy);
    pAmber->setEle14(ele14Energy);
    pAmber->setEle(eleEnergy);
    return nonBondedEnergy;
    */
    return 0.0;
}

// ============================================================
// Function : decompose()
// ------------------------------------------------------------
//
// ============================================================
double amberNonBonded::decompose()
{
    int nAtoms = pAmber->getNumAtoms();
    int nUniqueTypes = pAmber->getNumTypes();
    int* types = pAmber->getIntTypes();
    int* nonBondedParameterIndex = pAmber->getNonBondedParameterIndex();
    int* outFlags = pAmber->getOutputFlags();

    double *xyz = pAmber->getCoords();
    double *charges = pAmber->getCharges();
    double *r6Params = pAmber->getR6Params();
    double *r12Params = pAmber->getR12Params();

    int *numExcluded = pAmber->getNumExcluded();
    int *excluded = pAmber->getExcluded();

    int *numExcluded14 = pAmber->getNumExcluded14();
    int *excluded14 = pAmber->getExcluded14();

    int *atomFlags = pAmber->getAtomFlags();

    double nonBondedEnergy = 0;
    double atom1Charge = 0.0;
    double atom1_14Charge = 0.0;
    double atom2Charge = 0.0;
    double r = 0.0;
    double r2 = 0.0;
    double r6 = 0.0;
    double r12 = 0.0;
    double vdw_energy = 0.0;
    double ele_energy = 0.0;

    double R6 = 0.0;
    double R12 = 0.0;
    double vdW14Energy = 0.0;
    double ele14Energy= 0.0;
    double vdWEnergy = 0.0;
    double eleEnergy = 0.0;

    double r6Param = 0.0;
    double r12Param = 0.0;
    //int nij = 0;
    int excludedAtomCounter = 0;
    int excluded14AtomCounter = 0;
    for (int i = 0; i < nAtoms; i++) {
      int atomIType = types[i];
      setToOne(atomFlags, nAtoms);

      int nExcludedAtoms = numExcluded[i];
      for (int e = 0; e < nExcludedAtoms; e++) {
        atomFlags[excluded[excludedAtomCounter]] = 0;
        excludedAtomCounter++;
      }

      int nExcludedAtoms14 = numExcluded14[i];
      for (int e = 0; e < nExcludedAtoms14; e++) {
        atomFlags[excluded14[excluded14AtomCounter]] = -1;
        excluded14AtomCounter++;
      }

      atom1Charge = charges[i];
      atom1_14Charge = charges[i] * 0.83333; // scale factor, ELE == 1.2

      for (int j = i+1; j < nAtoms; j++) {
        if (atomFlags[j] == 1) {
          //nij++;
          int atomJType = types[j];
          int paramIndex = nonBondedParameterIndex[atomIType * nUniqueTypes + atomJType];

          r6Param = r6Params[paramIndex];   // 2 * epsilon * sigma
          r12Param = r12Params[paramIndex]; // epsilon *sigma^2
          r = sqrt ( pow((xyz[i*3  ] - xyz[j*3  ]),2) +
                     pow((xyz[i*3+1] - xyz[j*3+1]),2) +
                     pow((xyz[i*3+2] - xyz[j*3+2]),2) );
          r2 = r*r;

          // van der Waals
          r6 = r6Param / pow(r2, 3);
          r12 = r12Param / pow(r2, 6);
          vdw_energy = r12 - r6;

          // Electrostatic
          atom2Charge = charges[j];
          ele_energy = atom2Charge / r;

          // van der Waals
          R6  -= r6;
          R12 += r12;
          vdWEnergy += vdw_energy;

          // Electrostatic
          eleEnergy += atom1Charge * ele_energy;

          //if (outFlags[i] or outFlags[j]) {
          if ((outFlags[i]&&!outFlags[j]) or (!outFlags[i]&&outFlags[j])) {
/*
            gzprintf(pAmber->pwdGZFileStream, 
                "%7d%7d %10.7f %10.7f %10.7f %12.7f %10.4f\n",
                i+1, j+1, -r6, r12, vdw_energy, atom1Charge * ele_energy, r);
*/
            pAmber->printPWD(i+1, j+1, -r6, r12, vdw_energy, atom1Charge * ele_energy, r);
          }
        }
        else if (atomFlags[j] == -1) { // 1-4 Interaction
          //nij++;
          int atomJType = types[j];
          int paramIndex = nonBondedParameterIndex[atomIType * nUniqueTypes + atomJType];

          r6Param = r6Params[paramIndex];   // 2 * epsilon * sigma
          r12Param = r12Params[paramIndex]; // epsilon *sigma^2
          r = sqrt ( pow((xyz[i*3  ] - xyz[j*3  ]),2) +
                     pow((xyz[i*3+1] - xyz[j*3+1]),2) +
                     pow((xyz[i*3+2] - xyz[j*3+2]),2) );
          r2 = r*r;

          // van der Waals
          r6 = r6Param / pow(r2, 3);
          r12 = r12Param / pow(r2, 6);
          vdW14Energy += 0.5 * (r12 - r6);

          R6  -= 0.5 * r6;
          R12 += 0.5 * r12;

          // Electrostatic
          atom2Charge = charges[j];
          ele_energy = (atom2Charge / r);
          ele14Energy += atom1_14Charge * ele_energy;

          //if (outFlags[i] or outFlags[j]) {
          if ((outFlags[i]&&!outFlags[j]) or (!outFlags[i]&&outFlags[j])) {
/*
            gzprintf(pAmber->pwdGZFileStream, 
                "%7d%7d %10.7f %10.7f %10.7f %12.7f %10.4f\n",
                i+1, -(j+1), -0.5*r6, 0.5*r12, 0.5 * (r12 - r6), atom1_14Charge * ele_energy, r);
*/
            pAmber->printPWD(i+1, j+1, -r6, r12, vdw_energy, atom1Charge * ele_energy, r);
          }
        }
      }
    }

    //std::cout << " Total Number of Lines = " << nij << std::endl;

    nonBondedEnergy = vdWEnergy + eleEnergy + vdW14Energy + ele14Energy;

    pAmber->setR6(R6);
    pAmber->setR12(R12);
    pAmber->setVDW(vdWEnergy);
    pAmber->setEle(eleEnergy);

    pAmber->setVDW14(vdW14Energy);
    pAmber->setEle14(ele14Energy);

    return nonBondedEnergy;
}

// ============================================================
// Function : setToOne()
// ------------------------------------------------------------
//
// ============================================================
void amberNonBonded::setToOne(int v[], int size)
{
    for (int i = 0; i < size; i++) {
        v[i] = 1;
    }
}

} // MTKpp namespace

