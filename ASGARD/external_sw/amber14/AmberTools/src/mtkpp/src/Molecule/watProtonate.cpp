/*!
   \file watProtonate.cpp
   \brief Protonates water molecule

    This code was adapted from the gwh.f code in AMBER9
    Written by : Algorithm by : Shuichi Miyamoto, D.A. Case
               : Converted to C++ : Martin Peters

   $Date: 2010/08/19 11:33:30 $
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
#include "watProtonate.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "parameters.h"
#include "utility.h"

// Utils
#include "Utils/vector3d.h"
#include "Utils/constants.h"
#include "Utils/idObject.h"

#include "Log/errorHandler.h"
#include "Diagnostics/MTKException.h"

#include <sstream>
#include <math.h>
#include <cmath>

namespace MTKpp
{

// ============================================================
// Function : watProtonate()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
watProtonate::watProtonate()
{
    nSphereGridPts = 0;
    biggestNumberOfHits = 0;
    sphereGridCoords = 0;
    maxPairs = 0;
    Hs = 0;
    lonePairs = 0;
    this->setup();
}

// ============================================================
// Function : ~watProtonate()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
watProtonate::~watProtonate()
{
    delete sphereGridCoords;
    delete maxPairs;
    delete Hs;
    delete lonePairs;
}

// ============================================================
// Function : setup()
// ------------------------------------------------------------
//
// ============================================================
void watProtonate::setup()
{
    std::string eMess = "\n";
    int ijk = 0;
    // Create cubic grid
    //  If distance is less that or equal to that of a O-H bond store point
    //  Here we are just getting the total number of grid points for dynamic
    //  memory allocation
    for (int i = 1; i < 22; i++) {
      for (int j = 1; j < 22; j++) {
        for (int k = 1; k < 22; k++) {
          ijk = (i-11)*(i-11) + (j-11)*(j-11) + (k-11)*(k-11);
          if ((ijk > 81) and (ijk <= 100)) {
            nSphereGridPts++;
          }
        }
      }
    }

    eMess += " Number of sphere grid points = " + i2s(nSphereGridPts) + "\n";

    // Allocate memory
    try {
      sphereGridCoords = new int [nSphereGridPts*3];
      eMess += " Allocating sphereGridCoords: " + i2s(nSphereGridPts*3) + " ints \n";
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("watProtonate::run", "Memory Allocation Failure", MTK_ERROR);
      throw MTKException(" Memory Allocation Failure ");
      //exit(0);
    }

    // assign spherical grid coordinates
    int nsph = 0;
    for (int i = 1; i < 22; i++) {
      for (int j = 1; j < 22; j++) {
        for (int k = 1; k < 22; k++) {
          ijk = (i-11)*(i-11) + (j-11)*(j-11) + (k-11)*(k-11);
          if ((ijk > 81) and (ijk <= 100)) {
            sphereGridCoords[nsph*3  ] = i - 11;
            sphereGridCoords[nsph*3+1] = j - 11;
            sphereGridCoords[nsph*3+2] = k - 11;
            nsph++;
          }
        }
      }
    }

    // maxPairs --> For each grid point how many other grid points can be
    //              used to create a water molecule with the centroid
    try {
      maxPairs = new int [nSphereGridPts]; // -1
      eMess += " Allocating maxPairs: " + i2s(nSphereGridPts) + " ints \n";
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("watProtonate::run", "Memory Allocation Failure", MTK_ERROR);
      throw MTKException(" Memory Allocation Failure ");
      //exit(0);
    }

    for (int i = 0; i < nSphereGridPts; i++) {
      int nHits = 0;
      for (int j = i+1; j < nSphereGridPts; j++) {
        double kx = pow((sphereGridCoords[i*3  ] - sphereGridCoords[j*3  ]),2.0);
        double ky = pow((sphereGridCoords[i*3+1] - sphereGridCoords[j*3+1]),2.0);
        double kz = pow((sphereGridCoords[i*3+2] - sphereGridCoords[j*3+2]),2.0);
        double kxyz = kx + ky + kz;
        if ((kxyz > 213) and (kxyz < 243)) {
          nHits++;
        }
      }
      if (nHits > biggestNumberOfHits) biggestNumberOfHits = nHits;
      maxPairs[i] = nHits;
    }

    try {
      // Stores the matching grid points which form a pair of hydrogens
      Hs = new int [nSphereGridPts*biggestNumberOfHits];
      eMess += " Allocating Hs: " + i2s(nSphereGridPts*biggestNumberOfHits) + " ints \n";

      // For each pair of hydrogens we need to store the corresponding lone
      //  pairs
      lonePairs = new int [2*nSphereGridPts*biggestNumberOfHits];
      eMess += " Allocating lonePairs: " + i2s(2*nSphereGridPts*biggestNumberOfHits) + " ints";
    }
    catch (std::bad_alloc) {
      errorLogger.throwError("watProtonate::run", "Memory Allocation Failure", MTK_ERROR);
      //exit(0);
      throw MTKException(" Memory Allocation Failure ");
    }

    for (int i = 0; i < nSphereGridPts*biggestNumberOfHits; i++) {
     Hs[i] = 0;
    }

    for (int i = 0; i < 2*nSphereGridPts*biggestNumberOfHits; i++) {
      lonePairs[i] = 0;
    }

    for (int i = 0; i < nSphereGridPts; i++) {
      int nHits = 0;
      for (int j = i+1; j < nSphereGridPts; j++) {
        double kx = pow((sphereGridCoords[i*3  ] -
                         sphereGridCoords[j*3  ]),2.0);
        double ky = pow((sphereGridCoords[i*3+1] -
                         sphereGridCoords[j*3+1]),2.0);
        double kz = pow((sphereGridCoords[i*3+2] -
                         sphereGridCoords[j*3+2]),2.0);
        double kxyz = kx + ky + kz;
        if ((kxyz > 213) and (kxyz < 243)) {
          nHits++;
          Hs[i*biggestNumberOfHits + nHits] = j;
          vector3d a(sphereGridCoords[i*3  ]/10.0,
                     sphereGridCoords[i*3+1]/10.0,
                     sphereGridCoords[i*3+2]/10.0);
          vector3d b(sphereGridCoords[j*3  ]/10.0,
                     sphereGridCoords[j*3+1]/10.0,
                     sphereGridCoords[j*3+2]/10.0);
          vector3d c = a + b;
          vector3d z = cross(c,b);
          //vector3d y = cross(z,c);
          vector3d ux = c.unit();
          //vector3d uy = y.unit();
          vector3d uz = z.unit();
          vector3d lp1 (-ux[0] * 5.80 + uz[0] * 8.15,
                        -ux[1] * 5.80 + uz[1] * 8.15,
                        -ux[2] * 5.80 + uz[2] * 8.15);
          vector3d lp2 (-ux[0] * 5.80 - uz[2] * 8.15,
                        -ux[1] * 5.80 - uz[1] * 8.15,
                        -ux[2] * 5.80 - uz[2] * 8.15);
          double lp[6] = {lp1[0], lp1[1], lp1[2], lp2[0], lp2[1], lp2[2]};

          for (int p = 0; p < 2; p++) {
            double kmin = 99999;
            for (int n = 0; n < nSphereGridPts; n++) {
              double nx = pow((lp[p*3  ]- sphereGridCoords[n*3  ]),2.0);
              double ny = pow((lp[p*3+1]- sphereGridCoords[n*3+1]),2.0);
              double nz = pow((lp[p*3+2]- sphereGridCoords[n*3+2]),2.0);
              if (kmin > (nx+ny+nz)) {
                kmin = nx+ny+nz;
                try {
                  if (p*nSphereGridPts*biggestNumberOfHits +
                              i*biggestNumberOfHits + nHits < 2*nSphereGridPts*biggestNumberOfHits) {
                    lonePairs[p*nSphereGridPts*biggestNumberOfHits +
                              i*biggestNumberOfHits + nHits] = n;
                  }
                }
                catch (...) {
                  eMess += "\n lonePairs out-of-bounds";
                  errorLogger.throwError("watProtonate::setup", eMess, MTK_ERROR);
                  throw MTKException(" out-of-bounds ");
                }
              }
            }
          }
        }
      }
    }
    errorLogger.throwError("watProtonate::setup", eMess, INFO);
}

// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
void watProtonate::run(collection* pCol)
{
    std::string eMess = "";
    std::vector<atom*> nonWaterAtoms;
    std::vector<atom*> waterOs;
    std::vector<molecule*> molList = pCol->getMoleculeList();
    for (unsigned int i = 0; i < molList.size(); i++) {
      std::string molName = molList[i]->getName();
      if (molName == "HOH" or molName == "WAT" or molName == "MOH") {
        std::vector<atom*> curAtList = molList[i]->getAtomList();
        if (curAtList.size() == 3) continue;
        if (molName == "MOH" and curAtList.size() == 2) continue;
        for (unsigned int j = 0; j < curAtList.size(); j++) {
          if (!curAtList[j]->getElement()) {
            eMess += "Atom has no element assigned";
            errorLogger.throwError("watProtonate::run", eMess, MTK_ERROR);
            //exit(0);
            throw MTKException(" watProtonate::run: Error: Atom has no element assigned ... exiting ");
          }
          if (curAtList[j]->getElement()->symbol == "O") {
            waterOs.push_back(curAtList[j]);
          }
        }
      }
      else {
        std::vector<atom*> curAtList = molList[i]->getAtomList();
        for (unsigned int j = 0; j < curAtList.size(); j++) {
          nonWaterAtoms.push_back(curAtList[j]);
        }
      }
    }

    unsigned int nWaters = waterOs.size();

    eMess += "\nwatProtonate: Total Number of Water Molecules = " + i2s(nWaters);
    eMess += "\nwatProtonate: Total Number of Non Water Atoms = " + i2s(nonWaterAtoms.size()) + "\n";

    if (nWaters == 0) {
      return;
    }

    std::vector<vector3d*> newHcoords;
    for (unsigned int i = 0; i < nWaters*2; i++) {
      vector3d* newCoord = new vector3d(0.0);
      newHcoords.push_back(newCoord);
    }

    double* gridEnergiesP;
    double* gridEnergiesW;

    // Calculate the electrostatic energy for each grid point on all the
    // water oxygens
    try {
      // gridEnergiesP
      gridEnergiesP = new double [nWaters * nSphereGridPts];

      // gridEnergiesW
      gridEnergiesW = new double [nWaters * nSphereGridPts];

      eMess += " Allocating Memory: gridEnergiesP " + i2s(nWaters * nSphereGridPts) + "\n";
      eMess += " Allocating Memory: gridEnergiesW " + i2s(nWaters * nSphereGridPts) + "\n";

    }
    catch (std::bad_alloc) {
      errorLogger.throwError("watProtonate::run", "Memory Allocation Failure", MTK_ERROR);
      //exit(0);
      std::stringstream ss;
      ss << " watProtonate::Memory Allocation Failure ... exiting " << std::endl;
      std::cout << ss.str();
      throw MTKException(ss.str());
    }

    vector3d gridPt(0.0);
    double r2 = 0.0;
    double gridPtEnergy = 0.0;
    double enemin = BIGNUM;
    double stdAccCharge = 0.0;

    std::vector<idObject*> watEnergies;
    for (unsigned int i = 0; i < nWaters; i++) {
      eMess += "\nWater: " + i2s(i+1);
      enemin = BIGNUM;
      vector3d waterOcoords = *(waterOs[i]->getCoords());
      for (int j = 0; j < nSphereGridPts; j++) {
        gridPt.setX(double(sphereGridCoords[j*3  ] / 10.0) + waterOcoords[0]);
        gridPt.setY(double(sphereGridCoords[j*3+1] / 10.0) + waterOcoords[1]);
        gridPt.setZ(double(sphereGridCoords[j*3+2] / 10.0) + waterOcoords[2]);

        gridPtEnergy = 0.0;
        for (unsigned int k = 0; k < nonWaterAtoms.size(); k++) {
          stdAtom* pStdAcc = nonWaterAtoms[k]->getStdAtom();
          if (pStdAcc) {
            //std::cout << pStdAcc->atmCharge << " " << E2KCAL << " " << pStdAcc->atmCharge * E2KCAL << std::endl;
            double lCharge = 1.0;
            if (nonWaterAtoms[k]->getElement()->group > 2 and
                nonWaterAtoms[k]->getElement()->group < 13) lCharge = 5.0;
            stdAccCharge = lCharge * pStdAcc->atmCharge * E2KCAL;
          }
          else {
            eMess += "\n" + nonWaterAtoms[k]->getName() + " charge is zero";
            stdAccCharge = 0.0;
          }
          r2 = gridPt.dist(*(nonWaterAtoms[k]->getCoords()));
          gridPtEnergy += stdAccCharge / r2;
        }
        gridEnergiesP[i * nSphereGridPts + j] = gridPtEnergy;

        if (enemin > gridPtEnergy) enemin = gridPtEnergy;
      }

      idObject* idO = new idObject(i, enemin);
      watEnergies.push_back(idO);
      eMess += " E = " + d2s(enemin);
    }
    errorLogger.throwError("watProtonate::run", eMess, INFO);

    // sort watEnergies in ascending order
    std::sort(watEnergies.begin(), watEnergies.end(), idObject::less);
    errorLogger.throwError("watProtonate::run", "sorted water energies", INFO);

    double weigh = -1.0;

    unsigned int nCycles = 10;
    double gridPt_water2H1 = 0.0;
    double gridPt_water2H2 = 0.0;

    eMess = "";
    for (unsigned int i = 0; i < nCycles; i++) { // 230
      double etot = 0.0;
      for (unsigned int j = 0; j < nWaters; j++) { // 210
        int curWater = watEnergies[j]->getI();
        vector3d waterOcoords = *(waterOs[curWater]->getCoords());
        if (j != 0) { // 170
          for (int k = 0; k < nSphereGridPts; k++) { // 160
            gridPt.setX(double(sphereGridCoords[k*3  ]/10.0)+waterOcoords[0]);
            gridPt.setY(double(sphereGridCoords[k*3+1]/10.0)+waterOcoords[1]);
            gridPt.setZ(double(sphereGridCoords[k*3+2]/10.0)+waterOcoords[2]);
            double esum = 0.0;
            for (unsigned int j2 = 0; j2 < nWaters; j2++) { // 140
              if (i == 0 and j == j2) break;
              if (j == j2) continue;
              int curWater2 = watEnergies[j2]->getI();
              vector3d water2Ocoords = *(waterOs[curWater2]->getCoords());
              r2 = gridPt.dist(water2Ocoords);
              esum += (-0.834 * E2KCAL) / r2;

              // Calc elec. energy between gridPt of water1 and Hs on water2
              vector3d* water2_H1 = newHcoords[curWater2*2  ];
              vector3d* water2_H2 = newHcoords[curWater2*2+1];
              gridPt_water2H1 = gridPt.dist(*(water2_H1));
              gridPt_water2H2 = gridPt.dist(*(water2_H2));
              esum += (0.417 * E2KCAL) / gridPt_water2H1;
              esum += (0.417 * E2KCAL) / gridPt_water2H2;
            } // 140
            gridEnergiesW[j * nSphereGridPts + k] = esum;
          } // 160
        } // 170

        double emint = BIGNUM;
        int nmini = 0;
        int nminj = 0;
        int nmink = 0;
        for (int k = 0; k < nSphereGridPts; k++) { // 190
          if (maxPairs[k] == 0) continue;
          double eminj = BIGNUM;
          for (int l = 0; l < maxPairs[k]; l++) { // 180
            int lp1 = lonePairs[0*nSphereGridPts*biggestNumberOfHits +
                                k*biggestNumberOfHits + l];
            int lp2 = lonePairs[1*nSphereGridPts*biggestNumberOfHits +
                                k*biggestNumberOfHits + l];
            double elp =  gridEnergiesP[lp1] + gridEnergiesP[lp2] +
                          gridEnergiesW[lp1] + gridEnergiesW[lp2];

            if (eminj > gridEnergiesP[Hs[k*biggestNumberOfHits + l]] +
                        gridEnergiesW[Hs[k*biggestNumberOfHits + l]] -
                        (elp * weigh)) {
              eminj = gridEnergiesP[Hs[k*biggestNumberOfHits + l]] +
                      gridEnergiesW[Hs[k*biggestNumberOfHits + l]] -
                      (elp * weigh);
              nminj = Hs[k*biggestNumberOfHits + l];
            }
          } // 180
          if (emint > gridEnergiesP[curWater * nSphereGridPts + k] +
                      gridEnergiesW[curWater * nSphereGridPts + k] + eminj) {
            emint = gridEnergiesP[curWater * nSphereGridPts + k] +
                    gridEnergiesW[curWater * nSphereGridPts + k] + eminj;
            nmini = k;
            nmink = nminj;
          }
        } // 190

        etot += emint;

        if (nmini == nmink) {
          std::cout << " 1 Error in watProtonate ... exiting " << std::endl;
          //exit(0);
          throw MTKException(" Error in watProtonate ... exiting ");
        }

        vector3d vH1(0.0);
        vector3d vH2(0.0);
        vH1.setX(double(sphereGridCoords[nmini*3  ] / 10.0));
        vH1.setY(double(sphereGridCoords[nmini*3+1] / 10.0));
        vH1.setZ(double(sphereGridCoords[nmini*3+2] / 10.0));
        vH2.setX(double(sphereGridCoords[nmink*3  ] / 10.0));
        vH2.setY(double(sphereGridCoords[nmink*3+1] / 10.0));
        vH2.setZ(double(sphereGridCoords[nmink*3+2] / 10.0));
        vector3d c = vH1 + vH2;
        vector3d z = cross(c,vH2);
        //vector3d y = cross(z,c);
        vector3d ux = c.unit();
        //vector3d uy = y.unit();
        vector3d uz = z.unit();
        vector3d H1 (waterOcoords[0] + ux[0] * 0.586077 + uz[0] * 0.7568,
                     waterOcoords[1] + ux[1] * 0.586077 + uz[1] * 0.7568,
                     waterOcoords[2] + ux[2] * 0.586077 + uz[2] * 0.7568);
        vector3d H2 (waterOcoords[0] + ux[0] * 0.586077 - uz[0] * 0.7568,
                     waterOcoords[1] + ux[1] * 0.586077 - uz[1] * 0.7568,
                     waterOcoords[2] + ux[2] * 0.586077 - uz[2] * 0.7568);

        vector3d* aH1 = newHcoords[curWater*2  ];
        aH1->setX(H1[0]);
        aH1->setY(H1[1]);
        aH1->setZ(H1[2]);
        vector3d* aH2 = newHcoords[curWater*2+1];
        aH2->setX(H2[0]);
        aH2->setY(H2[1]);
        aH2->setZ(H2[2]);
      } // 210
      eMess += "\nnCycle: " + i2s(i+1) + " Average E = " + d2s(etot/double(nWaters));
    } // 230

    errorLogger.throwError("watProtonate::run", eMess, INFO);

    // Add H atoms to water molecules
    for (unsigned int i = 0; i < waterOs.size(); i++) {
      //std::cout << " Adding Hs to water " << i+1 << " ";
      submolecule* waterSMol = waterOs[i]->getParent();
      molecule* waterMol = waterSMol->getParent();
      //std::cout << waterMol->getName() << " " << waterSMol->getSubMolId() << std::endl;
      vector3d* vH1 =  newHcoords[i*2  ];
      vector3d* vH2 =  newHcoords[i*2+1];
      int toAdd[2] = {1,1};

      if (waterMol->getName() == "MOH") {
        double dHEnergy1 = 0.0;
        double dHEnergy2 = 0.0;

        for (unsigned int k = 0; k < nonWaterAtoms.size(); k++) {
          stdAtom* pStdAcc = nonWaterAtoms[k]->getStdAtom();
          if (pStdAcc) {
            stdAccCharge = pStdAcc->atmCharge * E2KCAL;
          }
          double r1 = vH1->dist(*(nonWaterAtoms[k]->getCoords()));
          r2 = vH2->dist(*(nonWaterAtoms[k]->getCoords()));
          dHEnergy1 += stdAccCharge / r1;
          dHEnergy2 += stdAccCharge / r2;
        }
        if (dHEnergy1 > dHEnergy2) {
          toAdd[0] = 0;
        }
        else {
          toAdd[1] = 0;
        }
/*
        std::cout << " MOH: " << dHEnergy1 << " " << dHEnergy2 << std::endl;
        std::cout << " MOH: " <<toAdd[0] << " " << toAdd[1] << std::endl;
*/
      }

      if (toAdd[0]) {
        atom* pAtomH1 = waterSMol->addAtom();
        pAtomH1->setElement(pCol->pElements->getElement("H"));
        pAtomH1->setName(" H1 ");
        pAtomH1->setCoords(vH1->getX(), vH1->getY(), vH1->getZ());

        // at1, at2, type, stereo, topology, length
        Bond* pBond1 = waterMol->addBond(waterOs[i], pAtomH1, 1, 0, 0, 0.0);
        if (!pBond1) {
          std::cout << " 2 Error in watProtonate ... exiting " << std::endl;
          //exit(0);
          throw MTKException(" Error in watProtonate ... exiting ");
        }
        waterOs[i]->addBondedAtom(pAtomH1);
        pAtomH1->addBondedAtom(waterOs[i]);

/*
        if (waterMol->getName() == "MOH") {
          std::cout << " MOH: " << pBond1->atom1->getName() << "-" << pBond1->atom2->getName() << std::endl;
        }
*/
      }

      if (toAdd[1]) {
        atom* pAtomH2 = waterSMol->addAtom();
        pAtomH2->setElement(pCol->pElements->getElement("H"));
        pAtomH2->setName(" H2 ");
        if (waterMol->getName() == "MOH") pAtomH2->setName(" H1 ");
        pAtomH2->setCoords(vH2->getX(), vH2->getY(), vH2->getZ());

        Bond* pBond2 = waterMol->addBond(waterOs[i], pAtomH2, 1, 0, 0, 0.0);
        if (!pBond2) {
          std::cout << " 3 Error in watProtonate ... exiting " << std::endl;
          //exit(0);
          throw MTKException(" Error in watProtonate ... exiting ");
        }
        waterOs[i]->addBondedAtom(pAtomH2);
        pAtomH2->addBondedAtom(waterOs[i]);
      }
    }

    newHcoords.erase(newHcoords.begin(), newHcoords.end());
    newHcoords.clear();

    watEnergies.erase(watEnergies.begin(), watEnergies.end());
    watEnergies.clear();

    delete gridEnergiesP;
    delete gridEnergiesW;
}

} // MTKpp namespace
