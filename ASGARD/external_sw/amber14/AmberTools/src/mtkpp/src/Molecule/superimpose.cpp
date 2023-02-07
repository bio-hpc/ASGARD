/*! 
   \file superimpose.cpp
   \brief Superimposes molecules
   \author Martin Peters

   $Date: 2010/05/04 20:53:54 $
   $Revision: 1.17 $

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

#include "superimpose.h"

// Molecule
#include "molecule.h"
#include "atom.h"
#include "Utils/vector3d.h"
#include "functionalize.h"

// Boost
//#include "Utils/diagonalize.h"

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

#include "Utils/diagonalize_eigen.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Class : superimpose()
// ------------------------------------------------------------
//  Constructor for the class.
// ============================================================
superimpose::superimpose()
{
    this->nAtoms = 0;
    this->dRMSD = 0.0;
    this->correspondenceType = -1;
}

// ============================================================
// Class : superimpose()
// ------------------------------------------------------------
//  Constructor for the class.
// ============================================================
superimpose::superimpose(int i)
{
    this->nAtoms = i;
    this->dRMSD = 0.0;
    this->correspondenceType = -1;
}

// ============================================================
// Function : ~superimpose()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
superimpose::~superimpose() {}

// ============================================================
// Function : fit()
// ------------------------------------------------------------
//
// ============================================================
double superimpose::fit(molecule* pMoleculeA, molecule* pMoleculeB)
{
    // - Check size of molecule A and B - //
    int nAtomsA = pMoleculeA->getNumAtoms();
    int nAtomsB = pMoleculeB->getNumAtoms();
    if (nAtomsA != nAtomsB) {
      std::string error = "ERROR = Molecule A contains: ";
      std::stringstream ss1;
      ss1 << nAtomsA;
      std::string str_error = ss1.str().c_str();
      std::string err2 = " Atoms, while Molecule B contains: ";
      std::stringstream ss2;
      ss2 << nAtomsB;
      std::string str_error2 = ss2.str().c_str();
      std::string err3 = " Atoms.";
      error = error + str_error + err2 + str_error2 + err3;
      throw MTKException(error);
    }
    else {
      nAtoms = nAtomsA;
    }
    double coordsA[nAtoms][3];
    double coordsB[nAtoms][3];

    pMoleculeA->getCoordinates(coordsA);
    pMoleculeB->getCoordinates(coordsB);
    dRMSD = this->fit(coordsA, coordsB, nAtoms);

    std::vector<atom*> atomListB = pMoleculeB->getAtomList();
    for (unsigned int i = 0; i < atomListB.size(); i++) {
      atomListB[i]->setCoords(coordsB[i][0], coordsB[i][1], coordsB[i][2]);
    }
    return dRMSD;
}

// ============================================================
// Function : fit()
// ------------------------------------------------------------
//
// ============================================================
double superimpose::fit(double coordsA[][3], double coordsB[][3], int n)
{
    this->nAtoms = n;

    //Eigen::VectorXd R(4);
    Eigen::MatrixXd quaternion(4,4);

    //ublas::vector<double> R(4);
    //ublas::matrix<double, ublas::column_major> quaternion(4,4);
    double itsDXM[nAtoms][3];
    double itsDXP[nAtoms][3];

    this->center(coordsA, centerA);
    this->center(coordsB, centerB);

    this->coordinateDiff(coordsA, coordsB, itsDXM, itsDXP);

    this->buildQuaternion(quaternion, itsDXM, itsDXP);

    SelfAdjointEigenSolver<MatrixXd> eigensolver(quaternion);
    MatrixXd evectors = eigensolver.eigenvectors();
    VectorXd R = eigensolver.eigenvalues();

    eigenValueSort(evectors, R, 1);

    this->buildRotation(rotMat, evectors);

/*
    // Boost
    // std::cout << " quaternion \n" << quaternion << std::endl;

    int r = diagonalize(quaternion,R);

    //std::cout << " quaternion \n" << quaternion << std::endl;
    //std::cout << " R \n" << R << std::endl;

    if (r != 0) {
      std::string errMessage = " Diagonalization Failed ";
      errorLogger.throwError("superimpose::fit", errMessage, MTK_ERROR);
      return -1.0;
    }

    eigenValueSort(quaternion,R,1);

    //std::cout << " After eigenValueSort " << std::endl;
    //std::cout << " quaternion \n" << quaternion << std::endl;
    //std::cout << " R \n" << R << std::endl;

    this->buildRotation(rotMat, quaternion);
*/

/*
    for (int j = 0; j < 3; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        std::cout << rotMat[i][j] << " ";
      }
      std::cout << " " << std::endl;
    }
*/
    // Move coordsB to the origin
    for (int j = 0; j < nAtoms; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        coordsB[j][i] = coordsB[j][i] - centerB[i];
      }
    }

    // rotate each coordinate of coordsB
    for (int j = 0; j < nAtoms; j++) {
      double tempCoord[3];
      tempCoord[0] = coordsB[j][0];
      tempCoord[1] = coordsB[j][1];
      tempCoord[2] = coordsB[j][2];
      coordsB[j][0] = 0.0; coordsB[j][1] = 0.0; coordsB[j][2] = 0.0;
      for (unsigned int a = 0; a < 3; a++) {
        for (unsigned int b = 0; b < 3; b++) {
          coordsB[j][a] = coordsB[j][a] + tempCoord[b] * rotMat[b][a];
        }
      }
    }

    // Move coordsB to coordsA frame of reference
    for (int j = 0; j < nAtoms; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        coordsB[j][i] = coordsB[j][i] + centerA[i];
      }
    }

    dRMSD = calculateRMSD(coordsA,coordsB);

    return dRMSD;
}

// ============================================================
// Function : fit()
// ------------------------------------------------------------
//
// ============================================================
int superimpose::fit(molecule* pMoleculeA, molecule* pMoleculeB,
                     std::vector<int> cor, int cT)
{
    this->correspondenceType = cT;
    double lRotMat[3][3];

    int lnAtomsA = pMoleculeA->getNumAtoms();
    int lnAtomsB = pMoleculeB->getNumAtoms();

    double allCoordsB[lnAtomsB][3];

    std::vector< vector3d > molACoords;
    pMoleculeA->getCoordinates(molACoords);

    std::vector< vector3d > molBCoords;
    pMoleculeB->getCoordinates(molBCoords);

    if (lnAtomsA == lnAtomsB) {
      nAtoms = lnAtomsA;
    }
    else {
      std::cout << " Error in superimpose ... exiting " << std::endl;
      return 1;
    }

    for (int i = 0; i < lnAtomsB; i++) {
      for (int j = 0; j < 3; j++) {
        allCoordsB[i][j] = molBCoords[i][j];
      }
    }

    if (this->correspondenceType == 2 or this->correspondenceType == 3) {
      this->nHeavyAtomsA = pMoleculeA->getNumHeavyAtoms();
      this->heavyAtomIndicesA = pMoleculeA->getHeavyAtomIndices();
      this->nHeavyAtomsB = pMoleculeB->getNumHeavyAtoms();
      this->heavyAtomIndicesB = pMoleculeB->getHeavyAtomIndices();
      if (this->nHeavyAtomsA != this->nHeavyAtomsB) {
        throw MTKException("ERROR4");
      }
      nAtoms = nHeavyAtomsA;
    }

    double coordsA[nAtoms][3];
    double coordsB[nAtoms][3];
    double cB[nAtoms][3];

    if (this->correspondenceType == 0 or this->correspondenceType == 1) {
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          coordsA[i][j] = molACoords[i][j];
          coordsB[i][j] = molBCoords[i][j];
        }
      }
    }
    else if (this->correspondenceType == 2 or this->correspondenceType == 3) {
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          coordsA[i][j] = molACoords[heavyAtomIndicesA[i]][j];
          coordsB[i][j] = molBCoords[heavyAtomIndicesB[i]][j];
        }
      }
    }

    for (int t = 0; t < nAtoms; t++) {
      for (int t2 = 0; t2 < nAtoms; t2++) {
        if (cor[t*nAtoms+t2]) {
          cB[t][0] = coordsB[t2][0];
          cB[t][1] = coordsB[t2][1];
          cB[t][2] = coordsB[t2][2];
        }
      }
    }

    //Eigen::VectorXd R(4);
    Eigen::MatrixXd quaternion(4,4);

    //ublas::vector<double> R(4);
    //ublas::matrix<double, ublas::column_major> quaternion(4,4);
    double itsDXM[nAtoms][3];
    double itsDXP[nAtoms][3];

    this->center(coordsA, centerA);
    this->center(cB, centerB);

    this->coordinateDiff(coordsA, cB, itsDXM, itsDXP);

    this->buildQuaternion(quaternion, itsDXM, itsDXP);

    SelfAdjointEigenSolver<MatrixXd> eigensolver(quaternion);
    MatrixXd evectors = eigensolver.eigenvectors();
    VectorXd R = eigensolver.eigenvalues();

    eigenValueSort(evectors, R, 1);

    this->buildRotation(rotMat, evectors);

/*
    // Boost
    int r = diagonalize(quaternion,R);

    if (r != 0) return 1;

    eigenValueSort(quaternion,R,1);

    this->buildRotation(lRotMat,quaternion);
*/

/*
    for (int j = 0; j < 3; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        std::cout << lRotMat[i][j] << " ";
      }
      std::cout << " " << std::endl;
    }
    this->updateCoords(cB, nAtoms, centerB, centerA, lRotMat);
    double asdfdRMSD = calculateRMSD(coordsA,cB);
*/

    this->updateCoords(allCoordsB, lnAtomsB, centerB, centerA, lRotMat);

    std::vector<atom*> atomListB = pMoleculeB->getAtomList();
    for (unsigned int i = 0; i < atomListB.size(); i++) {
      atomListB[i]->setCoords(allCoordsB[i][0], allCoordsB[i][1], allCoordsB[i][2]);
    }

    // Set centers and rotation matrix for later use
    for (unsigned int a = 0; a < 3; a++) {
      for (unsigned int b = 0; b < 3; b++) {
        rotMat[b][a] = lRotMat[b][a];
      }
    }

    return 0;
}

// ============================================================
// Function : rsmd()
// ------------------------------------------------------------
// Calculate the coordinate rmsd
// ============================================================
double superimpose::rmsd(molecule* pMoleculeA, molecule* pMoleculeB)
{
    // - Check size of molecule A and B - //
    int nAtomsA = pMoleculeA->getNumAtoms();
    int nAtomsB = pMoleculeB->getNumAtoms();
    if (nAtomsA != nAtomsB) {
      std::string error = "ERROR = Molecule A contains: ";
      std::stringstream ss1;
      ss1 << nAtomsA;
      std::string str_error = ss1.str().c_str();
      std::string err2 = " Atoms, while Molecule B contains: ";
      std::stringstream ss2;
      ss2 << nAtomsB;
      std::string str_error2 = ss2.str().c_str();
      std::string err3 = " Atoms.";
      error = error + str_error + err2 + str_error2 + err3;
      throw MTKException(error);
    }
    else {
      nAtoms = nAtomsA;
    }
    double coordsA[nAtoms][3];
    double coordsB[nAtoms][3];

    pMoleculeA->getCoordinates(coordsA);
    pMoleculeB->getCoordinates(coordsB);
    dRMSD = this->fit(coordsA, coordsB, nAtoms);
    return dRMSD;
}

// ============================================================
// Function : rsmd()
// ------------------------------------------------------------
// Calculate the coordinate rmsd
// ============================================================
double superimpose::rmsd(molecule* pMoleculeA,
                         std::vector< vector3d > &molBCoords,
                         std::vector<std::vector<int> > &correspondenceMatrices,
                         int &cor)
{
    // - Check size of molecule A and B - //
    int nAtomsA = pMoleculeA->getNumAtoms();
    if (nAtomsA != static_cast<int>(molBCoords.size())) {
      std::cout << " superimpose::rmsd Error [number of atoms don't match] A:"
                << nAtomsA << " B:" << molBCoords.size() << " ... exiting " << std::endl;
      throw MTKException(" Error in superimpose ... exiting ");
    }
    else {
      nAtoms = nAtomsA;
    }

    double minRMSD = BIGNUM;
    double lRMSD = 0.0;

    if (this->correspondenceType == 2 or this->correspondenceType == 3) {
      nAtoms = nHeavyAtomsA;
    }

    double coordsA[nAtoms][3];
    double coordsB[nAtoms][3];
    double cB[nAtoms][3];
    std::vector<int> correspondenceMatrix;

    std::vector< vector3d > molACoords;
    pMoleculeA->getCoordinates(molACoords);

    if (this->correspondenceType == 0 or this->correspondenceType == 1) {
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          coordsA[i][j] = molACoords[i][j];
          coordsB[i][j] = molBCoords[i][j];
        }
      }
    }
    else if (this->correspondenceType == 2 or this->correspondenceType == 3) {
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          coordsA[i][j] = molACoords[heavyAtomIndicesA[i]][j];
          coordsB[i][j] = molBCoords[heavyAtomIndicesB[i]][j];
        }
      }
    }

    for (unsigned int k = 0; k < correspondenceMatrices.size(); k++) {
      correspondenceMatrix = correspondenceMatrices[k];
      int sdf = 0;
      for (int t = 0; t < nAtoms; t++) {
        //for (int t2 = t; t2 < nAtoms; t2++)
        for (int t2 = 0; t2 < nAtoms; t2++) {
          if (correspondenceMatrix[t*nAtoms+t2]) {
            cB[t][0] = coordsB[t2][0];
            cB[t][1] = coordsB[t2][1];
            cB[t][2] = coordsB[t2][2];
            sdf++;
          }
        }
      }
      if (sdf != nAtoms) {
        throw MTKException(" Error in superimpose ... exiting ");
      }

      lRMSD = this->fit(coordsA, cB, nAtoms);
      if (lRMSD < 0.0) continue;
      if (lRMSD < minRMSD) {
        minRMSD = lRMSD;
        cor = k;
      }
    }

    if (minRMSD < 0.000001) {
      if (minRMSD < 0.0) {
        minRMSD = -1.0;
      }
      else {
        minRMSD = 0.0;
      }
    }
    return minRMSD;
}

// ============================================================
// Function : rsmd()
// ------------------------------------------------------------
// Calculate the coordinate rmsd
// ============================================================
double superimpose::rmsd(molecule* pMoleculeA, molecule* pMoleculeB,
                    std::vector<std::vector<int> > &correspondenceMatrices)
{
    double minRMSD = BIGNUM;
/*
    double lRMSD = 0.0;

    // - Check size of molecule A and B - //
    int nAtomsA = pMoleculeA->getNumAtoms();
    int nAtomsB = pMoleculeB->getNumAtoms();
    if (nAtomsA != nAtomsB) {
      std::string error = "ERROR = Molecule A contains: ";
      std::stringstream ss1;
      ss1 << nAtomsA;
      std::string str_error = ss1.str().c_str();
      std::string err2 = " Atoms, while Molecule B contains: ";
      std::stringstream ss2;
      ss2 << nAtomsB;
      std::string str_error2 = ss2.str().c_str();
      std::string err3 = " Atoms.";
      error = error + str_error + err2 + str_error2 + err3;
      std::cout << error << std::endl;
      throw MTKException(error);
    }
    else {
      nAtoms = nAtomsA;
    }

    double coordsA[nAtoms][3];
    pMoleculeA->getCoordinates(coordsA);

    double coordsB[nAtoms][3];
    pMoleculeB->getCoordinates(coordsB);

    double cB[nAtoms][3];

    std::vector<int> correspondenceMatrix;

    for (unsigned int k = 0; k < correspondenceMatrices.size(); k++) {
      correspondenceMatrix = correspondenceMatrices[k];
      for (int t = 0; t < nAtoms; t++) {
        for (int t2 = t; t2 < nAtoms; t2++) {
          if (correspondenceMatrix[t*nAtoms+t2]) {
            cB[t][0] = coordsB[t2][0];
            cB[t][1] = coordsB[t2][1];
            cB[t][2] = coordsB[t2][2];
          }
        }
      }

      lRMSD = this->fit(coordsA, cB, nAtoms);
      //std::cout << " Current rmsd = " << lRMSD << std::endl;
      if (lRMSD < minRMSD) minRMSD = lRMSD;
    }
    if (minRMSD < 0.000001) {
      minRMSD = 0.0;
    }
*/
    return minRMSD;
}

// ============================================================
// Function : rmsdNoFit()
// ------------------------------------------------------------
// Calculate the coordinate rmsd
// ============================================================
double superimpose::rmsdNoFit(molecule* pMoleculeA,
                    std::vector< vector3d > &molBCoords,
                    std::vector<std::vector<int> > &correspondenceMatrices)
{
    double minRMSD = BIGNUM;
    double lRMSD = 0.0;

    double coordsA[this->nHeavyAtomsA][3];
    double coordsB[this->nHeavyAtomsA][3];
    double cB[this->nHeavyAtomsA][3];
    std::vector<int> correspondenceMatrix;

    std::vector< vector3d > molACoords;
    pMoleculeA->getCoordinates(molACoords);

    for (int i = 0; i < this->nHeavyAtomsA; i++) {
      for (int j = 0; j < 3; j++) {
        coordsA[i][j] = molACoords[heavyAtomIndicesA[i]][j];
        coordsB[i][j] = molBCoords[heavyAtomIndicesB[i]][j];
      }
    }

    nAtoms = this->nHeavyAtomsA;

    for (unsigned int k = 0; k < correspondenceMatrices.size(); k++) {
      correspondenceMatrix = correspondenceMatrices[k];
      for (int t = 0; t < this->nHeavyAtomsA; t++) {
        for (int t2 = 0; t2 < this->nHeavyAtomsA; t2++) {
          if (correspondenceMatrix[t*this->nHeavyAtomsA+t2]) {
            cB[t][0] = coordsB[t2][0];
            cB[t][1] = coordsB[t2][1];
            cB[t][2] = coordsB[t2][2];
          }
        }
      }

      lRMSD = this->calculateRMSD(coordsA, cB);
      if (lRMSD < minRMSD) minRMSD = lRMSD;
    }
    if (minRMSD < 0.000001) {
      minRMSD = 0.0;
    }
    return minRMSD;
}

// ============================================================
// Function : getRotationMatrix()
// ------------------------------------------------------------
//
// ============================================================
void superimpose::getRotationMatrix(double lCenterA[3], double lCenterB[3],
                                   double lRotMat[][3])
{
    for (unsigned int a = 0; a < 3; a++) {
      lCenterA[a] = centerA[a];
      lCenterB[a] = centerB[a];
      for (unsigned int b = 0; b < 3; b++) {
        lRotMat[b][a] = rotMat[b][a];
      }
    }
}

// ============================================================
// Function : getRotationMatrix()
// ------------------------------------------------------------
//
// ============================================================
int superimpose::getRotationMatrix(double coordsA[][3], double coordsB[][3],
                                   int n, double lRotMat[][3])
{
    this->nAtoms = n;

    //Eigen::VectorXd R(4);
    Eigen::MatrixXd quaternion(4,4);

    //ublas::vector<double> R(4);
    //ublas::matrix<double, ublas::column_major> quaternion(4,4);
    double itsDXM[nAtoms][3];
    double itsDXP[nAtoms][3];

    this->center(coordsA, centerA);
    this->center(coordsB, centerB);

    this->coordinateDiff(coordsA, coordsB, itsDXM, itsDXP);

    this->buildQuaternion(quaternion, itsDXM, itsDXP);

    SelfAdjointEigenSolver<MatrixXd> eigensolver(quaternion);
    MatrixXd evectors = eigensolver.eigenvectors();
    VectorXd R = eigensolver.eigenvalues();

    eigenValueSort(evectors, R, 1);

    this->buildRotation(rotMat, evectors);

/*
    // Boost
    int r = diagonalize(quaternion,R);

    if (r != 0) return 1;

    eigenValueSort(quaternion,R,1);

    this->buildRotation(lRotMat,quaternion);
*/

    // Move coordsB to the origin
    for (int j = 0; j < nAtoms; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        coordsB[j][i] = coordsB[j][i] - centerB[i];
      }
    }

    // rotate each coordinate of coordsB
    for (int j = 0; j < nAtoms; j++) {
      double tempCoord[3];
      tempCoord[0] = coordsB[j][0];
      tempCoord[1] = coordsB[j][1];
      tempCoord[2] = coordsB[j][2];
      coordsB[j][0] = 0.0; coordsB[j][1] = 0.0; coordsB[j][2] = 0.0;
      for (unsigned int a = 0; a < 3; a++) {
        for (unsigned int b = 0; b < 3; b++) {
          coordsB[j][a] = coordsB[j][a] + tempCoord[b] * lRotMat[b][a];
        }
      }
    }

    // Move coordsB to coordsA frame of reference
    for (int j = 0; j < nAtoms; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        coordsB[j][i] = coordsB[j][i] + centerA[i];
      }
    }

    return 0;
}

// ============================================================
// Function : getRotationMatrix()
// ------------------------------------------------------------
//
// ============================================================
int superimpose::getRotationMatrix(double coordsA[][3], double coordsB[][3],
                                   double lCenterA[3], double lCenterB[3],
                                   int n, double lRotMat[][3])
{
    this->nAtoms = n;

    Eigen::MatrixXd quaternion(4,4);

    //ublas::vector<double> R(4);
    //ublas::matrix<double, ublas::column_major> quaternion(4,4);
    double itsDXM[nAtoms][3];
    double itsDXP[nAtoms][3];

    for (int t = 0; t < 3; t++) {
      this->centerA[t] = lCenterA[t];
      this->centerB[t] = lCenterB[t];
    }

    this->coordinateDiff(coordsA, coordsB, itsDXM, itsDXP);

    this->buildQuaternion(quaternion, itsDXM, itsDXP);
/*
    // Boost
    int r = diagonalize(quaternion,R);

    if (r != 0) return 1;

    eigenValueSort(quaternion,R,1);

    this->buildRotation(lRotMat,quaternion);
*/

    SelfAdjointEigenSolver<MatrixXd> eigensolver(quaternion);
    MatrixXd evectors = eigensolver.eigenvectors();
    VectorXd R = eigensolver.eigenvalues();

    eigenValueSort(evectors, R, 1);

    this->buildRotation(rotMat, evectors);

    return 0;
}

// ============================================================
// Function : updateCoords()
// ------------------------------------------------------------
//
// ============================================================
void superimpose::updateCoords(double coords[][3], int n, double c1[3],
                 double c2[3], double lRotMat[][3])
{
    // Move coords to the origin
    for (int j = 0; j < n; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        coords[j][i] = coords[j][i] - c1[i];
      }
    }

    // rotate each coordinate of coords
    for (int j = 0; j < n; j++) {
      double tempCoord[3];
      tempCoord[0] = coords[j][0];
      tempCoord[1] = coords[j][1];
      tempCoord[2] = coords[j][2];
      coords[j][0] = 0.0; coords[j][1] = 0.0; coords[j][2] = 0.0;
      for (unsigned int a = 0; a < 3; a++) {
        for (unsigned int b = 0; b < 3; b++) {
          coords[j][a] = coords[j][a] + tempCoord[b] * lRotMat[b][a];
        }
      }
    }

    // Move coordsB to coordsA frame of reference
    for (int j = 0; j < n; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        coords[j][i] = coords[j][i] + c2[i];
      }
    }
}

// ============================================================
// Function : center()
// ------------------------------------------------------------
// Calculate the center of mass of a coordinate set
// ============================================================
void superimpose::center(double coords[][3], double center[3])
{

    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;

    for (int j = 0; j < nAtoms; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        center[i] = center[i] + coords[j][i];
      }
    }

    //std::cout << " center: ";
    for (unsigned int i = 0; i < 3; i++) {
      center[i] = center[i] / nAtoms;
      //std::cout << center[i] << " ";
    }
    //std::cout << " " << std::endl;
}

// ============================================================
// Function : coordinateDiff()
// ------------------------------------------------------------
// Calculates Coordinate Differences itsDXM, itsDXP.
// ============================================================
void superimpose::coordinateDiff(double coordsA[][3], double coordsB[][3],
                                 double itsDXM[][3], double itsDXP[][3])
{
    for (int j = 0; j < nAtoms; j++) {
      for (int i = 0; i < 3; i++) {
        itsDXM[j][i] =  coordsB[j][i] - this->centerB[i] -
                       (coordsA[j][i] - this->centerA[i]);

        itsDXP[j][i] =  coordsB[j][i] - this->centerB[i] +
                       (coordsA[j][i] - this->centerA[i]);
      }
/*
#ifdef DEBUG
      std::cout << " itsDXM : " << itsDXM[j][0] << " " << itsDXM[j][1] << " " << itsDXM[j][2] << std::endl;
      std::cout << " itsDXP : " << itsDXP[j][0] << " " << itsDXP[j][1] << " " << itsDXP[j][2] << std::endl;
#endif
*/
    }
}

// ============================================================
// Function : buildQuaternion()
// ------------------------------------------------------------
// Builds Quaternion Matrix
// ============================================================
/*
void superimpose::buildQuaternion(ublas::matrix<double, ublas::column_major>& quaternion,
                  double itsDXM[][3], double itsDXP[][3])
*/
void superimpose::buildQuaternion(Eigen::MatrixXd& quaternion,
                  double itsDXM[][3], double itsDXP[][3])

{
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        quaternion(i,j) = 0.0;
      }
    }

    for (int i = 0; i < nAtoms; i++) {
      quaternion(0,0) = quaternion(0,0) + itsDXM[i][0]*itsDXM[i][0] +
                        itsDXM[i][1]*itsDXM[i][1] + itsDXM[i][2]*itsDXM[i][2];
      quaternion(1,1) = quaternion(1,1) + itsDXP[i][1]*itsDXP[i][1] +
                        itsDXP[i][2]*itsDXP[i][2] + itsDXM[i][0]*itsDXM[i][0];
      quaternion(2,2) = quaternion(2,2) + itsDXP[i][0]*itsDXP[i][0] +
                        itsDXP[i][2]*itsDXP[i][2] + itsDXM[i][1]*itsDXM[i][1];
      quaternion(3,3) = quaternion(3,3) + itsDXP[i][0]*itsDXP[i][0] +
                        itsDXP[i][1]*itsDXP[i][1] + itsDXM[i][2]*itsDXM[i][2];

      quaternion(0,1) = quaternion(0,1) + itsDXP[i][1]*itsDXM[i][2] -
                        itsDXM[i][1]*itsDXP[i][2];
      quaternion(0,2) = quaternion(0,2) + itsDXM[i][0]*itsDXP[i][2] -
                        itsDXP[i][0]*itsDXM[i][2];
      quaternion(0,3) = quaternion(0,3) + itsDXP[i][0]*itsDXM[i][1] -
                        itsDXM[i][0]*itsDXP[i][1];
      quaternion(1,2) = quaternion(1,2) + itsDXM[i][0]*itsDXM[i][1] -
                        itsDXP[i][0]*itsDXP[i][1];
      quaternion(1,3) = quaternion(1,3) + itsDXM[i][0]*itsDXM[i][2] -
                        itsDXP[i][0]*itsDXP[i][2];
      quaternion(2,3) = quaternion(2,3) + itsDXM[i][1]*itsDXM[i][2] -
                        itsDXP[i][1]*itsDXP[i][2];
    }
    quaternion(1,0) = quaternion(0,1);
    quaternion(2,0) = quaternion(0,2);
    quaternion(3,0) = quaternion(0,3);
    quaternion(2,1) = quaternion(1,2);
    quaternion(3,1) = quaternion(1,3);
    quaternion(3,2) = quaternion(2,3);
}

// ============================================================
// Function : buildRotation()
// ------------------------------------------------------------
// Builds Rotation Matrix
// ============================================================
/*
void superimpose::buildRotation(double t[][3],
                  ublas::matrix<double, ublas::column_major> R)
*/
void superimpose::buildRotation(double t[][3],
                  Eigen::MatrixXd R)
{
     t[0][0] =      R(0,3) * R(0,3) + R(1,3) * R(1,3) -
                    R(2,3) * R(2,3) - R(3,3) * R(3,3);
     t[0][1] = 2 * (R(1,3) * R(2,3) + R(0,3) * R(3,3));
     t[0][2] = 2 * (R(1,3) * R(3,3) - R(0,3) * R(2,3));

     t[1][0] = 2 * (R(1,3) * R(2,3) - R(0,3) * R(3,3));
     t[1][1] =      R(0,3) * R(0,3) + R(2,3) * R(2,3) -
                    R(1,3) * R(1,3) - R(3,3) * R(3,3);
     t[1][2] = 2 * (R(2,3) * R(3,3) + R(0,3) * R(1,3));

     t[2][0] = 2 * (R(1,3) * R(3,3) + R(0,3) * R(2,3));
     t[2][1] = 2 * (R(2,3) * R(3,3) - R(0,3) * R(1,3));
     t[2][2] =      R(0,3) * R(0,3) + R(3,3) * R(3,3) -
                    R(1,3) * R(1,3) - R(2,3) * R(2,3);
}

// ============================================================
// Function : calculateRMSD()
// ------------------------------------------------------------
// Calculates the rmsd between two coordinate sets
// ============================================================
double superimpose::calculateRMSD(double coordsA[][3], double coordsB[][3])
{
//std::cout << "superimpose::calculateRMSD " << std::endl;
    dRMSD = 0.0;
    for (int j = 0; j < nAtoms; j++) {
/*
std::cout << coordsA[j][0] << ":" << coordsB[j][0] << " "
          << coordsA[j][1] << ":" << coordsB[j][1] << " "
          << coordsA[j][2] << ":" << coordsB[j][2] << std::endl;
*/
      dRMSD += ((coordsB[j][0] - coordsA[j][0])*(coordsB[j][0] - coordsA[j][0])
              + (coordsB[j][1] - coordsA[j][1])*(coordsB[j][1] - coordsA[j][1])
              + (coordsB[j][2] - coordsA[j][2])*(coordsB[j][2] - coordsA[j][2])
               );
    }
    if (dRMSD > 0.0) {
      dRMSD = sqrt( dRMSD / double(nAtoms));
    }
    else {
      dRMSD = 0.0;
    }

    //std::cout << " RMSD = " << dRMSD << std::endl;
    return dRMSD;
}

// ============================================================
// Function : minRMSD()
// ------------------------------------------------------------
// Calculates the rmsd between two coordinate sets
// ============================================================
double superimpose::minRMSD(molecule* pMoleculeA, molecule* pMoleculeB,
                            int type, int &nFittedAtoms)
{
    // Molecule A
    int nAtomsA = pMoleculeA->getNumAtoms();
    double coordsA[nAtomsA][3];
    pMoleculeA->getCoordinates(coordsA);

    // Molecule B
    int nAtomsB = pMoleculeB->getNumAtoms();
    double coordsB[nAtomsB][3];
    pMoleculeB->getCoordinates(coordsB);

    char* molA_Types = NULL;
    char* molB_Types = NULL;

    if (type == 0) {
      molA_Types = pMoleculeA->getAtomTypes();
      if (!molA_Types) {
        std::cout << " Error in superimpose " << std::endl;
        return 1;
      }
      molB_Types = pMoleculeB->getAtomTypes();
      if (!molB_Types) {
        std::cout << " Error in superimpose " << std::endl;
        return 1;
      }
    }
    else if (type == 1) {
      molA_Types = pMoleculeA->getAtomSymbols();
      if (!molA_Types) {
        std::cout << " Error in superimpose " << std::endl;
        return 1;
      }
      molB_Types = pMoleculeB->getAtomSymbols();
      if (!molB_Types) {
        std::cout << " Error in superimpose " << std::endl;
        return 1;
      }
    }

    double dRMSD = 0.0;
    double lRMSD = 0.0;

    if (type == 2) {
      std::vector<atom*> atomListA = pMoleculeA->getAtomList();
      std::vector<atom*> atomListB = pMoleculeB->getAtomList();

      for (int t = 0; t < nAtomsA; t++) {
        std::string nameA = atomListA[t]->getName();
        for (int t2 = 0; t2 < nAtomsB; t2++) {
          if (atomListB[t2]->getName() == nameA) {
            dRMSD += ( (coordsB[t2][0] - coordsA[t][0]) *
                        (coordsB[t2][0] - coordsA[t][0]) +
                        (coordsB[t2][1] - coordsA[t][1]) *
                        (coordsB[t2][1] - coordsA[t][1]) +
                        (coordsB[t2][2] - coordsA[t][2]) *
                        (coordsB[t2][2] - coordsA[t][2]) );
            nFittedAtoms++;
            break;
          }
        }
      }
    }
    else {
      bool gotOne = false;
      int tIndex = 0;
      for (int t = 0; t < nAtomsA; t++) {
        double d = BIGNUM;
        gotOne = false;
        int t2Index = 0;
        for (int t2 = 0; t2 < nAtomsB; t2++) {
          if (molA_Types[tIndex  ] == molB_Types[t2Index  ] and
              molA_Types[tIndex+1] == molB_Types[t2Index+1] ) {
            lRMSD = ( (coordsB[t2][0] - coordsA[t][0]) *
                      (coordsB[t2][0] - coordsA[t][0]) +
                      (coordsB[t2][1] - coordsA[t][1]) *
                      (coordsB[t2][1] - coordsA[t][1]) +
                      (coordsB[t2][2] - coordsA[t][2]) *
                      (coordsB[t2][2] - coordsA[t][2]) );
            gotOne = true;
            d = std::min(lRMSD, d);
          }
          t2Index+=2;
        }
        tIndex+=2;
        if (gotOne) {
          dRMSD += d;
          nFittedAtoms++;
        }
      }
    }

    if (dRMSD > 0.0) {
      dRMSD = sqrt( dRMSD / nFittedAtoms);
    }
    return dRMSD;
}

// ============================================================
// Function : initializeCorrespondences()
// ------------------------------------------------------------
//
// ============================================================
int superimpose::initializeCorrespondences(molecule* pMolA, int type,
                 std::vector<std::vector<int> > &subGraphs)
{
    //std::cout << "superimpose::initializeCorrespondences" << std::endl;
    this->correspondenceType = type;

    pMoleculeA = pMolA;
    int molAtoms = 0;

    char* molAtomTypes = NULL;
    int* molAdjMatrix = 0;

    if (type == 0 or type == 1) {
      molAdjMatrix = pMolA->getAdjMatrix();
      molAtoms = pMolA->getNumAtoms();
      if (type == 0) {
        molAtomTypes = pMoleculeA->getAtomTypes();
        if (!molAtomTypes) {
          std::cout << " Error in superimpose " << std::endl;
          return 1;
        }
      }
      else if (type == 1) {
        molAtomTypes = pMolA->getAtomSymbols();
        if (!molAtomTypes) {
          std::cout << " Error in superimpose " << std::endl;
          return 1;
        }
      }
    }
    else if (type == 2 or type == 3) {
      pMolA->generateHeavyAdjMatrix();
      molAdjMatrix = pMolA->getHeavyAdjMatrix();

      molAtoms = pMolA->getNumHeavyAtoms();
      this->nHeavyAtomsA = pMolA->getNumHeavyAtoms();
      this->heavyAtomIndicesA = pMolA->getHeavyAtomIndices();
      this->nHeavyAtomsB = this->nHeavyAtomsA;
      this->heavyAtomIndicesB = this->heavyAtomIndicesA;

      if (type == 2) {
        molAtomTypes = pMolA->getHeavyAtomTypes();
        if (!molAtomTypes) {
          std::cout << " Error in superimpose " << std::endl;
          return 1;
        }
      }
      else if (type == 3) {
        molAtomTypes = pMolA->getHeavyAtomSymbols();
        if (!molAtomTypes) {
          std::cout << " Error in superimpose " << std::endl;
          return 1;
        }
      }
    }

    // Keep all bond as single since COO- moeities have a single and double bond in MTK++
    for (int y = 0; y < molAtoms; y++) {
      for (int x = 0; x < molAtoms; x++) {
        if (molAdjMatrix[y*molAtoms+x]) {
          molAdjMatrix[y*molAtoms+x] = 1;
        }
      }
    }

    try {
      genMatchMatrix   = new int [molAtoms*molAtoms];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    for (int i = 0; i < molAtoms; i++) {
      for (int j = 0; j < molAtoms; j++) {

        //int currAtomKind = molAtomKinds[j];
        //if (fragAtomKinds[i] == 7) {
        //  currAtomKind = 7;
        //}

        if ((molAtomTypes[j*2  ] == molAtomTypes[i*2  ]) &&
            (molAtomTypes[j*2+1] == molAtomTypes[i*2+1])) {

          std::vector<int> matched;
          unsigned int nStdBonds = 0;
          for (int k = 0; k < molAtoms; k++) {
            if (i == k) continue;
            int pos = i*molAtoms+k;
            if (molAdjMatrix[pos]) {
              nStdBonds++;
              for (int l = 0; l < molAtoms; l++) {
                if (j == l) continue;
                //int pos2 = j*molAtoms+l;
                if (((molAtomTypes[l*2  ] == molAtomTypes[k*2  ]) &&
                     (molAtomTypes[l*2+1] == molAtomTypes[k*2+1])) //&&
                     //(
                     //  (molAtomSymbols[pos] == molAdjMatrix[pos2]) or
                     //  ((molAtomSymbols[pos] == 4) and (molAdjMatrix[pos2] == 6)) or
                     //  ((molAtomSymbols[pos] == 4) and (molAdjMatrix[pos2] == 7)) or
                     //  ((molAtomSymbols[pos] == 5) and (molAdjMatrix[pos2] == 1)) or
                     //  ((molAtomSymbols[pos] == 5) and (molAdjMatrix[pos2] == 2))
                     //)
                   ) {
                  std::vector<int>::iterator result;
                  result = std::find(matched.begin(), matched.end(), l);
                  if (result == matched.end()) {
                    matched.push_back(l);
                  }
                }
              }
            }
          }
          if (matched.size() >= nStdBonds) {
            genMatchMatrix[i*molAtoms+j] = 1;
          }
          else {
            genMatchMatrix[i*molAtoms+j] = 0;
          }
        }
        else {
          genMatchMatrix[i*molAtoms+j] = 0;
        }
      }
    }

/*
#ifdef DEBUG
          std::cout << "   superimpose::initializeCorrespondences\n    genMatchMatrix = " << std::endl;
          for (int y = 0; y < molAtoms; y++) {
            for (int x = 0; x < molAtoms; x++) {
              std::cout << genMatchMatrix[y*molAtoms+x] << " ";
            }
            std::cout << " " << std::endl;
          }
          std::cout << " " << std::endl;
#endif
*/

    //std::cout << " ullmann" << std::endl;
    std::vector<int> subGraph;
    functionalize* pFunc = new functionalize();
    pFunc->ullmann(0, molAtoms, molAdjMatrix, molAtomTypes,
                      molAtoms, molAdjMatrix, molAtomTypes,
                      genMatchMatrix, subGraph, subGraphs);
    delete pFunc;
    //std::cout << " ullmann done" << std::endl;

#ifdef DEBUG
    for (unsigned int k = 0; k < subGraphs.size(); k++) {
      std::cout << "   superimpose::initializeCorrespondences\n    subgraph : " << k+1 << std::endl;
      subGraph = subGraphs[k];
      for (int t = 0; t < molAtoms; t++) {
        for (int t2 = 0; t2 < molAtoms; t2++) {
          std::cout << subGraph[t*molAtoms+t2] << " ";
        }
        std::cout << " " << std::endl;
      }
    }
    std::cout << " " << std::endl;
#endif

    return 0;
}

// ============================================================
// Function : initializeCorrespondences()
// ------------------------------------------------------------
//
// ============================================================
int superimpose::initializeCorrespondences(molecule* pMolA, molecule* pMolB,
                 int type, std::vector<std::vector<int> > &subGraphs)
{
    if (type < 2 and type > 3) {
      return 1;
    }

    this->correspondenceType = type;

    pMoleculeA = pMolA;
    pMoleculeB = pMolB;
    std::vector<atom*> atomsA = pMoleculeA->getHeavyAtomList();
    std::vector<atom*> atomsB = pMoleculeB->getHeavyAtomList();
/*
    std::cout << " superimpose::initializeCorrespondences"
              << "\n  nAtomsA: " << atomsA.size()
              << "\n  nAtomsB: " << atomsB.size()
              << "\n  nBondsA: " << pMoleculeA->numBonds()
              << "\n  nBondsB: " << pMoleculeB->numBonds() << "\n";
*/
    int molAtomsA = 0;
    int molAtomsB = 0;

    char* molAtomTypesA = NULL;
    char* molAtomTypesB = NULL;
    int* molAdjMatrixA = 0;
    int* molAdjMatrixB = 0;

    pMoleculeA->generateHeavyAdjMatrix();
    pMoleculeB->generateHeavyAdjMatrix();

    molAdjMatrixA = pMoleculeA->getHeavyAdjMatrix();
    molAtomsA = pMoleculeA->getNumHeavyAtoms();
    this->nHeavyAtomsA = pMoleculeA->getNumHeavyAtoms();
    this->heavyAtomIndicesA = pMoleculeA->getHeavyAtomIndices();

    molAdjMatrixB = pMoleculeB->getHeavyAdjMatrix();
    molAtomsB = pMoleculeB->getNumHeavyAtoms();
    this->nHeavyAtomsB = pMoleculeB->getNumHeavyAtoms();
    this->heavyAtomIndicesB = pMoleculeB->getHeavyAtomIndices();

    std::string errMessage = "";
    if (nHeavyAtomsA != nHeavyAtomsB) {
      errMessage = "Number of heavy atoms don't match";
      MTKpp::errorLogger.throwError("superimpose::initializeCorrespondences", errMessage, MTK_ERROR);
      //return 1;
      throw MTKException(errMessage);
    }

    if (this->correspondenceType == 2) {
      molAtomTypesA = pMolA->getHeavyAtomTypes();
      molAtomTypesB = pMolB->getHeavyAtomTypes();

      if (!molAtomTypesA or !molAtomTypesB) {
        errMessage = " Can't determine heavy atom types";
        MTKpp::errorLogger.throwError("superimpose::initializeCorrespondences", errMessage, MTK_ERROR);
        //return 1;
        throw MTKException(errMessage);
      }
    }
    else if (this->correspondenceType == 3) {
      molAtomTypesA = pMolA->getHeavyAtomSymbols();
      molAtomTypesB = pMolB->getHeavyAtomSymbols();

      if (!molAtomTypesA or !molAtomTypesB) {
        errMessage = " Can't determine heavy atom symbols";
        MTKpp::errorLogger.throwError("superimpose::initializeCorrespondences", errMessage, MTK_ERROR);
        return 1;
      }
    }

    // Keep all bond as single since COO- moeities have a single and double bond in MTK++
    for (int y = 0; y < molAtomsA; y++) {
      for (int x = 0; x < molAtomsA; x++) {
        if (molAdjMatrixA[y*molAtomsA+x]) {
          molAdjMatrixA[y*molAtomsA+x] = 1;
        }
      }
    }

    for (int y = 0; y < molAtomsB; y++) {
      for (int x = 0; x < molAtomsB; x++) {
        if (molAdjMatrixB[y*molAtomsB+x]) {
          molAdjMatrixB[y*molAtomsB+x] = 1;
        }
      }
    }

    try {
      genMatchMatrix   = new int [molAtomsA*molAtomsB];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    for (int i = 0; i < molAtomsA; i++) {
      for (int j = 0; j < molAtomsB; j++) {

        if ((molAtomTypesB[j*2  ] == molAtomTypesA[i*2  ]) &&
            (molAtomTypesB[j*2+1] == molAtomTypesA[i*2+1])) {

          std::vector<int> matched;
          unsigned int nStdBonds = 0;
          for (int k = 0; k < molAtomsA; k++) {
            if (i == k) continue;
            int pos = i*molAtomsA+k;
            if (molAdjMatrixA[pos]) {
              nStdBonds++;
              for (int l = 0; l < molAtomsB; l++) {
                if (j == l) continue;
                int pos2 = j*molAtomsB+l;
                if (((molAtomTypesB[l*2  ] == molAtomTypesA[k*2  ]) &&
                     (molAtomTypesB[l*2+1] == molAtomTypesA[k*2+1])) &&
                     (
                       (molAdjMatrixA[pos] == molAdjMatrixB[pos2])
                     )
                   ) {
                  std::vector<int>::iterator result;
                  result = std::find(matched.begin(), matched.end(), l);
                  if (result == matched.end()) {
                    matched.push_back(l);
                  }
                }
              }
            }
          }
          if (matched.size() >= nStdBonds) {
            genMatchMatrix[i*molAtomsB+j] = 1;
          }
          else {
            genMatchMatrix[i*molAtomsB+j] = 0;
          }
        }
        else {
          genMatchMatrix[i*molAtomsB+j] = 0;
        }
      }
    }

    std::vector<int> subGraph;
    functionalize* pFunc = new functionalize();
    pFunc->ullmann(0, molAtomsA, molAdjMatrixA, molAtomTypesA,
                      molAtomsB, molAdjMatrixB, molAtomTypesB,
                      genMatchMatrix, subGraph, subGraphs);

    if (subGraphs.size() < 1) {
      //return 1;
      throw MTKException("error");
    }

    // Print
    for (unsigned int k = 0; k < subGraphs.size(); k++) {
      subGraph = subGraphs[k];
      int r = 0;
      //std::cout << "subgraph " << k+1 << " ";
      for (int t = 0; t < molAtomsA; t++) {
        for (int t2 = 0; t2 < molAtomsA; t2++) {
          if (subGraph[t*molAtomsA+t2]) {
            //std::cout << t+1 << ":" << t2+1 << ";";
            r++;
          }
        }
      }
      //std::cout << " " << std::endl;
      if (molAtomsA != r) {
        errMessage = " SubGraph fault ";
        //MTKpp::errorLogger.throwError("superimpose::initializeCorrespondences", errMessage, MTK_ERROR);
        //return 1;

        std::stringstream ss;
        ss << " superimpose::initializeCorrespondences error " << errMessage << molAtomsA << " " << r << std::endl;
        std::cout << ss.str();
        throw MTKException(ss.str());
      }
    }

    delete pFunc;
    return 0;
}

}  // MTKpp namespace
