/*!
   \file pca.cpp
   \brief Principal Component Analysis
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
   $Revision: 1.6 $

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

#include "sheet.h"
#include "table.h"

#include "pca.h"
// #include "Utils/diagonalize.h"

#include "Utils/diagonalize_eigen.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Class : pca()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pca::pca() {}

// ============================================================
// Class : pca()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pca::pca(table<double>* X, sheet* output):BaseStats()
{
    itsX = X;
    nRows = X->getNumRows();
    nColumns = X->getNumColumns();
    outModel = output;
}

// ============================================================
// Function : ~pca()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
//pca::~pca() {}

// ============================================================
// Function : run
// ------------------------------------------------------------
// Principal Component Analysis
// ============================================================
int pca::run(int nKeep)
{
    errorLogger.throwError("pca::run", "Principal Component Analysis Start", INFO);
    std::string errMessage = "";

    // Setup
    if (!this->itsX) {
      errorLogger.throwError("pca::run", "Can't find X ", MTK_ERROR);
      return 1;
    }

    //ublas::matrix<double> X = this->itsX->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> X = this->itsX->getMatrix();

    table<double>* xCenterTable = this->outModel->addTable();
    xCenterTable->setName("X-Centers");
    xCenterTable->setSizes(1,nColumns);

    //ublas::matrix<double> &xCenters = xCenterTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &xCenters = xCenterTable->getMatrix();

    int failure = getColumnCenters(X, xCenters);
    if (failure) {
      errorLogger.throwError("pca::run", " Error in getColumnCenters ", MTK_ERROR);
      return 1;
    }

    //ublas::matrix<double> B(nRows,nColumns);
    Eigen::Matrix<double, Dynamic, Dynamic> B(nRows,nColumns);

    // todo: create function in baseStats where I pass the xCenters
    // Calculate the deviations from the mean
    centerColumns(X, B);

    //ublas::matrix<double, ublas::column_major> C(nColumns,nColumns);
    Eigen::Matrix<double, Dynamic, Dynamic> C(nColumns,nColumns);

    // Find the covariance matrix
    covarianceMatrix(B, C);

/////// Eigen
    SelfAdjointEigenSolver<MatrixXd> eigensolver(C);
    MatrixXd eigenVectorsC = eigensolver.eigenvectors();
    VectorXd eigenValuesC = eigensolver.eigenvalues();
    eigenValueSort(eigenVectorsC, eigenValuesC, 1);
    std::cout << eigenValuesC <<std::endl;
/////// Eigen


/////// Boost
/*
    // Set up work vector
    ublas::vector<double> eValuesW;
    eValuesW.resize(nColumns);
    for (unsigned int i = 0; i < nColumns; i++) {
      eValuesW(i) = 0.0;
    }

    // Find the eigenvectors and eigenvalues of the covariance matrix
    int r = diagonalize(C, eValuesW);
    if (r != 0) {
      std::cout << " pca, diagonalize failed " << std::endl;
      return r;
    }

    // Sort in order of decreasing eigenvalues
    eigenValueSort(C, eValuesW, 1);
*/

    //ublas::matrix<double> newXAll(nRows,nColumns);
    //newXAll = prod(B,C);

    Eigen::Matrix<double, Dynamic, Dynamic> newXAll(nRows,nColumns);
    newXAll = B * eigenVectorsC;

    // Calculate the Cumulative energy
    //ublas::vector<double> cumEnergy;
    //cumEnergy.resize(nColumns);
    Eigen::Matrix<double, Dynamic, Dynamic> cumEnergy(nColumns, 1);

    double eValueTotal = 0.0;
    for (unsigned int j = 0; j < nColumns; j++) {
      //eValueTotal+=eValuesW[j];
      eValueTotal+=eigenValuesC(j);
    }

    for (unsigned int j = 0; j < nColumns; j++) {
      //cumEnergy(j) = eValuesW(j)/eValueTotal * 100;
      cumEnergy(j,0) = eigenValuesC(j)/eValueTotal * 100;
    }

//#ifdef DEBUG
    std::cout << " EigenValues (";
    for (unsigned int j = 0; j < nColumns; j++) {
      //std::cout << eValuesW(j) << " ";
      std::cout << eigenValuesC(j) << " ";
    }
    std::cout << ")\n" << std::endl;

    std::cout << " sdev (";
    for (unsigned int j = 0; j < nColumns; j++) {
      //std::cout << sqrt(eValuesW(j)) << " ";
      std::cout << sqrt(eigenValuesC(j)) << " ";
    }
    std::cout << ")\n" << std::endl;

    double t = 0.0;
    std::cout << " cumulative (";
    for (unsigned int j = 0; j < nColumns; j++) {
      std::cout << cumEnergy(j) << " ";
      t+=cumEnergy(j);
    }
    std::cout << ") cumTotal = " << t << " \n" << std::endl;

    for (int i = 0; i < 5; i++) {
      std::cout << " EigenVectors [" << i+1 << "]: (";
      for (unsigned int j = 0; j < nColumns; j++) {
        //std::cout << C(j,i) << " ";
        std::cout << eigenVectorsC(j,i) << " ";
      }
      std::cout << ")" << std::endl;
    }

    for (int i = 0; i < 5; i++) {
      std::cout << " newX [" << i+1 << "]: (";
      for (unsigned int j = 0; j < nRows; j++) {
        std::cout << newXAll(j,i) << " ";
      }
      std::cout << ")" << std::endl;
    }
//#endif

    // Determine number of components to keep
    int nKeep2 = 0;
    for (unsigned int j = 0; j < cumEnergy.size(); j++) {
      if (cumEnergy(j) > 5.0) {
        nKeep2++;
      }
    }
    int nComponents = std::max(nKeep, nKeep2);

    // Save nComponents eigenvectors, eigenvalues, sdev, cumulative
    table<double>* eVectorsTable = this->outModel->addTable();
    eVectorsTable->setName("eigenvectors");
    eVectorsTable->setSizes(nColumns, nComponents);
    //ublas::matrix<double> &eVectors = eVectorsTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &eVectors = eVectorsTable->getMatrix();

    for (int i = 0; i < nComponents; i++) {
      for (unsigned int j = 0; j < nColumns; j++) {
        //eVectors(j,i) = C(j,i);
        eVectors(j,i) = eigenVectorsC(j,i);
      }
    }

    table<double>* newXTable = this->outModel->addTable();
    newXTable->setName("newX");
    newXTable->setSizes(nRows, nComponents);
    //ublas::matrix<double> &newX = newXTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &newX = newXTable->getMatrix();

    for (int i = 0; i < nComponents; i++) {
      for (unsigned int j = 0; j < nRows; j++) {
        newX(j,i) = newXAll(j,i);
      }
    }

    table<double>* eValuesTable = this->outModel->addTable();
    eValuesTable->setName("eigenvalues");
    eValuesTable->setSizes(1, nComponents);
    //ublas::matrix<double> &eValues = eValuesTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &eValues = eValuesTable->getMatrix();

    table<double>* sdevTable = this->outModel->addTable();
    sdevTable->setName("sdev");
    sdevTable->setSizes(1, nComponents);
    //ublas::matrix<double> &sdev = sdevTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &sdev = sdevTable->getMatrix();

    table<double>* cumulativeTable = this->outModel->addTable();
    cumulativeTable->setName("cumulative");
    cumulativeTable->setSizes(1, nComponents);
    //ublas::matrix<double> &cumulative = cumulativeTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &cumulative = cumulativeTable->getMatrix();

    for (int j = 0; j < nComponents; j++) {
      //eValues(0,j) = eValuesW(j);
      //sdev(0,j) = sqrt(eValuesW(j));
      eValues(0,j) = eigenValuesC(j);
      sdev(0,j) = sqrt(eigenValuesC(j));
      cumulative(0,j) = cumEnergy(j);
    }

    // Set row labels
    for (unsigned int j = 0; j < nColumns; j++) {
      eVectorsTable->setRowLabel(j, this->itsX->getColumnLabel(j));
    }

    for (unsigned int j = 0; j < nRows; j++) {
      newXTable->setRowLabel(j, this->itsX->getRowLabel(j));
    }

    // Set column labels
    for (int i = 0; i < nComponents; i++) {
      std::stringstream number;
      number << i+1;
      std::string colLabel = "PC" + number.str();
      eValuesTable->setColumnLabel(i, colLabel);
      sdevTable->setColumnLabel(i, colLabel);
      cumulativeTable->setColumnLabel(i, colLabel);
      eVectorsTable->setColumnLabel(i, colLabel);
      newXTable->setColumnLabel(i, colLabel);
    }
    return 0;
}
    //"sdev"     "loadings" "center"   "scale"    "n.obs"    "scores"   "call"

} // MTKpp namespace

