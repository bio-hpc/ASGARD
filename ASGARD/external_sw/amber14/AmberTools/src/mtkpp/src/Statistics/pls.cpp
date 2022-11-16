/*!
   \file pls.cpp
   \brief Partial Least Squares
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
   $Revision: 1.5 $

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

#include "pls.h"
#include "Utils/diagonalize.h"
#include "Utils/inverseMatrix.h"
#include "Utils/randomNumbers.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>

namespace MTKpp
{

// ============================================================
// Class : pls()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pls::pls():BaseStats()
{
    itsX = 0;
    itsY = 0;
    itsMethod = "KERNELPLS";
    maxIter = 0;
    epsilon = 0.0;
    CV = "NONE";
    nITER = 0;
    nTEST = 0;
    SEED = 0;
    nLV = 1;
    outModel = 0;
}

// ============================================================
// Class : pls()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pls::pls(table<double>* Y, table<double>* X, std::string method, int nlvs,
         std::string cv, bool& bError):BaseStats()
{
    itsX = X;
    itsY = Y;
    YRows = Y->getNumRows();
    XColumns = X->getNumColumns();
    itsMethod = method;
    maxIter = 0;
    epsilon = 0.0;

    if (method != "KERNELPLS") {
      maxIter = 200;
      epsilon = 0.000001;
    }

    CV = cv;
    nLV = nlvs;
    nITER = 0;
    SEED = 0;
    nTEST = 0;
    outModel = 0;
}

// ============================================================
// Class : pls()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pls::pls(table<double>* Y, table<double>* X, std::string method, int nlvs, std::string cv,
         sheet* output, bool& bError):BaseStats()
{
    itsX = X;
    itsY = Y;
    YRows = Y->getNumRows();
    XColumns = X->getNumColumns();
    itsMethod = method;

    maxIter = 0;
    epsilon = 0.0;

    if (method != "KERNELPLS") {
      maxIter = 200;
      epsilon = 0.000001;
    }

    CV = cv;
    nLV = nlvs;
    nITER = 0;
    SEED = 0;
    nTEST = 0;
    outModel = output;
}

// ============================================================
// Class : pls()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pls::pls(std::string Y, std::string X, sheet* S, std::string method, int nlvs, std::string cv,
         sheet* output, bool& bError):BaseStats()
{
    if (S) {
      itsX = S->getTable(X);
      itsY = S->getTable(Y);
    }
    if (S and itsX and itsY) {
      YRows = itsY->getNumRows();
      XColumns = itsX->getNumColumns();

      itsMethod = method;

      maxIter = 0;
      epsilon = 0.0;

      if (method != "KERNELPLS") {
        maxIter = 200;
        epsilon = 0.000001;
      }

      CV = cv;
      nLV = nlvs;
      nITER = 0;
      SEED = 0;
      nTEST = 0;
      outModel = output;
    }
    else {
      bError = true;
    }
}

// ============================================================
// Function : ~pls()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
//pls::~pls() {}

// ============================================================
// Function : run()
// ------------------------------------------------------------
// Run
// ============================================================
void pls::run(bool& bError)
{
    if (outModel) {
      if (this->itsMethod == "KERNELPLS") {
        bError = this->kernelPLS();
      }
    }
    if (!bError) {
      if (CV != "NONE") {
        this->runCV(bError);
      }
    }
}

// ============================================================
// Function : runCV()
// ------------------------------------------------------------
// Run
// ============================================================
void pls::runCV(bool& bError)
{
    //ublas::matrix<int> samples;
    Eigen::Matrix<int, Dynamic, Dynamic> samples;

    if (CV == "LOO") {
      this->nTEST = 1;
      this->nEXT = this->YRows / this->nTEST;

      table<int>* looSamplesTable = this->outModel->addIntTable();
      looSamplesTable->setName("LOO Samples");
      looSamplesTable->setSizes(this->nEXT,this->nTEST);
      //ublas::matrix<int> &looSamples = looSamplesTable->getMatrix();
      Eigen::Matrix<int, Dynamic, Dynamic> &looSamples = looSamplesTable->getMatrix();
      for (unsigned int i = 0; i < this->nEXT; i++) {
        for (unsigned int j = 0; j < this->nTEST; j++) {
          looSamples(i,j) = static_cast<int>(i);
        }
      }
      samples = looSamplesTable->getMatrix();
      //looSamplesTable->print();
    }
    else if (CV == "LNO") {
      if (this->nTEST > 0 and this->nTEST < this->YRows) {
        if (this->YRows % this->nTEST == 0) {
          this->nEXT = this->YRows / this->nTEST;
          table<int>* lnoSamplesTable = this->outModel->addIntTable();
          lnoSamplesTable->setName("LNO Samples");
          lnoSamplesTable->setSizes(this->nEXT,this->nTEST);
          //ublas::matrix<int> &lnoSamples = lnoSamplesTable->getMatrix();
          Eigen::Matrix<int, Dynamic, Dynamic> &lnoSamples = lnoSamplesTable->getMatrix();
          int c = 0;
          for (unsigned int i = 0; i < this->nEXT; i++) {
            for (unsigned int j = 0; j < this->nTEST; j++) {
              lnoSamples(i,j) = c;
              c++;
            }
          }
          samples = lnoSamplesTable->getMatrix();
          //lnoSamplesTable->print();
        }
        else {
          std::cout << " ERROR " << std::endl;
        }
      }
    }
    else if (CV == "RANDOM") {
      if ( (this->nTEST > 0) and
           (this->nTEST < this->YRows) and
           (this->nITER > 0) ) {
        int nYs = static_cast<int>(this->YRows);
        table<int>* randomSamplesTable = this->outModel->addIntTable();
        randomSamplesTable->setName("RANDOM Samples");
        randomSamplesTable->setSizes(this->nITER, this->nTEST);
        //ublas::matrix<int> &randomSamples = randomSamplesTable->getMatrix();
        Eigen::Matrix<int, Dynamic, Dynamic> &randomSamples = randomSamplesTable->getMatrix();
        for (unsigned int i = 0; i < this->nITER; i++) {
          for (unsigned int j = 0; j < this->nTEST; j++) {
            bool gotit = true;
            int s = 0;
            while (gotit) {
              s = randomIntegerBetweenZeroAndX(nYs);
              gotit = false;
              for (unsigned int j2 = 0; j2 < this->nTEST; j2++) {
                if (randomSamples(i,j2) == s) {
                  gotit = true;
                }
              }
            }
            randomSamples(i,j) = s;
          }
        }
        samples = randomSamplesTable->getMatrix();
        //randomSamplesTable->print();
      }
    }
    else {
      std::cout << " UNKNOWN CV METHOS " << std::endl;
      bError = true;
    }

    if (!this->itsY or !this->itsX) {
      bError = true;
      return;
    }
/* // boost/eigen
    //ublas::matrix<double> Ys = this->itsY->getMatrix();
    //ublas::matrix<double> Xs = this->itsX->getMatrix();

    Eigen::Matrix<double, Dynamic, Dynamic> Ys = this->itsY->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Xs = this->itsX->getMatrix();

    sheet* subSheet = new sheet();
    sheet* cvSheet = new sheet();

    table<double>* ySubTable = subSheet->addTable();
    ySubTable->setName("Y");
    ySubTable->setSizes(Ys.size1() - samples.size2(), 1);
    //ublas::matrix<double> &ySub = ySubTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &ySub = ySubTable->getMatrix();

    table<double>* xSubTable = subSheet->addTable();
    xSubTable->setName("X");
    xSubTable->setSizes(Ys.size1() - samples.size2(), Xs.size2());
    //ublas::matrix<double> &xSub = xSubTable->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &xSub = xSubTable->getMatrix();

    table<double>* ySubPred = 0;

    unsigned int rowIndex = 0;
    bool gotit = false;

    for (unsigned int i = 0; i < samples.size1(); i++) {
      //std::cout << " RUNNING CV : " << i+1 << std::endl;
      rowIndex = 0;
      for (unsigned int j = 0; j < Ys.size1(); j++) {
        gotit = false;
        for (unsigned int j2 = 0; j2 < samples.size2(); j2++) {
          if (samples(i,j2) == static_cast<int>(j)) {
            gotit = true;
          }
        }
        if (!gotit) {
          for (unsigned int k = 0; k < Xs.size2(); k++) {
            xSub(rowIndex, k) = Xs(j,k);
          }
          for (unsigned int k = 0; k < Xs.size2(); k++) {
            ySub(rowIndex, 0) = Ys(j,0);
          }
          rowIndex++;
        }
      }
      ySubTable->print();
      xSubTable->print();

      pls* cvPLS = new pls(ySubTable, xSubTable, this->itsMethod, this->nLV, "NONE", cvSheet, bError);
      cvPLS->run(bError);
      if (bError) {
        std::cout << " ERROR IN CV " << std::endl;
        return;
      }
      //cvPLS->predict(ySubPredTable, xSubPredTable, cvSheet);

      // CALCULATE Q^2
      ySubPred = cvSheet->getTable("Y-Pred");
      if (ySubPred) {
        ySubPred->print();
      }
      delete cvPLS;
    }

    delete cvSheet;
    delete subSheet;
*/ // boost/eigen
}

// ============================================================
// Function : setX()
// ------------------------------------------------------------
// Set X matrix
// ============================================================
void pls::setX(table<double>* x)
{
    this->itsX = x;
}

// ============================================================
// Function : setY()
// ------------------------------------------------------------
// Set Y matrix
// ============================================================
void pls::setY(table<double>* y)
{
    this->itsY = y;
}

// ============================================================
// Function : setMethod()
// ------------------------------------------------------------
// Set PLS Algorithm
// ============================================================
void pls::setMethod(std::string s)
{
    this->itsMethod = s;
}

// ============================================================
// Function : setMaxIter()
// ------------------------------------------------------------
// Set Maximum number of iterations for iterative methods
// ============================================================
void pls::setMaxIter(int i)
{
    this->maxIter = i;
}

// ============================================================
// Function : setEpsilon()
// ------------------------------------------------------------
// Set Convergence criteria for iterative methods
// ============================================================
void pls::setEpsilon(double d)
{
    this->epsilon = d;
}

// ============================================================
// Function : setCV()
// ------------------------------------------------------------
// Set Cross validation method
// ============================================================
void pls::setCV(std::string s)
{
    this->CV = s;
}

// ============================================================
// Function : setNITER()
// ------------------------------------------------------------
// Set number of samples to consider in RANDOM CV
// ============================================================
void pls::setNITER(int i)
{
    this->nITER = i;
}

// ============================================================
// Function : setNTEST()
// ------------------------------------------------------------
// Set size of test set in RANDOM and LNO CV
// ============================================================
void pls::setNTEST(int i)
{
    this->nTEST = i;
}

// ============================================================
// Function : setSEED()
// ------------------------------------------------------------
// Set random number generator seed in RANDOM CV
// ============================================================
void pls::setSEED(int s)
{
    this->SEED = s;
    setRandomNumberSeed(s);
}

// ============================================================
// Function : setNLV()
// ------------------------------------------------------------
// Set number of Latent Variables to be considered
// ============================================================
void pls::setNLV(int l)
{
    this->nLV = l;
}

// ============================================================
// Function : setOutModel()
// ------------------------------------------------------------
// The sheet where the model is stored
// ============================================================
void pls::setOutModel(sheet* s)
{
    this->outModel = s;
}

// ****************************************************************************
//
//
//                             M  E  T  H  O  D  S
//
//
// ****************************************************************************

// ============================================================
// Function : kernelPLS(double)
// ------------------------------------------------------------
// Partial Least Squares Regression
// ============================================================
int pls::kernelPLS()
{
/* // boost/eigen
    errorLogger.throwError("pls::kernelPLS", "Start", INFO);
    std::string errMessage = "";

    // Setup
    if (!this->itsY or !this->itsX) return 1;
    //ublas::matrix<double> Ys = this->itsY->getMatrix();
    //ublas::matrix<double> Xs = this->itsX->getMatrix();

    Eigen::Matrix<double, Dynamic, Dynamic> Ys = this->itsY->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Xs = this->itsX->getMatrix();

    int N = Ys.size1();
    int M = Ys.size2();
    int R = Xs.size2();
    int MaxComponents = N;

    if (M != 1) {
      errorLogger.throwError("pls::kernelPLS", "Need to implement multiple dependent variable", MTK_ERROR);
      return 1;
    }

    std::stringstream N_ss;
    N_ss << N;
    std::string N_str = N_ss.str();

    std::stringstream M_ss;
    M_ss << M;
    std::string M_str = M_ss.str();

    std::stringstream R_ss;
    R_ss << R;
    std::string R_str = R_ss.str();

    std::stringstream MaxComponents_ss;
    MaxComponents_ss << MaxComponents;
    std::string MaxComponents_str = MaxComponents_ss.str();

    // Calculate the Maximum number of components
    //
    //if (R <= N) {
    //  MaxComponents = R;
    //  if (nLV < MaxComponents) {
    //    MaxComponents = nLV;
    //  }
    //}
    //

    MaxComponents = std::min(std::min(R, N), static_cast<int>(nLV));

    errMessage = " Dimensions: \n ### Y[" + N_str + "][" + M_str + "]\n";
    errMessage += " ### X[" + N_str + "][" + R_str +"]\n";
    errMessage += " ### B[" + R_str + "][" + M_str + "]\n";
    errMessage += " ### Number of LVs = " + MaxComponents_str;
    errorLogger.throwError("pls::kernelPLS", errMessage, INFO);

    // Clear output sheet
    this->outModel->clear();

    // Create/Initialize Work arrays
    //ublas::matrix<double, ublas::column_major> kernel(R,R);
    Eigen::Matrix<double, Dynamic, Dynamic> kernel(R,R);

    for (int i = 0; i < R; i++) {
      for (int j = 0; j < R; j++) {
        kernel(i, j) = 0.0;
      }
    }

    //ublas::matrix<double> F(N,M);
    Eigen::Matrix<double, Dynamic, Dynamic> F(N,M);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        F(i, j) = 0.0;
      }
    }

    //ublas::matrix<double> F_tr(M,N);
    Eigen::Matrix<double, Dynamic, Dynamic> F_tr(M,N);
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        F_tr(i, j) = 0.0;
      }
    }

    //ublas::matrix<double> E(N,R);
    Eigen::Matrix<double, Dynamic, Dynamic> E(N,R);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < R; j++) {
        E(i, j) = 0.0;
      }
    }

    //ublas::matrix<double> E_tr(R,N);
    Eigen::Matrix<double, Dynamic, Dynamic> E_tr(R,N);
    for (int i = 0; i < R; i++) {
      for (int j = 0; j < N; j++) {
        E_tr(i, j) = 0.0;
      }
    }

    //ublas::matrix<double> E_trXF(R,M);
    Eigen::Matrix<double, Dynamic, Dynamic> E_trXF(R,M);
    for (int i = 0; i < R; i++) {
      for (int j = 0; j < M; j++) {
        E_trXF(i, j) = 0.0;
      }
    }

    //ublas::matrix<double> F_trXE(M,R);
    Eigen::Matrix<double, Dynamic, Dynamic> F_trXE(M,R);
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < R; j++) {
        F_trXE(i, j) = 0.0;
      }
    }

    // Create/Initialize work vectors
    ublas::vector<double> t(N);
    ublas::vector<double> u(N);

    for (int i = 0; i < N; i++) {
      t(i) = 0.0;
      u(i) = 0.0;
    }

    ublas::vector<double> c(M);
    for (int i = 0; i < M; i++) {
      c(i) = 0.0;
    }

    ublas::vector<double> eigenValues(R);
    ublas::vector<double> p(R);
    ublas::vector<double> w(R);
    for (int i = 0; i < R; i++) {
      eigenValues(i) = 0.0;
      p(i) = 0.0;
      w(i) = 0.0;
    }

    // Create output tables
    table<double>* wTable = this->outModel->addTable();
    wTable->setName("X-Weights");
    wTable->setSizes(R, MaxComponents);
    wTable->initialize(0.0);
    ublas::matrix<double> &W = wTable->getMatrix();

    table<double>* pTable = this->outModel->addTable();
    pTable->setName("X-Loadings");
    pTable->setSizes(R, MaxComponents);
    pTable->initialize(0.0);
    ublas::matrix<double> &P = pTable->getMatrix();

    table<double>* cTable = this->outModel->addTable();
    cTable->setName("Y-Weights");
    cTable->setSizes(M, MaxComponents);
    cTable->initialize(0.0);
    ublas::matrix<double> &C = cTable->getMatrix();

    table<double>* xCenterTable = this->outModel->addTable();
    xCenterTable->setName("X-Centers");
    xCenterTable->setSizes(1, R);
    xCenterTable->initialize(0.0);
    ublas::matrix<double> &xCenters = xCenterTable->getMatrix();

    table<double>* yCenterTable = this->outModel->addTable();
    yCenterTable->setName("Y-Centers");
    yCenterTable->setSizes(1, M);
    yCenterTable->initialize(0.0);
    ublas::matrix<double> &yCenters = yCenterTable->getMatrix();

    table<double>* xStdDevTable = this->outModel->addTable();
    xStdDevTable->setName("X-StdDevs");
    xStdDevTable->setSizes(1, R);
    xStdDevTable->initialize(0.0);
    ublas::matrix<double> &xStdDevs = xStdDevTable->getMatrix();

    table<double>* yStdDevTable = this->outModel->addTable();
    yStdDevTable->setName("Y-StdDevs");
    yStdDevTable->setSizes(1, M);
    yStdDevTable->initialize(0.0);
    ublas::matrix<double> &yStdDevs = yStdDevTable->getMatrix();

    // Calculate Column means
    int failure = getColumnCenters(Ys, yCenters);
    if (failure) {
      errorLogger.throwError("pls::kernelPLS", "Error calculating Y column centers", MTK_ERROR);
      return 1;
    }

    failure = getColumnCenters(Xs, xCenters);
    if (failure) {
      errorLogger.throwError("pls::kernelPLS", "Error calculating X column centers", MTK_ERROR);
      return 1;
    }

    // Calculate Column standard deviations
    failure = getStdDevColumns(Ys, yStdDevs);
    if (failure) {
      errorLogger.throwError("pls::kernelPLS", "Error calculating Y stdandard deviations", MTK_ERROR);
      return 1;
    }

    failure = getStdDevColumns(Xs, xStdDevs);
    if (failure) {
      errorLogger.throwError("pls::kernelPLS", "Error calculating X stdandard deviations", MTK_ERROR);
      return 1;
    }

    // Autoscale ==> z-score = (c[i] - mean[-] )/ stdDev[i]
    failure = autoScale(Ys, F, yCenters, yStdDevs);
    if (failure) {
      errorLogger.throwError("pls::kernelPLS", "Error scaling Y", MTK_ERROR);
      return 1;
    }

    failure = autoScale(Xs, E, xCenters, xStdDevs);
    if (failure) {
      errorLogger.throwError("pls::kernelPLS", "Error scaling X", MTK_ERROR);
      return 1;
    }

    errorLogger.throwError("pls::kernelPLS", "Setup finished", INFO);

    // Start Kernel PLS Algorithm
    E_tr = ublas::trans(E);
    // - Loop over all Latent Variables
    for (int k = 0; k < MaxComponents; k++) {
      E_tr.clear();
      F_tr.clear();

      E_tr = ublas::trans(E);
      F_tr = ublas::trans(F);

      E_trXF = ublas::prod(E_tr,F);
      F_trXE = ublas::prod(F_tr,E);
      kernel = ublas::prod(E_trXF,F_trXE);

      std::stringstream k_str;
      k_str << k+1;

      //errMessage = "Diagonalize " + k_str.str();
      //errorLogger.throwError("pls::kernelPLS", errMessage, INFO);

      int result = diagonalize(kernel, eigenValues);
      if (result != 0) {
        errorLogger.throwError("pls::kernelPLS", "Error in diagonization", MTK_ERROR);
        return result;
      }

      eigenValueSort(kernel, eigenValues, 1);

      //std::cout << " EigenVectors [" << k+1 << "]: (";
      for (unsigned int j = 0; j < kernel.size2(); j++) {
        w(j) = kernel(j,0);
        //std::cout << w(j) << " ";
      }
      //std::cout << ")" << std::endl;

      t = ublas::prod(E,w);

      double tt = ublas::inner_prod(t,t);

      c = ublas::prod(F_tr,t)/tt;

      double cc = ublas::inner_prod(c,c);

      u = ublas::prod(F,c)/cc;

      p = ublas::prod(E_tr,t)/tt;

      for (unsigned int j = 0; j < W.size1(); j++) {
        W(j,k) = w(j);
      }

      for (unsigned int j = 0; j < P.size1(); j++) {
        P(j,k) = p(j);
      }

      for (unsigned int j = 0; j < C.size1(); j++) {
        C(j,k) = c(j);
      }

      for (unsigned int i = 0; i < E.size1(); i++) {
        for (unsigned int j = 0; j < E.size2(); j++) {
          E(i,j) = E(i,j) - t(i) * p(j);
        }
      }

      for (unsigned int i = 0; i < F.size1(); i++) {
        for (unsigned int j = 0; j < F.size2(); j++) {
          F(i,j) = F(i,j) - t(i) * C(j,k);
        }
      }
    }

    //
    // B = W * inv(P' * W) * C'
    //   = [R, MaxComponents] * inv([MaxComponents, R] * [R, MaxComponents]) * [MaxComponents, M]
    //   = [R, MaxComponents] * [MaxComponents, MaxComponents] * [MaxComponents, M]
    //   = [R, MaxComponents] * [MaxComponents, M]
    //   = [R, M]
    //

    ublas::matrix<double> B(R,M);

    table<double>* regCoeffTable = this->outModel->addTable();
    regCoeffTable->setName("Regression Coefficients");
    regCoeffTable->setSizes(MaxComponents, R);
    regCoeffTable->setColumnLabels(this->itsX->getColumnLabels());

    for (int m = 0; m < MaxComponents; m++) {
      std::stringstream label;
      label << m+1;
      std::string str_label = ("r" + label.str()).c_str();
      regCoeffTable->setRowLabel(m, str_label);
    }

    ublas::matrix<double> &regCoeff = regCoeffTable->getMatrix();

    for (int m = 0; m < MaxComponents; m++) {
      ublas::matrix<double> PtW(m+1,m+1);
      ublas::matrix<double> invPtW(m+1,m+1);
      ublas::matrix<double> WinvPtW(R, m+1);

      ublas::matrix<double> tempP(R, m+1);
      ublas::matrix<double> tempW(R, m+1);
      ublas::matrix<double> tempC(M, m+1);
      for (int m2 = 0; m2 < m+1; m2++) {
        for (int r = 0; r < R; r++) {
          tempP(r,m2) = P(r,m2);
          tempW(r,m2) = W(r,m2);
        }
        for (int m3 = 0; m3 < M; m3++) {
          tempC(m3,m2) = C(m3,m2);
        }
      }
      PtW = ublas::prod(ublas::trans(tempP), tempW);
      InvertMatrix(PtW, invPtW);
      WinvPtW = ublas::prod(tempW, invPtW);
      B = ublas::prod(WinvPtW, ublas::trans(tempC));
      for (int r = 0; r < R; r++) {
        regCoeff(m,r) = B(r,0);
      }
    }

    // Create storage for Y-Pred
    table<double>* yPredTable = this->outModel->addTable();
    yPredTable->setName("Y-Pred");
    yPredTable->setSizes(N, MaxComponents);
    ublas::matrix<double> &yPred = yPredTable->getMatrix();

    for (int m = 0; m < MaxComponents; m++) {
      std::stringstream label;
      label << m+1;
      std::string str_label = ("r" + label.str()).c_str();
      yPredTable->setColumnLabel(m, str_label);
    }
    yPredTable->setRowLabels(this->itsY->getRowLabels());

    ublas::matrix<double> yPredc(N,1);

    failure = autoScale(Xs, E, xCenters, xStdDevs);
    if (failure) return 1;

    for (int m = 0; m < MaxComponents; m++) {
      for (int r = 0; r < R; r++) {
        B(r,0) = regCoeff(m,r);
      }
      yPredc = ublas::prod(E, B);
      yPredc = ublas::prod(yPredc, yStdDevs);

      for (unsigned int i = 0; i < yPredc.size1(); i++) {
        yPredc(i,0) += yCenters(0,0);
      }

      for (int i = 0; i < N; i++) {
        yPred(i,m) = yPredc(i,0);
      }
    }

    ///////////////////////////////////
    // Create storage for R^2
    ///////////////////////////////////
    table<double>* r2Table = this->outModel->addTable();
    r2Table->setName("R2");
    r2Table->setSizes(1, MaxComponents);
    ublas::matrix<double> &r2 = r2Table->getMatrix();
    r2Table->setColumnLabels(yPredTable->getColumnLabels());
    r2Table->setRowLabel(0, "R2");

    for (int m = 0; m < MaxComponents; m++) {
      double dR2 = rSquaredColumn(this->itsY->getMatrix(), 0, yPred, m);
      r2(0,m) = dR2;
    }

    ///////////////////////////////////
    // Create storage for Average Unsigned Error
    ///////////////////////////////////
    table<double>* aueTable = this->outModel->addTable();
    aueTable->setName("Average Unsigned Error");
    aueTable->setSizes(1, MaxComponents);
    ublas::matrix<double> &aue = aueTable->getMatrix();
    aueTable->setColumnLabels(yPredTable->getColumnLabels());
    aueTable->setRowLabel(0, "Average Unsigned Error");

    for (int m = 0; m < MaxComponents; m++) {
      double dAUE = UnsignedError(this->itsY->getMatrix(), 0, yPred, m);
      aue(0,m) = dAUE;
    }

    ///////////////////////////////////
    // Create storage for Root Mean Squared Error
    ///////////////////////////////////
    table<double>* rmseTable = this->outModel->addTable();
    rmseTable->setName("Root Mean Squared Error");
    rmseTable->setSizes(1, MaxComponents);
    ublas::matrix<double> &rmse = rmseTable->getMatrix();
    rmseTable->setColumnLabels(yPredTable->getColumnLabels());
    rmseTable->setRowLabel(0, "Root Mean Squared Error");

    for (int m = 0; m < MaxComponents; m++) {
      double dRMSE = RMSE(this->itsY->getMatrix(), 0, yPred, m);
      rmse(0,m) = dRMSE;
    }
*/
    return 0;
}

} // MTKpp namespace

