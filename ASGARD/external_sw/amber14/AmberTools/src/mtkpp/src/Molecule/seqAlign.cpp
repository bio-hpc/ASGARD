/*!
   \file seqAlign.cpp
   \brief Sequence Alignment
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
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

#include "seqAlign.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "utility.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : seqAlign()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
seqAlign::seqAlign()
{
    // Template
    tMol = 0;
    tNRes = 0;
    tCSeq = 0;
    tISeq = 0;

    // Query
    qMol = 0;
    qNRes = 0;
    qCSeq = 0;
    qISeq = 0;

    // Files and Parameters
    pamFile = "";
    algorithmType = 1;
    gapPenaltyType = 1;
    gapOpen = 2.0;
    gapExtend = 0.0;
    maxGap = 10;
    gapSize = 0;

    // Match Matrix
    matchMatrix = 0;
    maxIndices = 0;
    bestI = 0;
    bestJ = 0;
    bestMatchMatrixValue = 0.0;

    corrMap = 0;
}

// ============================================================
// Function : ~seqAlign()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
seqAlign::~seqAlign()
{
    delete [] tCSeq;
    delete [] tISeq;
    delete [] qCSeq;
    delete [] qISeq;
    delete [] pamTypes;
    delete [] pam;
    delete [] matchMatrix;
    delete [] maxIndices;
    delete [] corrMap;
}

// ============================================================
// Function :run()
// ------------------------------------------------------------
//
// ============================================================
int seqAlign::run()
{
    // Determine Match Matrix
    int failure = this->initMatchMatrix();
    if (failure) {
      std::string errMessage = " Initializing Match Matrix Failure";
      MTKpp::errorLogger.throwError("seqAlign::run", errMessage, 1);
      return 1;
    }

    // Determine the alignment
    failure = this->traceBack();

    // Print results
    if (failure) {
      std::string errMessage = " Trace Back Failure";
      MTKpp::errorLogger.throwError("seqAlign::run", errMessage, 1);
      return 1;
    }
    else {
      this->printResults();
    }
    return failure;
}

// ============================================================
// Function :setTemplate()
// ------------------------------------------------------------
//
// ============================================================
int seqAlign::setTemplate(molecule* pMolecule)
{
    this->tMol = pMolecule;
    tNRes = pMolecule->getNumSubMolecules();

    char* tCSeqTemp;

    try {
      tCSeqTemp = new char [this->tNRes];

      this->tCSeq     = new char [this->tNRes+1]; // character sequence
      this->tISeq     = new int  [this->tNRes+1]; // int sequence
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    pMolecule->getRes1LSymbols(tCSeqTemp);

    std::string errorMessage =  " Template Sequence \n";

    for (int i = 0; i < this->tNRes; i++) {
      //std::cout << tCSeqTemp[i];
      errorMessage += tCSeqTemp[i];
    }
    //std::cout << " " << std::endl;
    MTKpp::errorLogger.throwError(" seqAlign::setTemplate ", errorMessage, INFO);

    this->tCSeq[0] = 'd';
    this->tISeq[0] = -1;
    for (int i = 0; i < this->tNRes; i++) {
      this->tCSeq[i+1] = tCSeqTemp[i];
      for (int j = 0; j < this->pamSize; j++) {
        if (this->pamTypes[j] == tCSeqTemp[i]) {
          this->tISeq[i+1] = j;
          break;
        }
      }
    }

    delete [] tCSeqTemp;

    return 0;
}

// ============================================================
// Function :setQuery()
// ------------------------------------------------------------
//
// ============================================================
int seqAlign::setQuery(molecule* pMolecule)
{
    this->qMol = pMolecule;
    qNRes = pMolecule->getNumSubMolecules();

    char* qCSeqTemp;

    try {
      qCSeqTemp = new char [this->qNRes];

      this->qCSeq     = new char [this->qNRes+1]; // character sequence
      this->qISeq     = new int  [this->qNRes+1]; // int sequence
      this->corrMap   = new int  [this->qNRes  ]; // correspondence map
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    pMolecule->getRes1LSymbols(qCSeqTemp);

    std::string errorMessage =  " Query Sequence \n";
    for (int i = 0; i < this->qNRes; i++) {
        errorMessage += qCSeqTemp[i];
    }
    MTKpp::errorLogger.throwError(" seqAlign::setQuery ", errorMessage, INFO);
/*
    std::cout << " seqAlign::setQuery " << qNRes << std::endl;
    for (int i = 0; i < this->qNRes; i++) {
        std::cout << qCSeqTemp[i];
    }
    std::cout << " " << std::endl;
*/
    this->qCSeq[0] = 'd';
    this->qISeq[0] = -1;
    for (int i = 0; i < this->qNRes; i++) {
      this->qCSeq[i+1] = qCSeqTemp[i];
      this->corrMap[i] = 0;
      for (int j = 0; j < this->pamSize; j++) {
        if (this->pamTypes[j] == qCSeqTemp[i]) {
          this->qISeq[i+1] = j;
          break;
        }
      }
    }

    delete [] qCSeqTemp;

    return 0;
}

// ============================================================
// Function :initMatchMatrix()
// ------------------------------------------------------------
//
// ============================================================
int seqAlign::initMatchMatrix()
{
    std::string errorMessage;
    if (this->pamFile == "") {
      return 1;
    }
    errorMessage = "\n Number of residues in Target = " + i2s(this->tNRes) + "\n";
    errorMessage += " Number of residues in Query = " + i2s(this->qNRes) + "\n";

    int B2MB = 1024 * 1024;
    int doubleSize = sizeof(double);
    int intSize = sizeof(int);

    int mmSizeBytes = doubleSize * ((this->tNRes+1) * (this->qNRes+1));
    int mm2SizeBytes = intSize * ((this->tNRes+1) * (this->qNRes+1));

    int mmSize = (doubleSize * ((this->tNRes+1) * (this->qNRes+1))) / B2MB;
    int mm2Size = (intSize * ((this->tNRes+1) * (this->qNRes+1))) / B2MB;

    errorMessage +=  "\n Allocating " + i2s(mmSizeBytes + mm2SizeBytes) + " Bs or " + i2s(mmSize + mm2Size) + " MBs \n";

    try {
      this->matchMatrix = new double [(this->tNRes+1) * (this->qNRes+1)];
      this->maxIndices = new int [(this->tNRes+1) * (this->qNRes+1)];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    for (int i = 0; i < tNRes+1; i++) {
      for (int j = 0; j < qNRes+1; j++) {
        this->matchMatrix[(i*(qNRes+1))+j] = 0.0;
        this->maxIndices[(i*(qNRes+1))+j] = 0;
      }
    }

    double H_N = 0.0;
    double H_NW = 0.0;
    double H_W = 0.0;

    int H_N_index = 0;
    int H_NW_index = 0;
    int H_W_index = 0;
    int H_index = 0;

    double pamScore = 0.0;
    double bestScore = 0.0;
    int bestIndex = 0;

    // Set first column
    for (int i = 0; i < tNRes+1; i++) {
      H_index = i * (qNRes+1) + 0;
      H_N_index   = (i-1) * (qNRes+1) + 0;
      if (i > 0) {
        this->matchMatrix[H_index] = this->matchMatrix[H_N_index] - this->gapOpen;
        this->matchMatrix[H_index] = 0;
      }
      else {
        this->matchMatrix[H_index] = 0;
      }
    }

    // Set first row
    for (int j = 0; j < qNRes+1; j++) {
      H_index = j;
      H_W_index = j - 1;
      if (j > 0) {
        this->matchMatrix[H_index] = this->matchMatrix[H_W_index] - this->gapOpen;
        this->matchMatrix[H_index] = 0;
      }
      else {
        this->matchMatrix[H_index] = 0;
      }
    }

    if (this->algorithmType == 1) {
      errorMessage += " Using the Needleman-Wunsch Algorithm \n";
    }
    else if (this->algorithmType == 1) {
      errorMessage += " Using the Smith-Waterman Algorithm \n";
    }

    if (this->algorithmType == 1 or this->algorithmType == 2) { // Needleman-Wunsch
      for (int i = 1; i < tNRes+1; i++) {

        for (int j = 1; j < qNRes+1; j++) {
          // Calculate i-j PAM Score
          pamScore = this->pam[tISeq[i] * this->pamSize + qISeq[j]];

          H_index     = i     * (qNRes+1) + j;
          H_NW_index  = (i-1) * (qNRes+1) + (j-1);
          H_NW = this->matchMatrix[H_NW_index] + pamScore;
          int ki = std::max(1, i - maxGap);
          int kj = std::max(1, j - maxGap);

          bestScore = H_NW;
          bestIndex = H_NW_index;

          // Smith-Waterman
          if (this->algorithmType == 2) {
            if (bestScore < 0.0) {
              bestScore = 0.0;
              bestIndex = -1;
            }
          }

          // Constant gap penalty (gapPenalty = d)
          if (this->gapPenaltyType == 1) {
            for (int k = i-1; k > ki; k--) {
              H_N_index   = k * (qNRes+1) + j;
              H_N  = this->matchMatrix[H_N_index ] - this->gapOpen;
              if (H_N > bestScore) {
                bestScore = H_N;
              }
              if (this->matchMatrix[H_N_index ] > bestScore) {
                bestIndex = H_N_index;
              }
            }

            for (int k = j-1; k > kj; k--) {
              H_W_index   = i * (qNRes+1) + k;
              H_W  = this->matchMatrix[H_W_index ] - this->gapOpen;
              if (H_W > bestScore) {
                bestScore = H_W;
              }
              if (this->matchMatrix[H_W_index ] > bestScore) {
                bestIndex = H_W_index;
              }
            }
          }
          // Linear gap penalty   (gapPenalty = -g*d)
          else if (this->gapPenaltyType == 2) {
            for (int k = i-1; k > ki; k--) {
              gapSize = i - k;
              H_N_index   = k * (qNRes+1) + j;
              H_N  = this->matchMatrix[H_N_index ] - (this->gapOpen * gapSize);
              if (H_N > bestScore) {
                bestScore = H_N;
              }
              if (this->matchMatrix[H_N_index ] > bestScore) {
                bestIndex = H_N_index;
              }
            }

            for (int k = j-1; k > kj; k--) {
              gapSize = j - k;
              H_W_index   = i * (qNRes+1) + k;
              H_W  = this->matchMatrix[H_W_index ] - (this->gapOpen * gapSize);
              if (H_W > bestScore) {
                bestScore = H_W;
              }
              if (this->matchMatrix[H_W_index ] > bestScore) {
                bestIndex = H_W_index;
              }
            }
          }
          // Affine gap penalty   (gapPenalty = -d - (g-1) * e, where e < d)
          else if (this->gapPenaltyType == 3) {
            for (int k = i-1; k > ki; k--) {
              gapSize = i - k;
              H_N_index   = k * (qNRes+1) + j;
              H_N  = this->matchMatrix[H_N_index ] - (this->gapOpen - (gapSize-1) * this->gapExtend);
              if (H_N > bestScore) {
                bestScore = H_N;
                //bestIndex = H_N_index;
              }
              if (this->matchMatrix[H_N_index ] > bestScore) {
                  bestIndex = H_N_index;
              }
            }

            for (int k = j-1; k > kj; k--) {
              gapSize = j - k;
              H_W_index   = i * (qNRes+1) + k;
              H_W  = this->matchMatrix[H_W_index ] - (this->gapOpen - (gapSize-1) * this->gapExtend);
              if (H_W > bestScore) {
                bestScore = H_W;
              //  bestIndex = H_W_index;
              }
              if (this->matchMatrix[H_W_index ] > bestScore) {
                  bestIndex = H_W_index;
              }
            }
          }

          this->matchMatrix[H_index] = bestScore;
          this->maxIndices[H_index]  = bestIndex;

          if (this->algorithmType == 2 and (bestScore > bestMatchMatrixValue)) {
            bestMatchMatrixValue = bestScore;
            bestI = i;
            bestJ = j;
          }

        }
      }
      if (this->algorithmType == 1) {
        H_index = (tNRes+1) * (qNRes+1) - 1;
        bestMatchMatrixValue = this->matchMatrix[H_index];
        bestI = tNRes;
        bestJ = qNRes;
      }

      //errorMessage +=  " Success ";
      MTKpp::errorLogger.throwError(" seqAlign::initMatchMatrix ", errorMessage, 4);
    }
    else {
      std::string errorMessage =  " Unknown algorithm type: " + i2s(this->algorithmType) + ", options are 1 or 2.";
      MTKpp::errorLogger.throwError(" seqAlign::initialize ", errorMessage, 1);
      return 1;
    }
    return 0;
}

// ============================================================
// Function :traceBack()
// ------------------------------------------------------------
//
// ============================================================
int seqAlign::traceBack()
{
    std::string errorMessage =  " start ";
    MTKpp::errorLogger.throwError(" seqAlign::traceBack ", errorMessage, 4);
/*
    // Template
    std::cout << " TEMPLATE ";
    for (int i = 0; i < tNRes+1; i++) {
      std::cout << this->tCSeq[i];
    }
    std::cout << "\n" << std::endl;

    // Query
    std::cout << " QUERY    ";

    for (int j = 0; j < qNRes+1; j++) {
      std::cout << this->qCSeq[j];
    }
    std::cout << "\n" << std::endl;

    // Print  Indices
    std::cout << "INDICES\n";
    std::cout << std::setw(6) << "-";
    for (int j = 0; j < qNRes+1; j++) {
      std::cout << std::setw(7) << this->qCSeq[j];
    }
    std::cout << " " << std::endl;

    for (int i = 0; i < tNRes+1; i++) {
      std::cout << std::setw(6) << this->tCSeq[i] << " ";
      for (int j = 0; j < qNRes+1; j++) {
        std::cout << std::showpoint << std::setprecision(4) << std::setw(6)
                  << (i*(qNRes+1))+j << " ";
      }
      std::cout << " " << std::endl;
    }

    // Print Match Matrix
    std::cout << "MATCH MATRIX\n";
    std::cout << std::setw(6) << "-";
    for (int j = 0; j < qNRes+1; j++) {
      std::cout << std::setw(7) << this->qCSeq[j];
    }
    std::cout << " " << std::endl;

    for (int i = 0; i < tNRes+1; i++) {
      std::cout << std::setw(6) << this->tCSeq[i] << " ";
      for (int j = 0; j < qNRes+1; j++) {
        std::cout << std::showpoint << std::setprecision(4) << std::setw(6)
                  << this->matchMatrix[(i*(qNRes+1))+j] << " ";
      }
      std::cout << " " << std::endl;
    }

    // Print Max Indices
    std::cout << "MAX INDICES\n";
    std::cout << std::setw(6) << "-";
    for (int j = 0; j < qNRes+1; j++) {
      std::cout << std::setw(7) << this->qCSeq[j];
    }
    std::cout << " " << std::endl;

    for (int i = 0; i < tNRes+1; i++) {
      std::cout << std::setw(6) << this->tCSeq[i] << " ";
      for (int j = 0; j < qNRes+1; j++) {
        std::cout << std::showpoint << std::setprecision(4) << std::setw(6)
                  << this->maxIndices[(i*(qNRes+1))+j] << " ";
      }
      std::cout << " " << std::endl;
    }

    //for (unsigned int t = 0; t < this->qNRes; t++) {
    //  std::cout << this->corrMap[t] << " ";
    //}
    //std::cout << " " << std::endl;
*/
    alignmentA = "";
    alignmentB = "";

    int H_N_index = 0;
    int H_NW_index = 0;
    int H_W_index = 0;
    int H_index = 0;

    if (this->algorithmType == 1) { // Needleman-Wunsch
      // Start at H(m,n) where m is the number of rows and n is the number of columns
      int i = bestI;
      int j = bestJ;
      while (i > 0 and j > 0) {
        H_index     = i     * (qNRes+1) + j;
        H_NW_index  = (i-1) * (qNRes+1) + (j-1);

        if (maxIndices[H_index] == H_NW_index) {
          alignmentA = tCSeq[i] + alignmentA;
          alignmentB = qCSeq[j] + alignmentB;
          this->corrMap[j-1] = i;
          i--;
          j--;
          //std::cout << " i: " << i <<  " j: " << j << std::endl;
        }
        else {
          int ki = std::max(1, i - maxGap);
          int kj = std::max(1, j - maxGap);
          //std::cout << " ki: " << ki << " kj: " << kj << std::endl;
          int c = 0;
          bool found = false;
          std::string tempI = alignmentB;
          std::string tempJ = alignmentA;

          for (int k = i-1; k > ki; k--) {
            tempI = "-" + tempI;

//     --> todo test
//            tempJ = tCSeq[k] + tempJ;
//
            tempJ = tCSeq[k+1] + tempJ;
            c++;
            H_N_index   = k * (qNRes+1) + j;
            if (maxIndices[H_index] == H_N_index) {
              //tempJ = tCSeq[k] + tempJ;
              //tempJ = tCSeq[i] + tempJ;
              found = true;
              break;
            }
          }
          if (found) {
            alignmentB = tempI;
            alignmentA = tempJ;
            //std::cout << " old i " << i;
            i = i - c;
            //std::cout << " new " << i << std::endl;
          }
          else {
            tempI = alignmentB;
            tempJ = alignmentA;
            c = 0;
            for (int k = j-1; k > kj; k--) {
              tempJ = "-" + tempJ;
              tempI = qCSeq[k+1] + tempI;
              c++;
              H_W_index   = i * (qNRes+1) + k;
              if (maxIndices[H_index] == H_W_index) {
                //tempI = qCSeq[k] + tempI;
                //tempI = qCSeq[j] + tempI;
                found = true;
                break;
              }
            }
            if (found) {
              alignmentA = tempJ;
              alignmentB = tempI;
              j = j - c;
            }
            else {
              std::cout << " err " << std::endl;
              //exit(0);
              throw MTKException(" err ");
            }
          }

        }
      }
      while (i > 0) {
        alignmentA = tCSeq[i] + alignmentA;
        alignmentB = "-" + alignmentB;
        i--;
      }
      while (j > 0) {
        alignmentA = "-" + alignmentA;
        alignmentB = qCSeq[j] + alignmentB;
        j--;
      }
    }
    else if (this->algorithmType == 2) { // Smith-Waterman
      // Start at the highest value found in H
      std::cout << " IMPLEMENT ... exiting " << std::endl;
      //exit(1);
      throw MTKException(" IMPLEMENT Smith-Waterman Algorithm ");
    }
    else {
      return 1;
    }

    // Error trap
    if (alignmentA.size() != alignmentB.size()) {
      std::cout << alignmentA << "\n" << alignmentB << std::endl;
      std::cout << alignmentA.size() << "\n" << alignmentB.size() << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function :getSimScore()
// ------------------------------------------------------------
//
// ============================================================
double seqAlign::getSimScore(std::string s, std::string h)
{
    int s_index = 0;
    int h_index = 0;
    for (int j = 0; j < this->pamSize; j++) {
      if (this->pamTypes[j] == s[0]) {
        s_index = j;
        break;
      }
    }

    for (int j = 0; j < this->pamSize; j++) {
      if (this->pamTypes[j] == h[0]) {
        h_index = j;
        break;
      }
    }

    return pam[s_index*this->pamSize + h_index];
}

// ============================================================
// Function :printResults()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::printResults()
{
    std::string errorMessage = " Results \n ";

    errorMessage += " Options used: \n   Algorithm Type = " + i2s(this->algorithmType);
    errorMessage += " \n   Gap Penalty Type = " + i2s(this->gapPenaltyType);
    errorMessage += " \n   Gap Open Value = " + d2s(this->gapOpen);
    errorMessage += " \n   Gap Extend Value = " + d2s(this->gapExtend);
    errorMessage += " \n   Maximum Gap Value = " + i2s(this->maxGap);
    errorMessage += " \n   PAM File = " + this->pamFile;

    int nIdentical = 0;
    int nSimilar = 0;
    int nDissimilar = 0;
    int nGap = 0;
    int stepSize = 65;
    int x = 0;
    int i = 0;

    std::string gapString = "-";
    errorMessage += "\n\n KEY: \n   | == similar\n   : == similar\n   # == dissimilar\n   ^ == gap\n\n ALIGNMENT \n";

    while (i < static_cast<int>(alignmentA.size())) {
      if (x+stepSize < static_cast<int>(alignmentA.size())) {
        x += stepSize;
      }
      else {
        x = alignmentA.size();
      }

      errorMessage += "   ";
      for (int j = i; j < x; j++) {
        errorMessage += alignmentA[j];
      }
      errorMessage += "\n   ";

      for (int j = i; j < x; j++) {
        std::string lA = alignmentA.substr(j,1);
        std::string lB = alignmentB.substr(j,1);
        if (lA == lB) {
          errorMessage += "|";
          nIdentical++;
        }
        else if ((lA != gapString) and
                 (lB != gapString)) {
          double s = getSimScore(lA, lB);
          if (s > -0.00000001) {
            errorMessage += ":";
            nSimilar++;
          }
          else {
            errorMessage += "#";
            nDissimilar++;
          }
        }
        else {
          errorMessage += "^";
          nGap++;
        }
      }
      errorMessage += "\n   ";

      for (int j = i; j < x; j++) {
        errorMessage += alignmentB[j];
      }
      errorMessage += "\n\n";
      i += stepSize;
    }

    double precentageIdentical = (double(nIdentical)/double(alignmentA.size())) * 100.0;
    double precentageSimilar = (double(nIdentical + nSimilar)/double(alignmentA.size())) * 100.0;
    double precentageDissimilar = (double(nDissimilar)/double(alignmentA.size())) * 100.0;
    double precentageGap = (double(nGap)/double(alignmentA.size())) * 100.0;

    errorMessage += " STATS \n Identicals = " + i2s(nIdentical) + "/" + i2s(alignmentA.size()) +
                    " (" + d2s(precentageIdentical) + ") \n";

    errorMessage += " Similars = " + i2s(nIdentical + nSimilar) + "/" + i2s(alignmentA.size()) +
                    " (" + d2s(precentageSimilar) + ") \n";

    errorMessage += " Dissimilars = " + i2s(nDissimilar) + "/" + i2s(alignmentA.size()) +
                    " (" + d2s(precentageDissimilar) + ") \n";

    errorMessage += " Gaps = " + i2s(nGap) + "/" + i2s(alignmentA.size()) +
                    " (" + d2s(precentageGap) + ") \n";
    MTKpp::errorLogger.throwError(" seqAlign::printResults ", errorMessage, 0);
}

// ============================================================
// Function :setPAMFile()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setPAMFile(std::string s)
{
    this->pamFile = s;
}

// ============================================================
// Function :setPAMSize()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setPAMSize(int s)
{
    this->pamSize = s;
    try {
      this->pam = new double [(this->pamSize) * (this->pamSize)];
      this->pamTypes = new char [(this->pamSize)];
    }
    catch (std::bad_alloc) {
      errorLogger.throwError( "seqAlign::setPAMSize", " Memory Allocation Failed ", 1);
      //exit(1);
      std::stringstream ss;
      ss << "seqAlign::setPAMSize Memory Allocation Failed ";
      throw MTKException(ss.str());
    }
}

// ============================================================
// Function :setPAMType()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setPAMType(int i, std::string s)
{
    if (i < this->pamSize and s.size() == 1) {
      this->pamTypes[i] = s[0];
    }
}

// ============================================================
// Function :setPAMSize()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setPAMValue(int i, int j, double v)
{
    if (i < this->pamSize and i < this->pamSize) {
      pam[i*this->pamSize + j] = v;
    }
    else {
      errorLogger.throwError( "seqAlign::setPAMValue", " Out of bounds ", 1);
      //exit(1);
      std::stringstream ss;
      ss << "seqAlign::setPAMValue out-of-bounds";
      throw MTKException(ss.str());
    }
}

// ============================================================
// Function :setAlgorithmType()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setAlgorithmType(int i)
{
    this->algorithmType = i;
}

// ============================================================
// Function :setGapPenaltyType()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setGapPenaltyType(int i)
{
    this->gapPenaltyType = i;
}

// ============================================================
// Function :setGapOpen()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setGapOpen(double d)
{
    this->gapOpen = d;
}

// ============================================================
// Function :setGapExtend()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::setGapExtend(double d)
{
    this->gapExtend = d;
}

// ============================================================
// Function :printPAM()
// ------------------------------------------------------------
//
// ============================================================
void seqAlign::printPAM()
{
    std::string pamString = " \n PAM MATRIX\n";
    for (int j = 0; j < this->pamSize; j++) {
      std::stringstream n;
      n << std::setw(6) << this->pamTypes[j];
      pamString += n.str();
    }
    pamString += "\n";

    for (int i = 0; i < this->pamSize; i++) {
      std::cout << std::setw(6) << this->pamTypes[i] << " ";
      for (int j = 0; j < this->pamSize; j++) {
        std::cout << std::showpoint << std::setw(6) << this->pam[i*(this->pamSize) + j] << " ";
      }
      std::cout << " " << std::endl;
    }

    errorLogger.throwError( "seqAlign::setPAMValue", " Out of bounds ", 1);
    //exit(1);
    std::stringstream ss;
    ss << "seqAlign::pringPAM out-of-bounds";
    throw MTKException(ss.str());
}

// ============================================================
// Function :getCorrMap()
// ------------------------------------------------------------
//
// ============================================================
int* seqAlign::getCorrMap()
{
    return this->corrMap;
}

} // MTKpp namespace
