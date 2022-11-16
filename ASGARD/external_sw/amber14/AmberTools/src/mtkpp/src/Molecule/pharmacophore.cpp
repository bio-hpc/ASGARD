/*!
   \file pharmacophore.cpp
   \brief Determines the Maximum Common Pharmacophore (MCP) between two molecules
   \author Martin Peters

   $Date: 2007/09/14 11:35:24 $
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

#include "pharmacophore.h"

#include "molecule.h"
#include "atom.h"

namespace MTKpp
{

// ============================================================
// Function : clique()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
clique::clique()
{
    this->score = 0.0;
    this->dMax = 0.0;
}

// ============================================================
// Function : ~clique()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
clique::~clique() {}

// ============================================================
// Function : compare()
// ------------------------------------------------------------
// Compares two cliques
// ============================================================
bool clique::compare(clique* rhs)
{
    unsigned int clqSize = this->indicesA.size();
    //unsigned int counter = 0;
    if (clqSize == rhs->indicesA.size()) {
/*
      for (unsigned int i = 0; i < clqSize; i++) {
        counter = 0;
        for (unsigned int j = 0; j < clqSize; j++) {
          if (this->indicesA[i] == rhs->indicesA[j]) {
            if (this->indicesB[i] != rhs->indicesB[j]) {
              return true;
            }
          }
          else {
            counter++;
          }
        }
        if (counter == clqSize) return true;
      }
*/

      unsigned int nHits = 0;
      for (unsigned int i = 0; i < clqSize; i++) {
        for (unsigned int j = 0; j < clqSize; j++) {
          double coordDiff = this->coords[i].dist(rhs->coords[j]);
          if (std::abs(coordDiff) < 0.0001) {
            nHits++;
            break;
          }
        }
      }
      if (nHits != clqSize) return true;

    }
    else {
      return true;
    }
    return false;
}

// ============================================================
// Function : pharmacophore()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pharmacophore::pharmacophore(molecule *parent, double distMax)
{
    this->pParent = parent;
    this->deltaDistMax = distMax*distMax;
}

// ============================================================
// Function : ~pharmacophore()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
pharmacophore::~pharmacophore() {}

// ============================================================
// Function : run()
// ------------------------------------------------------------
// 
// ============================================================
int pharmacophore::run(molecule* pMol,
                   std::vector<std::vector<unsigned int> > &mcp,
                   std::vector<vector3d> &mcpCoords)
{
/*
#ifdef DEBUG
      std::cout << "   pharmacophore::run " << std::endl;
#endif
*/
    if (!pMol) {
      std::cout << "   Error in pharmacophore::run ... returning " << std::endl;
      return 1;
    }

    std::vector<std::string> featureLabelsA = this->pParent->getFeatureLabels();
    std::vector<vector3d> featureCoordsA = this->pParent->getFeatureCoordinates();
    unsigned int nFeaturesA = featureLabelsA.size();
    double* featureDistMatrixA = pParent->getFeatureDistMatrix();

    std::vector<std::string> featureLabelsB = pMol->getFeatureLabels();
    std::vector<vector3d> featureCoordsB = pMol->getFeatureCoordinates();
    unsigned int nFeaturesB = featureLabelsB.size();
    double* featureDistMatrixB = pMol->getFeatureDistMatrix();

    if ((featureLabelsA.size() < 3) or (featureLabelsB.size() < 3)) {
      std::cout << "   Error in pharmacophore::run ... returning " << std::endl;
      return 1;
    }

    // Generate the correspondence matrix between molecule A and B
    int* crpdeMatrix;
    try {
      crpdeMatrix   = new int [nFeaturesA*nFeaturesB];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    int failure = this->generateCorrMatrix(nFeaturesA, nFeaturesB,
                  featureLabelsA, featureLabelsB, crpdeMatrix);
    if (failure) {
      std::cout << "   Error generating correspondence matrix ... returning " << std::endl;
      return 1;
    }

    double cGood = 0.0;
    failure = this->getCGood(cGood, nFeaturesA, nFeaturesB,
                    featureLabelsA, featureDistMatrixA,
                    featureDistMatrixB, crpdeMatrix);
    if (failure) {
      std::cout << "   Error determining cGood ... returning " << std::endl;
      return 1;
    }

    // Starting with each possible pairing as a seed, add other pairs
    // such that the intramolecular feature distances are consistent
    // between the two structures

    std::vector<std::vector<std::vector<unsigned int> > > cliques;
    std::vector<double> scores;
    std::vector<double> scores2;
    std::vector<double> maxDistances;
    double score = 0.0;
    double score2 = 0.0;
    double dMax = 0.0;

    double bestScore = 0.0;
    unsigned int bestSize = 0;
    double bestDist = 0.0;
    int bestClique = 0;

    for (unsigned int i = 0; i < nFeaturesA; i++) {
      for (unsigned int j = 0; j < nFeaturesB; j++) {
        if (crpdeMatrix[i*nFeaturesB+j]) {

          //std::cout << " TESTING " << i << " " << j << " FROM CORRESPONDENCE MATRIX " << std::endl;
          std::vector<std::vector<unsigned int> > clique;
          failure = this->getClique(i, j, nFeaturesA, nFeaturesB,
                          featureLabelsA,
                          featureDistMatrixA, featureDistMatrixB,
                          crpdeMatrix, cGood, clique, score, score2, dMax);
          if (failure) {
            std::cout << "   Error in getClique ... returning " << std::endl;
            return 1;
          }
          if (clique.size() > 2) {
            cliques.push_back(clique);
            scores.push_back(score);
            scores2.push_back(score2);
            maxDistances.push_back(dMax);
          }

          int store = 0;
          // Store best clique
          if (score > bestScore+0.1) {
            store = 1;
          }
          else if (score > bestScore-0.1) {
            // score is more or less the same. use additional criteria
            if (clique.size() > bestSize) {
              store = 1;
            }
            else if (clique.size() == bestSize) {
              if (dMax > bestDist+0.1) {
                store = 1;
              }
/*
              else if (dMax == bestDist) {
                std::cout << "\n\n size and dMax are the same " << std::endl;
                for (unsigned int x = 0 ; x < bestSize; x++) {
                 std::cout << "     " << cliques[bestClique][x][0] << ":" << cliques[bestClique][x][1] << std::endl;
                }
                int t = cliques.size()-1;
                for (unsigned int x = 0 ; x < bestSize; x++) {
                 std::cout << "     " << cliques[t][x][0] << ":" << cliques[t][x][1] << std::endl;
                }
                std::cout << "\n\n\n" << std::endl;
              }
*/
            }
          }
          if (store) {
            bestScore = score;
            bestSize = clique.size();
            bestDist = dMax;
            bestClique = cliques.size()-1;
          }
        }
      }
    }
/*
#ifdef DEBUG
    std::cout << "\n\n MOL FEATURES A: ";
    for (unsigned int i = 0; i < nFeaturesA; i++) {
      std::cout << featureLabelsA[i] << " ";
    }
    std::cout << " " << std::endl;

    std::cout << " MOL FEATURES B: ";
    for (unsigned int i = 0; i < nFeaturesB; i++) {
      std::cout << featureLabelsB[i] << " ";
    }
    std::cout << " " << std::endl;

    std::cout << " CLIQUES " << std::endl;
    for (unsigned int x = 0 ; x < cliques.size(); x++) {
      for (unsigned int y = 0 ; y < cliques[x].size(); y++) {
        std::cout << "     " << cliques[x][y][0] << ":" << cliques[x][y][1] << std::endl;
      }
      std::cout << " score = " << scores[x] << " maxDist = " << maxDistances[x] << std::endl;
    }

    std::cout << " MCP RESULTS2" << std::endl;
    std::cout << "    BEST SCORE : " << bestScore << std::endl;
    std::cout << "    BEST SIZE  : " << bestSize << std::endl;
    std::cout << "    CLIQUE: " << std::endl;
    for (unsigned int x = 0 ; x < bestSize; x++) {
      std::cout << "     " << cliques[bestClique][x][0] << ":" << cliques[bestClique][x][1] << std::endl;
    }
#endif
*/
    for (unsigned int x = 0 ; x < bestSize; x++) {
      mcpCoords.push_back(featureCoordsA[cliques[bestClique][x][0]]);
      mcpCoords.push_back(featureCoordsB[cliques[bestClique][x][1]]);

      std::vector<unsigned int> pair;
      pair.push_back(cliques[bestClique][x][0]);
      pair.push_back(cliques[bestClique][x][1]);
      mcp.push_back(pair);
    }
    return 0;
}

// ============================================================
// Function : getCliques()
// ------------------------------------------------------------
// 
// ============================================================
int pharmacophore::getCliques(molecule* pMol, std::vector<clique*> &cliqueList)
{
/*
#ifdef DEBUG
      std::cout << "    pharmacophore::getCliques " << std::endl;
#endif
*/
    if (!pMol) {
      std::cout << "     Error in getCliques ... returning " << std::endl;
      return 1;
    }

    std::vector<std::string> featureLabelsA = this->pParent->getFeatureLabels();
    std::vector<vector3d> featureCoordsA = this->pParent->getFeatureCoordinates();
    unsigned int nFeaturesA = featureLabelsA.size();
    double* featureDistMatrixA = pParent->getFeatureDistMatrix();

    std::vector<std::string> featureLabelsB = pMol->getFeatureLabels();
    std::vector<vector3d> featureCoordsB = pMol->getFeatureCoordinates();
    unsigned int nFeaturesB = featureLabelsB.size();
    double* featureDistMatrixB = pMol->getFeatureDistMatrix();

    if ((featureLabelsA.size() < 3) or (featureLabelsB.size() < 3)) {
      std::cout << "     Error in getCliques ... returning " << std::endl;
      return 1;
    }

    // Generate the correspondence matrix between molecule A and B
    int* crpdeMatrix;
    try {
      crpdeMatrix   = new int [nFeaturesA*nFeaturesB];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    int failure = this->generateCorrMatrix(nFeaturesA, nFeaturesB,
                  featureLabelsA, featureLabelsB, crpdeMatrix);
    if (failure) {
      std::cout << "     Error in generateCorrMatrix ... returning " << std::endl;
      return 1;
    }

    double cGood = 0.0;
    failure = this->getCGood(cGood, nFeaturesA, nFeaturesB,
                    featureLabelsA, featureDistMatrixA,
                    featureDistMatrixB, crpdeMatrix);
    if (failure) {
      std::cout << "     Error in getCGood ... returning " << std::endl;
      return 1;
    }

    // Starting with each possible pairing as a seed, add other pairs
    // such that the intramolecular feature distances are consistent 
    // between the two structures

    std::vector<std::vector<std::vector<unsigned int> > > cliques;
    std::vector<double> scores;
    std::vector<double> scores2;
    std::vector<double> maxDistances;
    double score = 0.0;
    double score2 = 0.0;
    double dMax = 0.0;

    double bestScore = 0.0;
    unsigned int bestSize = 0;
    double bestDist = 0.0;
    int bestClique = 0;
    int nKept = 0;

    for (unsigned int i = 0; i < nFeaturesA; i++) {
      for (unsigned int j = 0; j < nFeaturesB; j++) {
        if (crpdeMatrix[i*nFeaturesB+j]) {
          std::vector<std::vector<unsigned int> > ijClique;
          failure = this->getClique(i, j, nFeaturesA, nFeaturesB,
                          featureLabelsA,
                          featureDistMatrixA, featureDistMatrixB,
                          crpdeMatrix, cGood, ijClique, score, score2, dMax);
          if (failure) {
            std::cout << " Error in getClique ... returning " << std::endl;
            return 1;
          }
          if (ijClique.size() > 2) {
            cliques.push_back(ijClique);
            scores.push_back(score);
            scores2.push_back(score2);
            maxDistances.push_back(dMax);

            clique* clq = new clique();
            clq->score = score;
            clq->score2 = score2;

            std::vector<vector3d> curCliqueCoords;
            int cClique = cliques.size()-1;
            for (unsigned int x = 0 ; x < ijClique.size(); x++) {
              curCliqueCoords.push_back(featureCoordsA[cliques[cClique][x][0]]);
              clq->indicesA.push_back(cliques[cClique][x][0]);
              curCliqueCoords.push_back(featureCoordsB[cliques[cClique][x][1]]);
              clq->indicesB.push_back(cliques[cClique][x][1]);
            }

            clq->coords = curCliqueCoords;
            bool gotIt = false;
            for (unsigned int d = 0; d < cliqueList.size(); d++) {
              if (*clq == *cliqueList[d]) gotIt = true;
              if (gotIt) break;
            }

            if (!gotIt) {
              cliqueList.push_back(clq);
              nKept++;
            }
            else {
              delete clq;
            }
          }

          int store = 0;
          // Store best clique
          if (score > bestScore+0.1) {
            store = 1;
          }
          else if (score > bestScore-0.1) {
            // score is more or less the same. use additional criteria
            if (ijClique.size() > bestSize) {
              store = 1;
            }
            else if (ijClique.size() == bestSize) {
              if (dMax > bestDist+0.1) {
                store = 1;
              }
            }
          }
          if (store) {
            bestScore = score;
            bestSize = ijClique.size();
            bestDist = dMax;
            bestClique = cliques.size()-1;
          }
        }
      }
    }

    if (nKept == 0) {
      std::cout << "     No cliques found, increase maxDist ... returning " << std::endl;
      return 1;
    }
/*
#ifdef DEBUG
    std::cout << "\n\n MOL FEATURES A: ";
    for (unsigned int i = 0; i < nFeaturesA; i++) {
      std::cout << featureLabelsA[i] << " ";
    }
    std::cout << " " << std::endl;

    std::cout << " MOL FEATURES B: ";
    for (unsigned int i = 0; i < nFeaturesB; i++) {
      std::cout << featureLabelsB[i] << " ";
    }
    std::cout << " " << std::endl;

    std::cout << " CLIQUES " << std::endl;
    for (unsigned int x = 0 ; x < cliques.size(); x++) {
      for (unsigned int y = 0 ; y < cliques[x].size(); y++) {
        std::cout << "     " << cliques[x][y][0] << ":" 
                  << cliques[x][y][1] << std::endl;
      }
      std::cout << " score = " << scores[x] << " score2 = " << scores2[x]
                << " maxDist = " << maxDistances[x] << std::endl;
    }

    std::cout << " MCP RESULTS1" << std::endl;
    std::cout << "    BEST SCORE : " << bestScore << std::endl;
    std::cout << "    BEST SIZE  : " << bestSize << std::endl;
    std::cout << "    CLIQUE: " << std::endl;
    for (unsigned int x = 0 ; x < bestSize; x++) {
      std::cout << "     " << cliques[bestClique][x][0]
                << ":" << cliques[bestClique][x][1] << std::endl;
    }
#endif
*/
    return 0;
}

    //-----------------------//
    // - PRIVATE FUNCTIONS - //
    //-----------------------//

// ============================================================
// Function : generateCorrMatrix()
// ------------------------------------------------------------
// 
// ============================================================
int pharmacophore::generateCorrMatrix(unsigned int nFeaturesA, unsigned int nFeaturesB,
                                       std::vector<std::string> featureLabelsA,
                                       std::vector<std::string> featureLabelsB,
                                       int crpdeMatrix[])
{
//    std::cout << "\n   pharmacophore::Correspondence Matrix:" << std::endl;
    for (unsigned int i = 0; i < nFeaturesA; i++) {
      for (unsigned int j = 0; j < nFeaturesB; j++) {
        if (featureLabelsA[i] == featureLabelsB[j]) {
          crpdeMatrix[i*nFeaturesB+j] = 1;
        }
        else {
          crpdeMatrix[i*nFeaturesB+j] = 0;
        }
        //std::cout << crpdeMatrix[i*nFeaturesB+j] << " ";
      }
      //std::cout << " " << std::endl;
    }
    return 0;
}

// ============================================================
// Function : getCGood()
// ------------------------------------------------------------
// 
// ============================================================
int pharmacophore::getCGood(double &cGood,
                   unsigned int nFeaturesA, unsigned int nFeaturesB,
                   std::vector<std::string> featureLabelsA,
                   double* featureDistMatrixA, double* featureDistMatrixB,
                   int crpdeMatrix[])
{
    double distA = 0.0;
    double distB = 0.0;
    double deltaDist = 0.0;
    double cBest = 0.0;
    double cScore = 0.0;
    double cAvg = 0.0;
    int nScore = 0;
    int nPairs = 0;

    for (unsigned int i = 0; i < nFeaturesA; i++) {
      for (unsigned int j = 0; j < nFeaturesB; j++) {
        cScore = 0.0;
        nScore = 0;

        if (crpdeMatrix[i*nFeaturesB+j]) {
          std::string f = featureLabelsA[i];
          nPairs++;
          for (unsigned int k = 0; k < nFeaturesA; k++) {
            if (k == i) continue;
            for (unsigned int l = 0; l < nFeaturesB; l++) {
              if (l == j) continue;
              if (crpdeMatrix[k*nFeaturesB+l]) {
                distA = featureDistMatrixA[k*nFeaturesA+i];
                distB = featureDistMatrixB[l*nFeaturesB+j];
                if ((distA < 0.1) or (distB < 0.1)) continue;
                deltaDist = pow((distA - distB),2);
                //std::cout << "     k=" << k << " i=" << i << " l=" 
                //                       << l << " j=" << j
                //          << " distKI = " << distA << " distLJ = " << distB
                //          << " cScore = " << exp(-deltaDist/deltaDistMax) << std::endl;

                cScore += exp(-deltaDist/deltaDistMax);

                nScore++;
              }
            }
          }
        }
        if (nScore != 0) cScore = cScore/nScore;
        cBest = std::max(cBest, cScore);
        cAvg += cScore;
      }
    }
    cAvg = cAvg / nPairs;
    cGood = 0.5*(cBest + cAvg);
    return 0;
}

// ============================================================
// Function : getClique()
// ------------------------------------------------------------
// 
// ============================================================
int pharmacophore::getClique(unsigned int iA, unsigned int jB,
                   unsigned int nFeaturesA, unsigned int nFeaturesB,
                   std::vector<std::string> featureLabelsA,
                   double* featureDistMatrixA, double* featureDistMatrixB,
                   int crpdeMatrix[], double cGood,
                   std::vector<std::vector<unsigned int> > &clique,
                   double &score, double &score2, double &dMax)
{
    bool getPair = true;
    double distA = 0.0;
    double distB = 0.0;
    double deltaDist = 0.0;
    double cScore = 0.0;
    double cTest = 0.0;
    double cMax = 0.0;
    double cTotal = 0.0;
    unsigned int iMax = nFeaturesA+1;
    unsigned int jMax = nFeaturesB+1;
    score = 0.0;
    dMax = 0.0;

    double cScore2 = 0.0;
    double cTest2 = 0.0;
    double cMax2 = 0.0;
    double cTotal2 = 0.0;

    clique.clear();
    std::vector<unsigned int> pair;
    pair.push_back(iA);
    pair.push_back(jB);
    clique.push_back(pair);

    while (getPair) {
      getPair = false;
      cMax = 0.0;
      iMax = nFeaturesA+1;
      jMax = nFeaturesB+1;
      for (unsigned int i = 0; i < nFeaturesA; i++) {
        if (i == iA) continue;
        for (unsigned int j = 0; j < nFeaturesB; j++) {
          if (j == jB) continue;
          if (crpdeMatrix[i*nFeaturesB+j]) {

            std::string f = featureLabelsA[i];
            cTest = 0.0;

            for (unsigned int k = 0; k < clique.size(); k++) {
              distA = featureDistMatrixA[clique[k][0]*nFeaturesA+i];
              distB = featureDistMatrixB[clique[k][1]*nFeaturesB+j];

              if ((distA < 0.1) or (distB < 0.1)) {
                cScore = 0.0;
                cTest = 0.0;
                break;
              }

              deltaDist = pow((distA - distB), 2);
              cScore = exp(-deltaDist/deltaDistMax);
              //cScore2 = deltaDist;
              cScore2 = std::abs(distA - distB);

              //std::cout << "      " << i << ":" << clique[k][0] << " distA = " << distA
              //          << "      " << j << ":" << clique[k][1] << " distB = " << distB;
              //std::cout << " cGood = " << cGood << " cScore = " << cScore << std::endl;

              if (cScore < cGood) {
                cScore = 0.0;
                cTest = 0.0;
                cScore2 = 0.0;
                cTest2 = 0.0;
                break;
              }
              cTest += cScore;
              cTest2 += cScore2;
            }
            if (cTest > cMax) {
              cMax = cTest;
              cMax2 = cTest2;
              iMax = i;
              jMax = j;
            }
          }
        }
      }

      if ((iMax < nFeaturesA+1) and (jMax < nFeaturesB+1)) {
        bool gotIt = false;
        for (unsigned int x = 0; x < clique.size(); x++) {
          if ((clique[x][0] == iMax) and (clique[x][1] == jMax)) gotIt = true;
        }
        if (!gotIt) {
          std::vector<unsigned int> newPair;
          newPair.push_back(iMax);
          newPair.push_back(jMax);
          clique.push_back(newPair);
          cTotal += cMax;
          cTotal2 += cMax2;
          getPair = true;
        }
      }
    }
    score = cTotal;

    int cSize = clique.size();

    int nDists = 0;
    for (unsigned int x = cSize-1; x > 0 ; x--) {
      nDists += x;
    }

    if (std::abs(nDists - cTotal) < 0.001) {
      score2 = double(-cSize);
    }
    else if (cTotal2 < 0.00001) {
      score2 = double(-cSize);
    }
    else {
      score2 = cTotal2;
    }
/*
#ifdef DEBUG
    std::cout << " CLIQUE SIZE = " << clique.size() << " cTotal = "
              << cTotal << " cTotal2 = " << cTotal2 << std::endl;
    for (unsigned int x = 0; x < clique.size(); x++) {
      std::cout << clique[x][0] << " <--Mapped To--> " << clique[x][1] << std::endl;
    }
#endif
*/
//    score2 = score2 / (double(clique.size())/double(std::min(nFeaturesA, nFeaturesB)));

    for (unsigned int x = 0; x < clique.size(); x++) {
      for (unsigned int y = x+1; y < clique.size(); y++) {
        double dA = featureDistMatrixA[clique[x][0]*nFeaturesA+clique[y][0]];
        if (dA > dMax) {
          dMax = dA;
        }
      }
    }
    return 0;
}

} // MTKpp namespace
