/*!
   \file functionalize.cpp
   \brief Breaks a molecule into its functional group
   \author Martin Peters

   $Date: 2010/04/29 18:59:18 $
   $Revision: 1.12 $

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

#include "functionalize.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"
#include "utility.h"

#include "fingerPrint.h"
#include "Log/errorHandler.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : functionalize()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
functionalize::functionalize(molecule *parent)
{
    this->pParent = parent;
    pStdGroup    = 0;
    pStdFrag     = 0;
    pStdAtom     = 0;
    pStdBond     = 0;
    pFingerPrint = 0;
}

// ============================================================
// Function : ~functionalize()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
functionalize::~functionalize() {}

// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
int functionalize::run(stdLibrary* pStdLib)
{
    errorLogger.throwError(" functionalize::run ", " Begin ", 4);

    if (!pStdLib) {
      errorLogger.throwError(" functionalize::run ", " Please provide a library ", 1);
      //exit(0);
      std::stringstream ss;
      ss << " functionalize::run "<< " Please provide a library ";
      throw MTKException(ss.str());
    }

    this->pFingerPrint = new fingerPrint();
    if (!this->pFingerPrint) {
      errorLogger.throwError(" functionalize::run ", " Error creating fingerPrint ", 4);
      //exit(0);
      std::stringstream ss;
      ss << " functionalize::run "<< " Error creating fingerPrint ";
      throw MTKException(ss.str());
    }

    std::vector<unsigned int> molSFP = pParent->getSimpleFP();

    int molAtoms = pParent->getNumAtoms();

    char* molAtomSymbols = pParent->getAtomSymbols();
    if (!molAtomSymbols) return 1;

    int molAdjMatSize = pParent->getAdjMatrixSize();
    int* molAdjMatrix = pParent->getAdjMatrix();
    int* molAtomKinds = pParent->getAtomKinds();
    if (!molAdjMatrix) return 1;

    stdGroupList = pStdLib->getStdGroupList();

    int match = 0;

    for (stdGroupIterator g = stdGroupList.begin(); g != stdGroupList.end(); g++) {
      this->pStdGroup = *g;

      stdFragList = pStdGroup->getStdFragList();
      for (stdFragIterator f = stdFragList.begin(); f != stdFragList.end(); f++) {
        this->pStdFrag = *f;
        std::vector<unsigned int> fragSFP = pStdFrag->getSimpleFP();

        if (fragSFP.size() < 105) return -1;
        std::string errMessage = " Check if the " + this->pStdFrag->getName()
                               + " functional group is in the molecule. \n ";

        match = this->pFingerPrint->compareSimpleFP(molSFP, fragSFP);

        if (match == 1) {
          int fragAtoms = this->pStdFrag->numStdAtoms();
          char* fragAtomSymbols = this->pStdFrag->getAtomSymbols();
          int* fragAdjMatrix = this->pStdFrag->getAdjMatrix();
          int fragAdjMatSize = this->pStdFrag->getAdjMatrixSize();
          int* fragAtomKinds = this->pStdFrag->getAtomKinds();

          int* genMatchMatrix;
          try {
            genMatchMatrix   = new int [molAtoms*fragAtoms];
          }
          catch (std::bad_alloc) {
            errorLogger.throwError(" functionalize::run ", " Memory Allocation Failure ", 1);
            return 1;
          }

          this->initializeMatchMatrix(genMatchMatrix,
          molAtoms, molAdjMatSize, molAdjMatrix, molAtomSymbols, molAtomKinds,
          fragAtoms, fragAdjMatSize, fragAdjMatrix, fragAtomSymbols, fragAtomKinds);
/*
#ifdef DEBUG
          std::cout << " functionalize::run\n   genMatchMatrix = "<<std::endl;
          for (int y = 0; y < fragAtoms; y++) {
            for (int x = 0; x < molAtoms; x++) {
              std::cout << genMatchMatrix[y*molAtoms+x] << " ";
            }
            std::cout << " " << std::endl;
          }
          std::cout << " " << std::endl;
#endif
*/
          for (int y = 0; y < fragAtoms; y++) {
            int ms = 0;
            for (int x = 0; x < molAtoms; x++) {
              ms += genMatchMatrix[y*molAtoms+x];
            }
            if (ms == 0) {
              std::stringstream number;
              number << y+1;

              std::stringstream numk;
              numk << fragAtomKinds[y];

              errMessage += " FRAGMENT ATOM " + number.str() + " NOT MATCHED \n";
              errMessage += "  FRAGMENT ATOM KIND " + numk.str() + "\n";
            }
          }

          std::vector<int> subGraph;
          std::vector<std::vector<int> > subGraphs;
          this->ullmann(0, fragAtoms, fragAdjMatrix, fragAtomSymbols,
                        molAtoms, molAdjMatrix, molAtomSymbols,
                        genMatchMatrix, subGraph, subGraphs);

          if (subGraphs.size() < 1) {
            errMessage += "   A fingerprint match but there is no " + pStdFrag->getName() + " "
                       + pStdFrag->getSymbol() + "\n";
          }
          else {
            errMessage += "  Found \n";
          }
          errorLogger.throwError(" functionalize::run ", errMessage, 4);

          for (unsigned int k = 0; k < subGraphs.size(); k++) {
            this->pParent->addFunctionalGroup(fragAtoms, pStdFrag,
                                              molAtoms, subGraphs[k]);
          }

/*
#ifdef DEBUG
          for (unsigned int k = 0; k < subGraphs.size(); k++) {
            std::cout << " functionalize::run \n subgraph : " << k+1 << std::endl;
            subGraph = subGraphs[k];
            for (int t = 0; t < fragAtoms; t++) {
              for (int t2 = 0; t2 < molAtoms; t2++) {
                std::cout << subGraph[t*molAtoms+t2] << " ";
              }
              std::cout << " " << std::endl;
            }
          }
          std::cout << " " << std::endl;
#endif
*/
          delete genMatchMatrix;
        }
        else {
          errMessage += " match = " + i2s(match) + "\n\n";
          errMessage += "      No fingerprint match for " + pStdFrag->getName() + " and "
                     + pStdFrag->getSymbol();
/*
          errMessage += "\n frag fingerprint = \n";
          std::vector<unsigned int> pf = pStdFrag->getSimpleFP();
          for (unsigned int x = 0; x < pf.size(); x++) {
            errMessage += i2s(pf[x]);
          }

          errMessage += "\n molecule fingerprint = \n";
          std::vector<unsigned int> mf = this->pParent->getSimpleFP();
          for (unsigned int x = 0; x < mf.size(); x++) {
            errMessage += i2s(mf[x]);
          }
*/
          errorLogger.throwError(" functionalize::run ", errMessage, 4);
        }
      }
    }
    delete this->pFingerPrint;
    return 0;
}

    //-----------------------//
    // - PRIVATE FUNCTIONS - //
    //-----------------------//

// ============================================================
// Function : initializeMatchMatrix()
// ------------------------------------------------------------
//
// ============================================================
void functionalize::initializeMatchMatrix(int genMatchMatrix[],
                    int molAtoms, int molAdjMatSize, int molAdjMatrix[],
                    char molAtomSymbols[], int molAtomKinds[],
                    int fragAtoms, int fragAdjMatSize,
                    int fragAdjMatrix[], char fragAtomSymbols[],
                    int fragAtomKinds[])
{
    for (int i = 0; i < fragAtoms; i++) {
      for (int j = 0; j < molAtoms; j++) {

        int currAtomKind = molAtomKinds[j];
        if (fragAtomKinds[i] == 7) {
          if (currAtomKind == 3 or currAtomKind == 4) {
            currAtomKind = 7;
          }
        }

        if (fragAtomKinds[i] == 8) {
          if (currAtomKind == 2 or currAtomKind == 3 or currAtomKind == 4) {
            currAtomKind = 8;
          }
        }

        if ((molAtomSymbols[j*2  ] == fragAtomSymbols[i*2  ]) &&
            (molAtomSymbols[j*2+1] == fragAtomSymbols[i*2+1]) &&
            ((currAtomKind == fragAtomKinds[i]) or
             (fragAtomKinds[i] == 0)) ) {

          std::vector<int> matched;
          unsigned int nStdBonds = 0;
          for (int k = 0; k < fragAtoms; k++) {
            if (i == k) continue;
            int pos = i*fragAtoms+k;
            if (fragAdjMatrix[pos]) {
              nStdBonds++;
              for (int l = 0; l < molAtoms; l++) {
                if (j == l) continue;
                int pos2 = j*molAtoms+l;
                if (((molAtomSymbols[l*2  ] == fragAtomSymbols[k*2  ]) &&
                     (molAtomSymbols[l*2+1] == fragAtomSymbols[k*2+1])) &&
                     (
                       (fragAdjMatrix[pos] == molAdjMatrix[pos2]) or
                       ((fragAdjMatrix[pos] == 4) and (molAdjMatrix[pos2] == 6)) or
                       ((fragAdjMatrix[pos] == 4) and (molAdjMatrix[pos2] == 7)) or
                       ((fragAdjMatrix[pos] == 5) and (molAdjMatrix[pos2] == 1)) or
                       ((fragAdjMatrix[pos] == 5) and (molAdjMatrix[pos2] == 2))
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
}

// ============================================================
// Function : ullmann()
// ------------------------------------------------------------
//
// ============================================================
void functionalize::ullmann(int pos, int fragAtoms, int fragAdjMatrix[], char fragAtomSymbols[],
                    int molAtoms, int molAdjMatrix[], char molAtomSymbols[], int matchMatrix[],
                    std::vector<int>& subGraph, std::vector<std::vector<int> >& subGraphs)
{
    // back up match matrix
    int matchMatrixBackUp[fragAtoms*molAtoms];
    for (int i = 0; i < fragAtoms; i++) {
      for (int j = 0; j < molAtoms; j++) {
        matchMatrixBackUp[i*molAtoms+j] = matchMatrix[i*molAtoms+j];
      }
    }

    bool mismatch = false;
    for (int i = 0; i < molAtoms; i++) {
      if (matchMatrix[pos*molAtoms+i]) {
        this->updateMatchMatrix(i, pos, fragAtoms, molAtoms, matchMatrix);
        this->refineMatchMatrix(fragAtoms, fragAdjMatrix, fragAtomSymbols,
              molAtoms, molAdjMatrix, molAtomSymbols, matchMatrix, mismatch);

        if (!mismatch) {
          if (pos == fragAtoms-1) {
            for (int f = 0; f < fragAtoms; f++) {
              for (int g = 0; g < molAtoms; g++) {
                subGraph.push_back(matchMatrix[f*molAtoms+g]);
                matchMatrix[f*molAtoms+g] = matchMatrixBackUp[f*molAtoms+g];
              }
            }
            subGraphs.push_back(subGraph);
            subGraph.clear();
          }
          else {
            ullmann(pos+1, fragAtoms, fragAdjMatrix, fragAtomSymbols,
                    molAtoms, molAdjMatrix, molAtomSymbols, matchMatrix,
                    subGraph, subGraphs);
          }
        }
        else {
          // reset match matrix
          for (int f = 0; f < fragAtoms; f++) {
            for (int j = 0; j < molAtoms; j++) {
              matchMatrix[f*molAtoms+j] = matchMatrixBackUp[f*molAtoms+j];
            }
          }
        }
      }
      // reset match matrix
      for (int f = 0; f < fragAtoms; f++) {
        for (int j = 0; j < molAtoms; j++) {
          matchMatrix[f*molAtoms+j] = matchMatrixBackUp[f*molAtoms+j];
        }
      }
    }
}

// ============================================================
// Function : updateMatchMatrix()
// ------------------------------------------------------------
//
// ============================================================
void functionalize::updateMatchMatrix(int posMolAtom, int posFragAtom, int fragAtoms, int molAtoms, int matchMatrix[])
{
    // zero out the row
    for (int i = 0; i < molAtoms; i++) {
      if (i == posMolAtom) {
        continue;
      }
      else {
        matchMatrix[posFragAtom*molAtoms+i] = 0;
      }
    }

    // zero out the column
    for (int i = 0; i < fragAtoms; i++) {
      if (i == posFragAtom) {
        continue;
      }
      else {
        matchMatrix[i*molAtoms+posMolAtom] = 0;
      }
    }
/*
#ifdef DEBUG
    std::cout << " functionalize::updateMatchMatrix = " << std::endl;
    for (int y = 0; y < fragAtoms; y++) {
      for (int x = 0; x < molAtoms; x++) {
        std::cout << matchMatrix[y*molAtoms+x] << " ";
      }
      std::cout << " " << std::endl;
    }
    std::cout << " " << std::endl;
#endif
*/
}

// ============================================================
// Function : refineMatchMatrix()
// ------------------------------------------------------------
//
// ============================================================
void functionalize::refineMatchMatrix(int fragAtoms, int fragAdjMatrix[], char fragAtomSymbols[],
                                      int molAtoms,  int molAdjMatrix[],  char molAtomSymbols[],
                                      int matchMatrix[], bool &mismatch)
{
    mismatch = false;
    bool change = true;

    while (change and !mismatch) {
      change = false;
      for (int i = 0; i < fragAtoms; i++) {
        for (int j = 0; j < molAtoms; j++) {
          if (matchMatrix[i*molAtoms+j]) {
            std::vector<int> matched;
            unsigned int nStdBonds = 0;
            for (int k = 0; k < fragAtoms; k++) {
              if (i == k) continue;
              int pos = i*fragAtoms+k;
              if (fragAdjMatrix[pos]) {
                nStdBonds++;
                for (int l = 0; l < molAtoms; l++) {
                  if (j == l) continue;
                  int pos2 = j*molAtoms+l;
                  if ( ((molAtomSymbols[l*2  ] == fragAtomSymbols[k*2  ]) &&
                        (molAtomSymbols[l*2+1] == fragAtomSymbols[k*2+1])) &&
                     ( (fragAdjMatrix[pos] == molAdjMatrix[pos2]) or
                       ((fragAdjMatrix[pos] == 4) and (molAdjMatrix[pos2] == 6)) or
                       ((fragAdjMatrix[pos] == 4) and (molAdjMatrix[pos2] == 7)) or
                       ((fragAdjMatrix[pos] == 5) and (molAdjMatrix[pos2] == 1)) or
                       ((fragAdjMatrix[pos] == 5) and (molAdjMatrix[pos2] == 2)) ) &&
                     (matchMatrix[k*molAtoms+l] > 0)) {
                    std::vector<int>::iterator result;
                    result = std::find(matched.begin(), matched.end(), l);
                    if (result == matched.end()) {
                      //std::cout << "   matched i-k" << i+1<<"-"<<k+1 << " j-l" <<j+1<<"-"<<l+1 << std::endl;
                      matched.push_back(l);
                    }
                  }
                }
              }
            }
            // if matched.size > nStdBonds means that a H in the mol is being mapped onto a number of Hs in the frag
            if (nStdBonds != matched.size()) {
              //std::cout << " frag:" << i+1 << " mol:" << j+1 << " nStdBonds = " << nStdBonds << " matched.size() = " << matched.size() << std::endl;
              matchMatrix[i*molAtoms+j] = 0;
              change = true;
            }
          }
        }
      }

      // determine mismatch
      for (int i = 0; i < fragAtoms; i++) {
        bool rowOk = false;
        for (int j = 0; j < molAtoms; j++) {
          if (matchMatrix[i*molAtoms+j] > 0) rowOk = true;
        }
        if (!rowOk) {
          mismatch = true;
        }
      }
    }
/*
#ifdef DEBUG
    std::cout << " functionalize::refineMatchMatrix = " << std::endl;
    for (int y = 0; y < fragAtoms; y++) {
      for (int x = 0; x < molAtoms; x++) {
        std::cout << matchMatrix[y*molAtoms+x] << " ";
      }
      std::cout << " " << std::endl;
    }
    std::cout << " mismatch = " << mismatch << std::endl;
#endif
*/
}



// ============================================================
// Function : determineSymmetricFrags()
// ------------------------------------------------------------
//
// ============================================================
int functionalize::determineSymmetricFrags(stdLibrary* pStdLib)
{
#ifdef DEBUG
      std::cout << "   functionalize::determineSymmetricFrags " << std::endl;
#endif

    if (!pStdLib) {
      std::cout << "    Error in functionalize::determineSymmetricFrags. " << std::endl;
      std::cout << "     Please provide a library ... exiting " << std::endl;
      return 1;
    }

    if (!pStdLib->hasSimpleFPGenerated() or
        !pStdLib->hasAdjMatricesGenerated() or
        !pStdLib->hasAtomKindsGenerated()) {
      std::cout << "    Error in functionalize::determineSymmetricFrags. ";
      std::cout << "... exiting " << std::endl;
      return 1;
    }

    stdGroupList = pStdLib->getStdGroupList();
    for (stdGroupIterator g = stdGroupList.begin(); g != stdGroupList.end(); g++) {
      this->pStdGroup = *g;
      stdFragList = pStdGroup->getStdFragList();
      for (stdFragIterator f = stdFragList.begin(); f != stdFragList.end(); f++) {
        this->pStdFrag = *f;
        if (this->pStdFrag->getNumStdConnPtsList() != 1) continue;

        int   fragAtoms       = this->pStdFrag->numStdAtoms();
        char* fragAtomSymbols = this->pStdFrag->getAtomSymbols();
        int*  fragAdjMatrix   = this->pStdFrag->getAdjMatrix();
        int   fragAdjMatSize  = this->pStdFrag->getAdjMatrixSize();
        int*  fragAtomKinds   = this->pStdFrag->getAtomKinds();

/*
        int fragAtoms = this->pStdFrag->numStdHeavyAtoms();
        int f = this->pStdFrag->generateHeavyAtomSymbols();
        if (f) {
          std::stringstream ss;
          ss << "   Error in functionalize::determineSymmetricFrags. ";
          ss << " ... exiting " << std::endl;
          std::cout << ss.str();
          throw MTKException(ss.str());
        }
        char* fragAtomSymbols = this->pStdFrag->getHeavyAtomSymbols();
        f = this->pStdFrag->generateHeavyAdjMatrix();
        if (f) {
          std::stringstream ss;
          ss << "   Error in functionalize::determineSymmetricFrags. ";
          ss << " ... exiting " << std::endl;
          std::cout << ss.str();
          throw MTKException(ss.str());
        }
        int* fragAdjMatrix = this->pStdFrag->getHeavyAdjMatrix();
        int fragAdjMatSize = this->pStdFrag->getHeavyAdjMatrixSize();
        f = this->pStdFrag->generateHeavyAtomKinds();
        if (f) {
          std::stringstream ss;
          ss << "   Error in functionalize::determineSymmetricFrags. ";
          ss << " ... exiting " << std::endl;
          std::cout << ss.str();
          throw MTKException(ss.str());
        }
        int* fragAtomKinds = this->pStdFrag->getHeavyAtomKinds();
*/
        int* genMatchMatrix;
        try {
          genMatchMatrix   = new int [fragAtoms*fragAtoms];
        }
        catch (std::bad_alloc) {
          std::cout << " Memory Allocation Failure " << std::endl;
          return 1;
        }

        // Initial Match Matrix
        this->initializeMatchMatrix(genMatchMatrix,
        fragAtoms, fragAdjMatSize, fragAdjMatrix, fragAtomSymbols, fragAtomKinds,
        fragAtoms, fragAdjMatSize, fragAdjMatrix, fragAtomSymbols, fragAtomKinds);

        // Run Ullmann
        std::vector<int> subGraph;
        std::vector<std::vector<int> > subGraphs;
        this->ullmann(0, fragAtoms, fragAdjMatrix, fragAtomSymbols,
                         fragAtoms, fragAdjMatrix, fragAtomSymbols,
                         genMatchMatrix, subGraph, subGraphs);

        if (subGraphs.size() >= 2) {
          std::cout  << "    " << this->pStdFrag->getName() << " is a symmetric functional group. "
                     << " Its got " << subGraphs.size() << " subgraphs. " << std::endl;
        }
        else {
          std::cout << "     " << this->pStdFrag->getName() << " has got "
                    << subGraphs.size() << " subgraphs. " << std::endl;
        }
        delete genMatchMatrix;
      }
    }
    return 0;
}

} // MTKpp namespace
