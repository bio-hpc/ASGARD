/*!
   \file hybridize.cpp
   \brief Determines hybridizations of atoms in a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:43:22 $
   $Revision: 1.18 $

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

#include "hybridize.h"

#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"
#include "Utils/vector3d.h"
#include "Utils/index.h"
#include "utility.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

// Graph
#include "Graph/graph.h"
#include "Graph/edge.h"
#include "Graph/vertex.h"

// Boost
// #include "Utils/diagonalize.h"

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

#include "Utils/diagonalize_eigen.h"

namespace MTKpp
{

// ============================================================
// Function : hybridize()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
hybridize::hybridize(molecule *parent, std::string f)
{
    this->pParent = parent;
    this->parameterFile = f;

    //errorMessage = "MAXATOMS: " + i2s(MAXATOMS);
    //errorLogger.throwError("hybridize::construct", errorMessage, MTK_ERROR);

    atomList = pParent->getAtomList();
    nAtoms = atomList.size();
    dim.assign(nAtoms,0);
    B.assign(nAtoms,0);
    Q.assign(nAtoms, 0);
    hybridizations.assign(nAtoms, 0);
    planarFlags.assign(nAtoms, 0);

    for (unsigned int i = 0; i < nAtoms; i++) {
      atNumbers.push_back(atomList[i]->getAtomicNum());
      atGroups.push_back(atomList[i]->getElement()->group);
      atPeriods.push_back(atomList[i]->getElement()->period);
      atSymbols.push_back(atomList[i]->getElementSymbol());
      atENs.push_back(atomList[i]->getElement()->paulingEN);
    }
}

// ============================================================
// Function : ~hybridize()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
hybridize::~hybridize() {}

// ============================================================
// Function : readParameters()
// ------------------------------------------------------------
//
// ============================================================
int hybridize::readParameters()
{
    // check if file exists
    std::fstream fin;
    fin.open(this->parameterFile.c_str(), std::ios::in);

    if (!fin.is_open()) {
      fin.close();
      return 1;
    }
    fin.close();

    std::ifstream iPF;
    iPF.open(this->parameterFile.c_str());
    if (!iPF) {
      errorMessage = "Unable to open parameter file: " + this->parameterFile;
      errorLogger.throwError("hybridize::readParameters", errorMessage, 1);
      return 1;
    }

    std::string message = " \n Reading Labute Parameter File: " + this->parameterFile + " \n";

    std::string fileline = "";
    std::string heading = "";

    // READ THE FILE
    while (iPF) {
      getline(iPF, fileline);
      heading = "SMALLEST OUT-OF-PLANE DIHEDRAL";
      if (fileline.substr(0, heading.size()) == heading) {
        message += heading + ": \n";
        getline(iPF, fileline);
        std::vector<std::string> words;
        lSplitString(fileline, " ", words, 0);
        this->smallOutOfPlaneDihedral = strtod(words[0].c_str(), 0);
        message += " " + words[0] + " \n";
      }

      heading = "SINGLE BONDS";
      if (fileline.substr(0, heading.size()) == heading) {
        message += heading + ": \n";
        getline(iPF, fileline);
        while (fileline != "") {
          std::vector<std::string> words;
          lSplitString(fileline, " ", words, 0);
          double dist = strtod(words[2].c_str(), 0);
          std::string s1 = words[0]+"-"+words[1];
          std::string s2 = words[1]+"-"+words[0];
          singleBonds[s1] = dist;
          singleBonds[s2] = dist;
          message += " " + s1 + " = " + words[2] + " \n";
          words.clear();
          getline(iPF, fileline);
        }
      }

      heading = "LOG LIKELIHOOD RATIOS";
      if (fileline.substr(0, heading.size()) == heading) {
        message += heading + ": \n";
        getline(iPF, fileline);
        while (fileline != "") {
          std::vector<std::string> words;
          lSplitString(fileline, " ", words, 0);
          std::vector<double> values;
          values.push_back(strtod(words[3].c_str(), 0));
          values.push_back(strtod(words[4].c_str(), 0));
          values.push_back(strtod(words[5].c_str(), 0));

          std::string lineName = "";
          if (words[1] != "X" and words[2] != "X") {
            lineName = words[0]+"-"+words[1]+"-"+words[2];
          }
          else if (words[1] != "X") {
            lineName = words[0]+"-"+words[1];
          }
          else {
            lineName = words[0];
          }

          logLikelihoodRatios[lineName] = values;
          message += " " + lineName + " = [" + words[3] + "," + words[4] + "," + words[5] + "]\n";
          words.clear();
          getline(iPF, fileline);
        }
      }

      heading = "EDGE WEIGHT PARAMETERS";
      if (fileline.substr(0, heading.size()) == heading) {
        message += heading + ": \n";
        getline(iPF, fileline);
        std::vector<std::string> words;
        lSplitString(fileline, " ", words, 0);
        this->edgeWeightSingle = strtod(words[0].c_str(), 0);
        this->edgeWeightDouble = strtod(words[1].c_str(), 0);
        message += " " + words[0] + "," + words[1] + "\n";
      }

      heading = "TRIPLE BOND PARAMETER";
      if (fileline.substr(0, heading.size()) == heading) {
        message += heading + ": \n";
        getline(iPF, fileline);
        std::vector<std::string> words;
        lSplitString(fileline, " ", words, 0);
        this->tripleBondParameter = strtod(words[0].c_str(), 0);
        message += " " + words[0];
      }
    }
    iPF.close();

    errorLogger.throwError("hybridize::readParameters", message, 4);

    return 0;
}

//=============================================================
// Function : lSplitString
// ------------------------------------------------------------
// splits up a string based on a separator and returns a vector
// ============================================================
void hybridize::lSplitString(std::string &text, const std::string separator,
                 std::vector<std::string> &words, int mystart)
{
    int n = text.length();
    int start, stop;

    start = text.find_first_not_of(separator, mystart);

    while ((start >= 0) && (start < n)) {
      stop = text.find_first_of(separator, start);
      if ((stop < 0) || (stop > n)) {
        stop = n;
      }

      words.push_back(text.substr(start, stop-start));
      start = text.find_first_not_of(separator, stop+1);
    }
}

// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
int hybridize::run()
{
    int failure = 0;

    bonds = pParent->getBondMap();

    failure = this->runLabute1();
    if (failure) return 1;

    failure = this->runLabute2();
    if (failure) return 1;

    failure = this->runLabute3();
    if (failure) return 1;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      if (!hybridizations[i]) {
        int bondOrder = 0;

        std::vector<atom*> bondedAtoms = pAtom->getBondedAtoms();
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          atom* pAt2 = bondedAtoms[j];
          pBond = pParent->getBond(pAtom, pAt2);
          if (pBond) {
            if (pBond->type > bondOrder) {
              bondOrder = pBond->type;
            }
          }
        }
        if (bondOrder == 1) {
          hybridizations[i] = 4;
        }
        else if (bondOrder == 2) {
          hybridizations[i] = 3;
        }
        else if (bondOrder == 3) {
          hybridizations[i] = 2;
        }
      }
    }
    this->setHybridizations();
    return 0;
}

// ============================================================
// Function : runLabute()
// ------------------------------------------------------------
//
// ============================================================
int hybridize::runLabute1()
{
    vector3d* com = new vector3d(0.0);

/* 
    // Boost
    ublas::matrix<double, ublas::column_major> gramMatrix(3,3);
    ublas::vector<double> R(3);
*/
    // Eigen
    Eigen::MatrixXd gramMatrix(3, 3);
    
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      pAtom->setHybridization(0); // initialize to undefined
      std::vector<atom*> bondedAtoms = pAtom->getBondedAtoms();
      std::vector<atom*> clusterAtoms;
      clusterAtoms.push_back(pAtom);
      for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
        clusterAtoms.push_back(bondedAtoms[j]);
      }

      if (bondedAtoms.size() > 1) {
        this->COM(clusterAtoms, com);

        for (unsigned int k = 0; k < 3; k++) {
          //R(k) = 0.0;
          for (unsigned int l = 0; l < 3; l++) {
            gramMatrix(k,l) = 0.0;
          }
        }

        for (unsigned int j = 0; j < clusterAtoms.size(); j++) {
          vector3d sub = (*clusterAtoms[j]->getCoords()) - (*com);
          for (unsigned int k = 0; k < 3; k++) {
            for (unsigned int l = 0; l < 3; l++) {
              gramMatrix(k,l) += sub[k] * sub[l];
            }
          }
        }
/*
        // Boost
        int r = diagonalize(gramMatrix,R);
        if (r != 0) return 1;
*/
        SelfAdjointEigenSolver<MatrixXd> eigensolver(gramMatrix);
        //MatrixXd evectors = eigensolver.eigenvectors();
        VectorXd R = eigensolver.eigenvalues();
        // eigenValueSort(evectors, evalues, 1);

        /*
        std::cout << " gramMatrix \n" << gramMatrix << std::endl;
        std::cout << " R \n" << R << std::endl;
        */
        for (unsigned int j = 0; j < 3; j++) {
          if (( R(j) > 0.0 ) and ( sqrt(R(j)) > 0.2 )) {
            dim[i]++;
          }
        }
      }
      else {
        dim[i] = bondedAtoms.size();
      }
    }
    delete com;
    return 0;
}

// ============================================================
// Function : runLabute2()
// ------------------------------------------------------------
//
// ============================================================
int hybridize::runLabute2()
{
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      int atNumber = pAtom->getAtomicNum();
      if (atNumber < 3) { // Hydrogen and Helium
        B[i] = 1;
        continue;
      }
      if (dim[i] == 0) { // disconnected atoms
        B[i] = 0;
      }
      else if (dim[i] == 1) {
        B[i] = 2; // sp hybridizations and linear geometries
      }
      else if (dim[i] == 2) {
        if (atNumber < 11) { // sp2 hybridizations for second row elements
          B[i] = 3;
        }
        else {
          B[i] = 4; // square planar or sp3 hybridizations
        }
      }
      else if (dim[i] == 3) {
        if (atNumber < 11) { // square planar or sp3 hybridizations
          B[i] = 4;
        }
        else {
          B[i] = 7;
        }
      }
    }

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      std::vector<atom*> bondedAtoms = pAtom->getBondedAtoms();
      int numBonds = static_cast<int>(bondedAtoms.size());
      if (numBonds > B[i]) {
        for (int j = 0; j < B[i]-numBonds; j++) {
          // delete longest bond
          std::cout << " need to delete bond ... implement" << std::endl;
          //exit(0);
          throw MTKException(" need to delete bond ... implement");
        }
        Q[i] = B[i];
      }
      else {
        Q[i] = numBonds;
      }
    }

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      std::vector<int> b;
      std::vector<atom*> bondedAtoms = pAtom->getBondedAtoms();
      for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
        b.push_back(bondedAtoms[j]->getIndex()-1);
      }
      bdAtoms.push_back(b);
    }
    return 0;
}

// ============================================================
// Function : runLabute3()
// ------------------------------------------------------------
//
// ============================================================
int hybridize::runLabute3()
{
    errorMessage = "\n";

    //    3.1
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atNumbers[i] < 3) {
        hybridizations[i] = 4; // Hydrogen and Helium (sp3)
      }
    }

    //    3.2
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (Q[i] > 4 and atGroups[i] == 5) {
          hybridizations[i] = 5; // dsp3
        }
        if (Q[i] == 5 and (atGroups[i] > 3 and atGroups[i] < 9)) {
          hybridizations[i] = 5; // sp3d
        }
      }
    }

    //    3.3
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (Q[i] > 4 and atGroups[i] == 6) {
          hybridizations[i] = 6; // sp3d2
        }
        if (Q[i] == 6 and (atGroups[i] > 3 and atGroups[i] < 9)) {
          hybridizations[i] = 6; // sp3d2
        }
      }
    }

    //    3.4
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (Q[i] > 4 and atGroups[i] == 7) {
          hybridizations[i] = 7; // sp3d3
        }
        if (Q[i] == 7 and (atGroups[i] > 3 and atGroups[i] < 9)) {
          hybridizations[i] = 7; // sp3d3
        }
      }
    }

    //    3.5
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (Q[i] == 4 and atNumbers[i] > 10 and dim[i] == 2) {
          hybridizations[i] = 6; // sp3d2, square planar
        }
      }
    }

    //    3.6
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (atNumbers[i] > 21 and atNumbers[i] < 31) {
          hybridizations[i] = 6; // sp3d2, transition metals
        }
        else if (atNumbers[i] > 38 and atNumbers[i] < 49) {
          hybridizations[i] = 6; // sp3d2, transition metals
        }
        else if (atNumbers[i] > 70 and atNumbers[i] < 81) {
          hybridizations[i] = 6; // sp3d2, transition metals
        }
        else if (atNumbers[i] > 102 and atNumbers[i] < 113) {
          hybridizations[i] = 6; // sp3d2, transition metals
        }
      }
    }

    //    3.7
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (atNumbers[i] > 10 and atNumbers[i] != 14 and
            atNumbers[i] != 15 and atNumbers[i] != 16 and
            atNumbers[i] != 34) { // Si, P, S, Se
          if (Q[i] > 4) {
            hybridizations[i] = 6; // sp3d2, square planar
          }
          else {
            hybridizations[i] = 4; // sp3
          }
        }
      }
    }

    //    3.8
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (Q[i] == 4) {
          hybridizations[i] = 4; // sp3
        }
        else if (Q[i] == 3 and dim[i] == 3) {
          hybridizations[i] = 4; // sp3
        }
      }
    }

    //    3.9
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (Q[i] > 2 and (atGroups[i] > 5 and atGroups[i] < 9)) {
          hybridizations[i] = 4; // sp3
        }
      }
    }

    //    3.10
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (!hybridizations[i]) {
        if (atNumbers[i] != 6  and atNumbers[i] != 7  and atNumbers[i] != 8 and
            atNumbers[i] != 14 and atNumbers[i] != 15 and atNumbers[i] != 16 and
            atNumbers[i] != 34) { // C, N, O, Si, P, S, Se
          hybridizations[i] = 4; // sp3
        }
      }
    }

    //    3.11
    bool changed = true;
    unsigned int assigned = 0;
    while (changed) {
      changed = false;
      for (unsigned int i = 0; i < nAtoms; i++) {
       if (!hybridizations[i]) {
         assigned = 0;
         for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
           if (hybridizations[bdAtoms[i][j]]) assigned++;
         }
         if (assigned == bdAtoms[i].size()) {
           hybridizations[i] = 4; // sp3
           changed = true;
           break;
         }
       }
      }
    }

    errorMessage += " Atom # (I), Atom Symbol (S), # BONDED (K), DIMENSION (D),";
    errorMessage += " CONNECTIONS (Q),\n ATOMIC # (Z), GROUP (G), HYBRIDIZATION (H) \n";
    for (unsigned int i = 0; i < nAtoms; i++) {
     atom* pAt = atomList[i];
     errorMessage += " I=" + i2s(pAt->getFileID())
                   + " S=" + atSymbols[i]
                   + " K=" + i2s(bdAtoms[i].size())
                   + " D=" + i2s(dim[i])
                   + " Q=" + i2s(Q[i])
                   + " Z=" + i2s(atNumbers[i])
                   + " G=" + i2s(atGroups[i])
                   + " H=" + i2s(hybridizations[i]) + "\n";
    }
    if (!bonds.empty()) {

      errorMessage += " Bond #, Atom # (I) - Atom # (J), Bond Type \n";

      int i = 1;
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        errorMessage += " " + i2s(i) + " "
                      + i2s(pBond->atom1->getFileID()) + "-"
                      + i2s(pBond->atom2->getFileID()) + " "
                      + i2s(pBond->type) + "\n";
        i++;
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 3", errorMessage, 4);
    }

    return 0;
}

// ============================================================
// Function : runLabute()
// ------------------------------------------------------------
//
// ============================================================
int hybridize::runLabute()
{
    int failure = 0;
    failure = this->readParameters();
    if (failure) {
      return 1;
    }

    // STEP 1
    failure = this->runLabute1();
    if (failure) {
      errorMessage = " Labute step 1 ";
      errorLogger.throwError("hybridize::runLabute", errorMessage, 1);
      return 1;
    }

    // STEP 2
    failure = this->runLabute2();
    if (failure) {
      errorMessage = " Labute step 2 ";
      errorLogger.throwError("hybridize::runLabute", errorMessage, 1);
      return 1;
    }

    // STEP 3
    bonds = pParent->getBondMap();
    if (!bonds.empty()) {
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        pBond->type = 0; // initialize to undefined
      }
    }

    failure = this->runLabute3();
    if (failure) {
      errorMessage = " Labute step 3 ";
      errorLogger.throwError("hybridize::runLabute", errorMessage, 1);
      return 1;
    }

    // STEP 4
    errorMessage = "\n";
    if (!bonds.empty()) {
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        int at1 = pBond->atom1->getIndex()-1;
        int at2 = pBond->atom2->getIndex()-1;
        if (hybridizations[at1] or hybridizations[at2]) {
          pBond->type = 1;
          errorMessage += " Bond " + i2s(pBond->atom1->getFileID()) + "-" +
                          i2s(pBond->atom2->getFileID()) + " is single \n";
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 4", errorMessage, 4);
    }

    // STEP 5
    errorMessage = "\n";
    std::map<int, int> bondPlanar;
    if (!bonds.empty()) {
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        int bondIndex = b->first;
        bondPlanar[bondIndex] = 0;
        pBond = b->second;

        int atI = pBond->atom1->getIndex()-1;
        int atJ = pBond->atom2->getIndex()-1;

        if (dim[atI] == 1 and dim[atJ] == 1) { // linear atoms
          continue;
        }

        atom* atomI = atomList[atI];
        atom* atomJ = atomList[atJ];

        //double smallestOutOfPlaneDihedral = 0.0;
        double smallestOutOfPlaneDihedral = 180.0;
        int nTs = 0;

        for (unsigned int i = 0; i < bdAtoms[atI].size(); i++) {
          int atA = bdAtoms[atI][i];
          if (atA == atJ) continue;
          atom* atomA = atomList[atA];

          for (unsigned int j = 0; j < bdAtoms[atJ].size(); j++) {
            int atB = bdAtoms[atJ][j];
            if (atB == atI) continue;
            atom* atomB = atomList[atB];
            if (atSymbols[atA] == "H" or atSymbols[atB] == "H") continue;
            nTs++;
            double dihedralAIJB = torsion(*(atomA->getCoords()), *(atomI->getCoords()),
                                          *(atomJ->getCoords()), *(atomB->getCoords()));
            double m = std::min(std::abs(dihedralAIJB),
                                std::min(std::abs(PI-dihedralAIJB),
                                         std::abs(-PI-dihedralAIJB)));

//            if (m > smallestOutOfPlaneDihedral) {
            if (m < smallestOutOfPlaneDihedral) {
              smallestOutOfPlaneDihedral = m;
            }
            /*
            errorMessage += "  " + i2s(atomA->getFileID()) + "-" + i2s(atomI->getFileID()) + "-" +
                                   i2s(atomJ->getFileID()) + "-" +  i2s(atomB->getFileID()) + " " +
                                   d2s(dihedralAIJB * RAD2DEG) + " " + d2s(m * RAD2DEG) + " " + d2s(smallestOutOfPlaneDihedral * RAD2DEG) + "\n";
            */
          }
        }

        if (nTs > 0) {
          // Set up a flag for this step -- it doesn't work for Phenyl-N=N-Phenyl systems
          if (smallestOutOfPlaneDihedral*RAD2DEG > this->smallOutOfPlaneDihedral) { // used to be 15
            errorMessage += "Setting bond " + i2s(pBond->atom1->getFileID()) + "-" +
                            i2s(pBond->atom2->getFileID()) + " to single.";
            errorMessage += " Smallest out of plane dihedral = " +
                            d2s(smallestOutOfPlaneDihedral*RAD2DEG) + "\n";

            pBond->type = 1;
          }
          else if (smallestOutOfPlaneDihedral*RAD2DEG < 10) {
            bondPlanar[bondIndex] = 1;
            errorMessage += "Setting bond " + i2s(pBond->atom1->getFileID()) + "-" +
                            i2s(pBond->atom2->getFileID()) + " to planar.";
            errorMessage += " Smallest out of plane dihedral = " +
                            d2s(smallestOutOfPlaneDihedral*RAD2DEG) + "\n";
          }
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 5", errorMessage, 4);
    }

    // STEP 5b
    errorMessage = "\n";
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      if (bdAtoms[i].size() == 3) {
        atom* pAt2 = atomList[bdAtoms[i][0]];
        atom* pAt3 = atomList[bdAtoms[i][1]];
        atom* pAt4 = atomList[bdAtoms[i][2]];
        double improperABCD = std::abs(torsion(*(pAt2->getCoords()),
                                               *(pAt1->getCoords()),
                                               *(pAt3->getCoords()),
                                               *(pAt4->getCoords()))) * RAD2DEG;

        if (improperABCD > 170.0) {
          planarFlags[i] = 1;
          errorMessage += "Setting atom " + i2s(pAt1->getFileID()) +
                          " to planar." +
                          " Improper dihedral = " +
                          d2s(improperABCD) + "\n";
        }
      }
    }
    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 5B", errorMessage, 4);
    }

    // STEP 6
    errorMessage = "\n";
    double sBondDist = 0.0;
    if (!bonds.empty()) {
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        int atI = pBond->atom1->getIndex()-1;
        atom* atomI = atomList[atI];

        int atJ = pBond->atom2->getIndex()-1;
        atom* atomJ = atomList[atJ];

        std::string bT1 = atSymbols[atI] + "-" + atSymbols[atJ];
        std::string bT2 = atSymbols[atJ] + "-" + atSymbols[atI];

        mapIterator sB1 = singleBonds.find(bT1);
        mapIterator sB2 = singleBonds.find(bT2);

        sBondDist = 0.0;
        if (sB1 != singleBonds.end()) {
          sBondDist = singleBonds[bT1];
        }
        else if (sB2 != singleBonds.end()) {
          sBondDist = singleBonds[bT2];
        }

        if (sBondDist > 0.0) {
          double distIJ = atomI->getCoords()->dist(*atomJ->getCoords());
          if (distIJ > sBondDist-0.05) {
            errorMessage += "Setting bond " + i2s(atomI->getFileID()) + "-" +
                         i2s(atomJ->getFileID()) + " to single." +
                         " Dist = " +
                         d2s(distIJ) + ", Std Dist = " + d2s(sBondDist-0.05) + "\n";

            pBond->type = 1;
          }
        }
        if (bT1 == "S-N" or bT1 == "N-S") {
          pBond->type = 1;

          errorMessage += "Setting bond (S-N) " + i2s(atomI->getFileID()) + "-" +
                          i2s(atomJ->getFileID()) + " to single. \n";

        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 6", errorMessage, 4);
    }

    // STEP 7
    errorMessage = "\n";
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      if (!hybridizations[i]) {
        unsigned int assigned = 0;
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
           atom* pAt2 = atomList[bdAtoms[i][j]];
           pBond = pParent->getBond(pAt1, pAt2);
           if (pBond) {
             if (pBond->type) assigned++;
           }
        }
        //if (assigned == bdAtoms[i].size()) {
        //  hybridizations[i] = 4; // sp3
        //  errorMessage += "Setting atom " + i2s(pAt1->getFileID()) +
        //                  " to sp3. \n";
        //}
        if (assigned) {
          hybridizations[i] = 4; // sp3
          errorMessage += "Setting atom " + i2s(pAt1->getFileID()) +
                          " to sp3. \n";
        }
      }
    }

    // Note: At this point a bond between two sp3 atoms can be left
    //       unassigned and so could be changed to double e.g. 1D7X

    // R-O-CH3
    //
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (hybridizations[i]) continue;
      atom* pAt1 = atomList[i];
      if ((pAt1->getElement()->symbol != "O") and
          (pAt1->getElement()->symbol != "S")) continue;
      if (bdAtoms[i].size() != 2) continue;

      atom* pAtA = atomList[bdAtoms[i][0]];
      atom* pAtB = atomList[bdAtoms[i][1]];
      if (bdAtoms[bdAtoms[i][0]].size() == 1) { // terminal
        pBond = pParent->getBond(pAt1, pAtA);
        if (pBond) {
          hybridizations[i] = 4; // sp3
          hybridizations[bdAtoms[i][0]] = 4;
          pBond->type = 1;
        }

        pBond = pParent->getBond(pAt1, pAtB);
        if (pBond) {
          pBond->type = 1;
        }
      }
      else {
        pBond = pParent->getBond(pAt1, pAtB);
        if (pBond) {
          hybridizations[i] = 4; // sp3
          hybridizations[bdAtoms[i][1]] = 4;
          pBond->type = 1;
        }

        pBond = pParent->getBond(pAt1, pAtA);
        if (pBond) {
          pBond->type = 1;
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 7", errorMessage, 4);
    }

    // STEP 8
    errorMessage = "\n";
    graph* molGraph = new graph();
    vertex* pVertexI = 0;
    vertex* pVertexJ = 0;
    std::vector<int> piAtoms;

    if (!bonds.empty()) {
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        int atI = pBond->atom1->getIndex();
        int atJ = pBond->atom2->getIndex();
        pVertexI = 0;
        pVertexJ = 0;

        if (!pBond->type) {
          std::stringstream stAtI;
          stAtI << atI;
          std::stringstream stAtJ;
          stAtJ << atJ;

          pVertexI = molGraph->getVertex(atI);
          pVertexJ = molGraph->getVertex(atJ);

          if (!pVertexI) {
            pVertexI = molGraph->addVertex(atI);
            pVertexI->setName(stAtI.str().c_str());
            piAtoms.push_back(atI);
          }
          if (!pVertexJ) {
            pVertexJ = molGraph->addVertex(atJ);
            pVertexJ->setName(stAtJ.str().c_str());
            piAtoms.push_back(atJ);
          }
          if (pVertexI and pVertexJ) {
            molGraph->addEdge(pVertexI, pVertexJ);
          }
        }
      }
    }

    // Find subgraphs
    std::vector<graph*> subGraphs;
    std::vector<vertex*> vertices = molGraph->getVertices();
    std::vector<vertex*> blockVertices;

    for (unsigned int i = 0; i < vertices.size(); i++) {
      if (!vertices[i]->isVisited()) {
        graph* subGraph = new graph();

        molGraph->dfs(vertices[i]);

        for (unsigned int j = 0; j < vertices.size(); j++) {
          if (vertices[j]->isVisited()) {
            bool taken = false;
            for (unsigned int k = 0; k < subGraphs.size(); k++) {
              if (subGraphs[k]->getVertex(vertices[j]->getIndex())) {
                taken = true;
              }
            }
            if (!taken) {
              subGraph->addVertex(vertices[j]);
            }
          }
        }
        subGraphs.push_back(subGraph);
      }
    }
    molGraph->reset();

    // Update edges in new subgraphs
    for (unsigned int k = 0; k < subGraphs.size(); k++) {
      std::vector<vertex*> subGV = subGraphs[k]->getVertices();
      for (unsigned int l = 0; l < subGV.size(); l++) {
        for (unsigned int p = l; p < subGV.size(); p++) {
          if (molGraph->hasEdge(subGV[l]->getIndex(), subGV[p]->getIndex())) {
            subGraphs[k]->addEdge(subGV[l], subGV[p]);
          }
        }
      }
    }

    typedef std::map<int, edge*>::iterator edgeIterator;
    typedef std::vector<vertex*>::iterator vertexIterator;

    std::vector<double> atomWeights; // weights for each atom
    atomWeights.assign(nAtoms, 0.0);

    // assign weights to all atoms
    for (unsigned int i = 0; i < subGraphs.size(); i++) {
      std::vector<vertex*> subGV = subGraphs[i]->getVertices();
      for (unsigned int j = 0; j < subGV.size(); j++) {
        int atIJ = subGV[j]->getIndex()-1;
        switch (atGroups[atIJ]) {
          case 14 : { // Carbon and Silicon
            bool bO = false;
            bool bN = false;
            for (unsigned int k = 0; k < bdAtoms[atIJ].size(); k++) {
              int atIJK = bdAtoms[atIJ][k];
              if (atSymbols[atIJK] == "O") bO = true;
              if (atSymbols[atIJK] == "N") bN = true;
            }

            if (bO) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["C-O"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["C-O"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["C-O"][2];
              }
            }
            else if (bN) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["C-N"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["C-N"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["C-N"][2];
              }
            }
            else {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["C"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["C"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["C"][2];
              }
            }
            if (atPeriods[atIJ] > 2) {
              atomWeights[atIJ] -= 0.1;
            }
            break;
          }
          case 15 : { // Nitrogen and Phosphorus
            bool bO = false;
            bool bN = false;
            bool bS = false;
            for (unsigned int k = 0; k < bdAtoms[atIJ].size(); k++) {
              int atIJK = bdAtoms[atIJ][k];
              if (atSymbols[atIJK] == "C") {
                for (unsigned int l = 0; l < bdAtoms[atIJK].size(); l++) {
                  int atIJKL = bdAtoms[atIJK][l];
                  if (atIJKL != atIJ) {
                    if (atSymbols[atIJKL] == "O") bO = true;
                    if (atSymbols[atIJKL] == "N") bN = true;
                  }
                }
              }
              else if (atSymbols[atIJK] == "S") {
                bS = true;
              }
            }
            if (bO) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["N-C-O"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["N-C-O"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["N-C-O"][2];
              }
            }
            else if (bN) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["N-C-N"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["N-C-N"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["N-C-N"][2];
              }
            }
            else {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["N"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["N"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["N"][2];
                if (bS) {
                  atomWeights[atIJ] = -20.0;
                }
              }
            }
            if (atPeriods[atIJ] > 2) {
              atomWeights[atIJ] -= 0.1;
            }
            break;
          }
          case 16 : { // Oxygen, Sulfur, and Selenium
            bool bO = false;
            bool bN = false;
            bool bS = false;
            for (unsigned int k = 0; k < bdAtoms[atIJ].size(); k++) {
              int atIJK = bdAtoms[atIJ][k];
              if (atSymbols[atIJK] == "C") {
                for (unsigned int l = 0; l < bdAtoms[atIJK].size(); l++) {
                  int atIJKL = bdAtoms[atIJK][l];
                  if (atIJKL != atIJ) {
                    if (atSymbols[atIJKL] == "O") bO = true;
                    if (atSymbols[atIJKL] == "N") bN = true;
                    if (atSymbols[atIJKL] == "S") bO = true;
                  }
                }
              }
            }
            if (bO) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-O"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-O"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-O"][2];
              }
            }
            else if (bN) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-N"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-N"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-N"][2];
              }
            }
            else if (bS) {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-O"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-O"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["O-C-O"][2];
              }
            }
            else {
              if (Q[atIJ] == 1) {
                atomWeights[atIJ] = logLikelihoodRatios["O"][0];
              }
              else if (Q[atIJ] == 2) {
                atomWeights[atIJ] = logLikelihoodRatios["O"][1];
              }
              else if (Q[atIJ] == 3) {
                atomWeights[atIJ] = logLikelihoodRatios["O"][2];
              }
            }
            if (atPeriods[atIJ] > 2) {
              atomWeights[atIJ] -= 0.1;
            }
            break;
          }
          default : {
            atomWeights[atIJ] = -20.0;
            atom* pAtTemp = pParent->getAtom(atIJ+1, 1, 0);
            errorMessage += " Using default weight for atom " + i2s(pAtTemp->getFileID()) + " ";
          }
        }
      }

      // compute the edge weights for the Maximum Weighted Matching Algorithm
      std::map<int, edge*>  subGE = subGraphs[i]->getEdges();
      if (!subGE.empty()) {
        for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
          edge* pEdge = e->second;
          int atI = pEdge->v1->getIndex()-1;
          int atJ = pEdge->v2->getIndex()-1;
          atom* atomI = atomList[atI];
          atom* atomJ = atomList[atJ];
          std::string bT1 = atSymbols[atI]+"-"+atSymbols[atJ];
          std::string bT2 = atSymbols[atJ]+"-"+atSymbols[atI];
          mapIterator sB1 = singleBonds.find(bT1);
          mapIterator sB2 = singleBonds.find(bT2);

          sBondDist = 0.0;
          if (sB1 != singleBonds.end()) {
            sBondDist = singleBonds[bT1];
          }
          else if (sB2 != singleBonds.end()) {
            sBondDist = singleBonds[bT2];
          }

          double doubleBond = 0.0;
          if (sBondDist > 0.0) {
            double distIJ = atomI->getCoords()->dist(*atomJ->getCoords());
            if (distIJ < sBondDist - edgeWeightSingle) {
              doubleBond += 2 * 1.0;
            }
            if (distIJ < sBondDist - edgeWeightDouble) {
              doubleBond += 1.0;
            }
          }

          double Wij = atomWeights[atI] + atomWeights[atJ] + doubleBond;
          pEdge->itsValue = Wij;
        }
      }
    }

    errorMessage = "\n";
    for (unsigned int k = 0; k < subGraphs.size(); k++) {
      std::vector<vertex*> subGV = subGraphs[k]->getVertices();
      std::map<int, edge*>  subGE = subGraphs[k]->getEdges();

      errorMessage += "Subgraph: " + i2s(k+1) + " \n";
      for (unsigned int l = 0; l < subGV.size(); l++) {
        atom* pAtTemp = pParent->getAtom(subGV[l]->getIndex(), 1, 0);
        errorMessage += i2s(pAtTemp->getFileID()) + " (" + d2s(atomWeights[subGV[l]->getIndex()-1]) + ") ";
      }
      errorMessage += "\n";

      if (!subGE.empty()) {
        for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
          edge* pEdge = e->second;
          atom* pAtTemp1 = pParent->getAtom(pEdge->v1->getIndex(), 1, 0);
          atom* pAtTemp2 = pParent->getAtom(pEdge->v2->getIndex(), 1, 0);
          errorMessage += i2s(pAtTemp1->getFileID()) + "-" + i2s(pAtTemp2->getFileID()) +
                          " (" + d2s(pEdge->itsValue) + ")"+ "; ";
        }
      }
      errorMessage += "\n";
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 8", errorMessage, 4);
    }

    // Maximum Weighted Matching Algorithm [Greedy Algorithm]
    bool finished = false;
    std::vector<graph*> bestMatches;
    for (unsigned int i = 0; i < subGraphs.size(); i++) {
      finished = false;
      graph* bestMatch = 0;
      double matchScore = 0.0;
      double bestMatchScore = 0.0;
      double bestEdgeWeight = 0.0;
      std::map<int, edge*>  subGE = subGraphs[i]->getEdges();
      if (!subGE.empty()) {
        for (edgeIterator e1 = subGE.begin(); e1 != subGE.end(); e1++) {
          edge* pEdge1 = e1->second;
          finished = false;
          graph* Match = new graph();
          edge* bestEdge = 0;
          bestEdgeWeight = 0.0;
          matchScore = 0.0;
          matchScore += pEdge1->itsValue;

          vertex* v11 = Match->addVertex(pEdge1->v1);
          vertex* v12 = Match->addVertex(pEdge1->v2);
          Match->addEdge(v11, v12);
          while (!finished) {
            finished = true;
            bestEdgeWeight = -100.0;

            for (edgeIterator e2 = subGE.begin(); e2 != subGE.end(); e2++) {
              edge* pEdge2 = e2->second;
              vertex* v21 = pEdge2->v1;
              vertex* v22 = pEdge2->v2;

              bool testing = false; ////
              std::vector<vertex*> matchEdges = Match->getVertices();
              for (vertexIterator vt = matchEdges.begin(); vt != matchEdges.end(); vt++) {
                vertex* verTex = *vt;
                if (subGraphs[i]->hasEdge(verTex->getIndex(), v21->getIndex()) or
                    subGraphs[i]->hasEdge(verTex->getIndex(), v22->getIndex())) {
                  testing = true;
                }
              }
              if (!testing) continue; /////

              unsigned int different = 0;
              if (!Match->hasEdge(v21->getIndex(), v22->getIndex())) {
                std::map<int, edge*> matchEdges = Match->getEdges();
                for (edgeIterator e3 = matchEdges.begin(); e3 != matchEdges.end(); e3++) {
                  edge* pEdge3 = e3->second;
                  vertex* v31 = pEdge3->v1;
                  vertex* v32 = pEdge3->v2;
                  if ((v21->getIndex() != v31->getIndex() and v21->getIndex() != v32->getIndex() and
                       v22->getIndex() != v31->getIndex() and v22->getIndex() != v32->getIndex())) {
                    different++;
                  }
                }
                if (different == matchEdges.size()) {
                  if (pEdge2->itsValue > bestEdgeWeight) {
                    bestEdgeWeight = pEdge2->itsValue;
                    bestEdge = pEdge2;
                    finished = false;
                  }
                }
              }
            }
            if (!finished and bestEdge) {
              vertex* vM1 = Match->addVertex(bestEdge->v1);
              vertex* vM2 = Match->addVertex(bestEdge->v2);
              Match->addEdge(vM1, vM2);
              matchScore += bestEdge->itsValue;
            }
          }
          if (matchScore > bestMatchScore) {
            bestMatchScore = matchScore;
            delete bestMatch;
            bestMatch = Match;
          }
          else if (bestMatchScore - matchScore < 0.1) {
            // HAS A RING BEEN FORMED?
            unsigned int nBs = 0;
            std::vector<vertex*> subMatchGV = Match->getVertices();
            std::map<int, edge*> subMatchGE = Match->getEdges();
            for (unsigned int j = 0; j < subMatchGV.size(); j++) {
              int atI = subMatchGV[j]->getIndex()-1;
              atom* atomI = atomList[atI];
              for (unsigned int j2 = j+1; j2 < subMatchGV.size(); j2++) {
                int atJ = subMatchGV[j2]->getIndex()-1;
                atom* atomJ = atomList[atJ];
                pBond = pParent->getBond(atomI, atomJ);
                if (pBond) {
                  nBs++;
                }
              }
            }
            if (nBs == 2*subMatchGE.size()) {
              bestMatchScore = matchScore;
              delete bestMatch;
              bestMatch = Match;
            }
            else {
              delete Match;
              Match = 0;
            }
          }
          else {
            delete Match;
            Match = 0;
          }
        }
      }
      if (bestMatch) {
        bestMatches.push_back(bestMatch);
      }
      //else {
      //  std::cout << "    best match is not found " << std::endl;
      //}
    }
    errorLogger.flush();

    errorMessage = "\n";
    for (unsigned int i = 0; i < bestMatches.size(); i++) {
      std::vector<vertex*> subGV = bestMatches[i]->getVertices();
      std::map<int, edge*>  subGE = bestMatches[i]->getEdges();
      errorMessage += "Match: " + i2s(i+1) + " \n";

      for (unsigned int j = 0; j < subGV.size(); j++) {
        atom* pAtTemp = pParent->getAtom(subGV[j]->getIndex(), 1, 0);
        errorMessage += i2s(pAtTemp->getFileID()) + " ";
      }
      errorMessage += "\n";
      if (!subGE.empty()) {
        for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
          edge* pEdge = e->second;

          atom* pAtTemp1 = pParent->getAtom(pEdge->v1->getIndex(), 1, 0);
          atom* pAtTemp2 = pParent->getAtom(pEdge->v2->getIndex(), 1, 0);
          errorMessage += i2s(pAtTemp1->getFileID()) + "-" + i2s(pAtTemp2->getFileID()) + "; ";
        }
      }
      errorMessage += "\n";
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 8 (Maximum Weighted Matching Algorithm)",
      errorMessage, 4);
    }

    for (unsigned int i = 0; i < subGraphs.size(); i++) {
      std::vector<vertex*> subGV = subGraphs[i]->getVertices();
      std::map<int, edge*>  subGE = subGraphs[i]->getEdges();

      if (i+1 > bestMatches.size()) {
        for (unsigned int l = 0; l < subGV.size(); l++) {
          if (planarFlags[subGV[l]->getIndex()-1]) {
            hybridizations[subGV[l]->getIndex()-1] = 3; // sp2
          }
          else {
            hybridizations[subGV[l]->getIndex()-1] = 4; // sp3
          }
        }
        if (!subGE.empty()) {
          for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
            edge* pEdge = e->second;
            int atI = pEdge->v1->getIndex()-1;
            int atJ = pEdge->v2->getIndex()-1;
            atom* atomI = atomList[atI];
            atom* atomJ = atomList[atJ];
            pBond = pParent->getBond(atomI, atomJ);
            pBond->type = 1; // single bond
          }
        }
      }
      else {
        for (unsigned int l = 0; l < subGV.size(); l++) {
          if (!bestMatches[i]->getVertex(subGV[l]->getIndex())) {
            if (planarFlags[subGV[l]->getIndex()-1]) {
              hybridizations[subGV[l]->getIndex()-1] = 3; // sp2
            }
            else {
              hybridizations[subGV[l]->getIndex()-1] = 4; // sp3
            }
          }
         }

        if (!subGE.empty()) {
          for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
            edge* pEdge = e->second;
            if (!bestMatches[i]->hasEdge(pEdge->v1->getIndex(), pEdge->v2->getIndex())) {
              int atI = pEdge->v1->getIndex()-1;
              int atJ = pEdge->v2->getIndex()-1;
              atom* atomI = atomList[atI];
              atom* atomJ = atomList[atJ];
              pBond = pParent->getBond(atomI, atomJ);
              pBond->type = 1; // single bond
            }
          }
        }
      }
    }

    errorMessage = "\n";
    for (unsigned int i = 0; i < bestMatches.size(); i++) {
      std::map<int, edge*>  subGE = bestMatches[i]->getEdges();
      if (!subGE.empty()) {
        for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
          edge* pEdge = e->second;
          int atI = pEdge->v1->getIndex()-1;
          int atJ = pEdge->v2->getIndex()-1;
          atom* atomI = atomList[atI];
          atom* atomJ = atomList[atJ];
          pBond = pParent->getBond(atomI, atomJ);

          std::string bT1 = atSymbols[atI] + "-" + atSymbols[atJ];
          std::string bT2 = atSymbols[atJ] + "-" + atSymbols[atI];

          mapIterator sB1 = singleBonds.find(bT1);
          mapIterator sB2 = singleBonds.find(bT2);

          sBondDist = 0.0;
          if (sB1 != singleBonds.end()) {
            sBondDist = singleBonds[bT1];
          }
          else if (sB2 != singleBonds.end()) {
            sBondDist = singleBonds[bT2];
          }
          else {
            //std::cout << " CAN'T FIND BOND " << atI+1 << "-" << atJ+1 << std::endl;
            //exit(0);
            std::stringstream ss;
            ss << " CAN'T FIND BOND " << atI+1 << "-" << atJ+1 << std::endl;
            std::cout << ss.str();
            throw MTKException(ss.str());
          }

          if (dim[atI] == 1 and dim[atJ] == 1) { // linear atoms
            double distIJ = atomI->getCoords()->dist(*atomJ->getCoords());
            if (distIJ < sBondDist - tripleBondParameter) {
              hybridizations[atI] = 2; // sp
              hybridizations[atJ] = 2; // sp
              if (pBond) {
                errorMessage += "Setting bond " + i2s(atomI->getFileID()) + "-" +
                                                  i2s(atomJ->getFileID()) + " to triple. \n";
                pBond->type = 3;
              }
            }
            else {
              errorMessage += "ERROR: NEED TO IMPLEMENT THIS \n";
            }
          }
          else {
            hybridizations[atI] = 3; // sp2
            hybridizations[atJ] = 3; // sp2
            if (pBond) {
              pBond->type = 2;
            }
            else {
              errorMessage += "Can't find bond " + i2s(atomI->getFileID()) + "-" +
                                                    i2s(atomJ->getFileID()) + ". \n";
            }
          }
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 8 ", errorMessage, 4);
    }

    // STEP 9
    errorMessage = "\n";
    errorMessage += "The formal charge of atom i, f_i, is calculated as follows:\n";
    errorMessage += "  f_i = c_i - o_i + b_i\n";
    errorMessage += "   where:\n";
    errorMessage += "    c_i = atom group in periodic table\n";
    errorMessage += "    o_i = nominal octet (2 for hydrogen, 6 for boron, 8 for carbon and \n";
    errorMessage += "          all other sp3 atoms in groups 5,6,7,8)\n";
    errorMessage += "    b_i = sum of the atom bond orders\n";

    formalCharges.assign(nAtoms, 0);
    int totalCharge = 0;

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      formalCharges[i] = 0;

      // c_i
      int c_i = 0;
      if (atGroups[i] > 12) {
        c_i = atGroups[i] - 10;
      }
      else {
        c_i = atGroups[i];
      }

      // o_i
      int o_i = 0;
      if (atNumbers[i] == 1) { // hydrogen
        o_i = 2;
      }
      else if (atNumbers[i] == 5) { // boron
        o_i = 6;
      }
      else if (atGroups[i] > 13) { // C, N, O groups
        o_i = 8;
      }
      else if (atGroups[i] > 2 and atGroups[i] < 13) { // transition metals
        if (Q[i] == 0) {
          o_i = 0; // isolated atom
        }
        else {
          o_i = atGroups[i]; /// PROBABLY WRONG
        }
      }

      // b_i
      int b_i = 0;
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        atom* pAt2 = atomList[bdAtoms[i][j]];
        pBond = pParent->getBond(pAt1, pAt2);
        if (pBond) {
          b_i += pBond->type;
        }
      }
      formalCharges[i] = c_i - o_i + b_i;

      totalCharge+=formalCharges[i];
      //std::cout << i << ":\t" << c_i << " " << " " << o_i << " " << b_i << " " << formalCharges[i] << std::endl;

      errorMessage += " " + i2s(pAt1->getFileID()) + " FC_i = " +
                      i2s(formalCharges[i]) + " c_i = " + i2s(c_i) + " o_i = "
                      + i2s(o_i) + " b_i = " + i2s(b_i) + "\n";
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Steps 9.0 to 9.4 ", errorMessage, 4);
    }

    //    STEP 9.5
    errorMessage = "\n";
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (formalCharges[i] < 0) {
        atom* pAt1 = atomList[i];
        bool positiveNeighbor = false;
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (formalCharges[bdAtoms[i][j]] > 0) {
            positiveNeighbor = true;
          }
        }
        if (!positiveNeighbor) {
          formalCharges[i] = 0;
          errorMessage += " Setting Formal Charge of " + i2s(pAt1->getFileID()) + " to 0 \n";
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 9.5 ", errorMessage, 4);
    }

    //   STEP 9.6
    errorMessage = "\n";
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (formalCharges[i] > 0) {
        atom* pAt1 = atomList[i];
        errorMessage += " Formal Charge of " + i2s(pAt1->getFileID()) + " is " + i2s(formalCharges[i]) + "\n";
        int changeI = 0;
        for (int j = 0; j < formalCharges[i]; j++) {
          // get most electronegative bonded atom
          int mostElectronegativeBdAtom = 0;
          double mostElectronegative = 0.0;
          int mostNegFC = 5;
/*
          for (unsigned int k = 0; k < bdAtoms[i].size(); k++) {
            atom* pBdAt = atomList[bdAtoms[i][k]];
            errorMessage += "  Bonded to " + i2s(pBdAt->getFileID()) +
                            " EN = " + d2s(atENs[bdAtoms[i][k]]) +
                            " FC = " + i2s(formalCharges[bdAtoms[i][k]]);

            if (atENs[bdAtoms[i][k]] >= mostElectronegative) {
              mostElectronegative = atENs[bdAtoms[i][k]];
              if (formalCharges[bdAtoms[i][k]] <= mostNegFC) {
                mostElectronegativeBdAtom = bdAtoms[i][k];
                mostNegFC = formalCharges[bdAtoms[i][k]];
                errorMessage += " mostElectronegativeBdAtom ";
              }
            }
            errorMessage += "\n";
          }
*/

          for (unsigned int k = 0; k < bdAtoms[i].size(); k++) {
            atom* pBdAt = atomList[bdAtoms[i][k]];
            errorMessage += "  Bonded to " + i2s(pBdAt->getFileID()) +
                            " EN = " + d2s(atENs[bdAtoms[i][k]]) +
                            " FC = " + i2s(formalCharges[bdAtoms[i][k]]) + "\n";
            //if (atENs[bdAtoms[i][k]] > mostElectronegative) {
            if (atENs[bdAtoms[i][k]] >= mostElectronegative) {
              mostElectronegative = atENs[bdAtoms[i][k]];
/*
            }
          }

          for (unsigned int k = 0; k < bdAtoms[i].size(); k++) {
            atom* pBdAt = atomList[bdAtoms[i][k]];
            if (atENs[bdAtoms[i][k]] == mostElectronegative) {
*/
              if (formalCharges[bdAtoms[i][k]] <= mostNegFC) {
                mostElectronegativeBdAtom = bdAtoms[i][k];
                mostNegFC = formalCharges[bdAtoms[i][k]];
                errorMessage += "    mostElectronegativeBdAtom " + i2s(pBdAt->getFileID());
              }
            }
          }

          atom* pAt2 = atomList[mostElectronegativeBdAtom];
          pBond = pParent->getBond(pAt1, pAt2);

          // This works for S(=O)(=O)(NH2)(C)
          if (pBond) {
            changeI++;
            formalCharges[mostElectronegativeBdAtom]++;
            pBond->type++;

            errorMessage += " Changing Bond " + i2s(pAt1->getFileID()) + "-" +
                            i2s(pAt2->getFileID()) + " to " +
                            i2s(pBond->type) +" \n";
          }
        }
        if (changeI) formalCharges[i] -= changeI;

        // Change least electronegative to zero
        std::vector<int> updated;
        typedef std::vector<int>::iterator intIterator;

        int nLEN = bdAtoms[i].size() - changeI;
        for (int j = 0; j < nLEN; j++) {
          int leastElectronegativeBdAtom = 0;
          double leastElectronegative = 0.0;

          for (unsigned int k = 0; k < bdAtoms[i].size(); k++) {
            intIterator u = std::find(updated.begin(), updated.end(), bdAtoms[i][k]);
            if (atENs[bdAtoms[i][k]] <= leastElectronegative and u == updated.end()) {
              updated.push_back(bdAtoms[i][k]);
              leastElectronegative = atENs[bdAtoms[i][k]];
              leastElectronegativeBdAtom = bdAtoms[i][k];
            }
          }
          formalCharges[leastElectronegativeBdAtom] = 0;
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 9.6 ", errorMessage, 4);
    }

    errorMessage = "\n";
    // Make sure amide nitrogens are sp2
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "N") {
        atom* pAt1 = atomList[i];
        // amide
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (atSymbols[bdAtoms[i][j]] == "C") {
            int cAtom = bdAtoms[i][j];
            atom* pAt2 = atomList[cAtom];
            Bond* pBdNC = pParent->getBond(pAt1, pAt2);
            if (pBdNC->type == 1) {
              for (unsigned int k = 0; k < bdAtoms[cAtom].size(); k++) {
                if (atSymbols[bdAtoms[cAtom][k]] == "O") {
                  int oAtom = bdAtoms[cAtom][k];
                  atom* pAt3 = atomList[oAtom];
                  Bond* pBdCO = pParent->getBond(pAt2, pAt3);
                  if (pBdCO->type == 2) {
                    hybridizations[i] = 3; // sp2

                    errorMessage += " Amide Bond (N-C=O): " + i2s(pAt1->getFileID()) + "-" +
                                    i2s(pAt2->getFileID()) + "-" + i2s(pAt3->getFileID()) +
                                    ". Setting N to sp2  \n";
                  }
                }
              }
            }
          }
        }
      }
    }

    /*
             O-
           /
       R-N+
          \\
           O
    */
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "N") {
        atom* pAt1 = atomList[i];
        int Oxygen = 0;
        atom* pO1 = 0;
        Bond* pNOBond = 0;
        int nODoubleBonds = 0;
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (atSymbols[bdAtoms[i][j]] == "O") {
            int cAtom = bdAtoms[i][j];
            atom* pAt2 = atomList[cAtom];
            Bond* pBd = pParent->getBond(pAt1, pAt2);
            if (pBd->type == 2) {
              pO1 = pAt2;
              Oxygen = bdAtoms[i][j];
              pNOBond = pBd;
              nODoubleBonds++;
            }
          }
        }
        if (nODoubleBonds == 2) {
          pNOBond->type = 1;
          formalCharges[i] = 1;
          formalCharges[Oxygen] = -1;

          errorMessage += " R-N+(=O)-O- Found: (N) " + i2s(pAt1->getFileID()) + " \n";
        }
      }
    }

    // make sure one bond in PO4 is double
    double smallestBond = 100.0;
    Bond* smallestPOBond = 0;
    int nPOBonds = 0;
    for (unsigned int i = 0; i < nAtoms; i++) {
      smallestBond = 100.0;
      smallestPOBond = 0;
      if (atSymbols[i] == "P") {
        atom* pAt1 = atomList[i];
        bool doublePOFound = false;
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (atSymbols[bdAtoms[i][j]] == "O") {
            nPOBonds++;
            int cAtom = bdAtoms[i][j];
            atom* pAt2 = atomList[cAtom];
            Bond* pBdPO = pParent->getBond(pAt1, pAt2);

            errorMessage += " P (" + i2s(pAt1->getFileID()) + ") - O ("
                          + i2s(pAt2->getFileID()) + ") type = "
                          + i2s(pBdPO->type) + "\n";


            if (pBdPO->type == 1) {
              double distPO = pAt1->getCoords()->dist(*pAt2->getCoords());
              if (distPO < smallestBond) {
                smallestBond = distPO;
                smallestPOBond = pBdPO;
              }
            }
            else if (pBdPO->type == 2) {
              doublePOFound = true;
            }
          }
        }
        if (nPOBonds >= 2 and !doublePOFound) {
          smallestPOBond->type = 2;
          errorMessage += " R-PO" +i2s(nPOBonds) + " Found: (P) " + i2s(pAt1->getFileID()) + " \n";
          errorMessage += "  Setting " + i2s(smallestPOBond->atom1->getFileID()) + "-" + i2s(smallestPOBond->atom2->getFileID())
                        + " to double \n"; 
        }
      }
    }

    //    The following step is not in the published algorithm,
    //      required to perceive sulfoxide or sulfinyl functional groups
    //    STEP 9.7 .. set all atoms (fi < 0) without a positive neighbor to 0
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (formalCharges[i] < 0) {
        bool positiveNeighbor = false;
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (formalCharges[bdAtoms[i][j]] > 0) {
            positiveNeighbor = true;
          }
        }
        if (!positiveNeighbor) {
          formalCharges[i] = 0;
          atom* pAt1 = atomList[i];
          errorMessage += " Setting Charge of " + i2s(pAt1->getFileID()) + " to 0 \n";
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 9.7 ", errorMessage, 4);
    }

    // STEP 9.8
    errorMessage = "\n";
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (planarFlags[i]) {
        atom* pAt1 = atomList[i];
        bool bOK = false;
        if (bdAtoms[i].size() == 3) {
          for (unsigned int j = 0; j < 3; j++) {
            atom* pAt2 = atomList[bdAtoms[i][j]];
            pBond = pParent->getBond(pAt1, pAt2);
            if (pBond->type == 2) {
              bOK = true;
            }
          }

          int mostLikelyAtom = -1;
          double mostLikelyAtomDist = 10.0;
          if (!bOK) {
            for (unsigned int j = 0; j < 3; j++) {
              int a2 = bdAtoms[i][j];
              atom* pAt2 = atomList[a2];

              std::string bT1 = atSymbols[i] + "-" + atSymbols[a2];
              std::string bT2 = atSymbols[a2] + "-" + atSymbols[i];

              mapIterator sB1 = singleBonds.find(bT1);
              mapIterator sB2 = singleBonds.find(bT2);

              double sBondDist = 0.0;
              if (sB1 != singleBonds.end()) {
                sBondDist = singleBonds[bT1];
              }
              else if (sB2 != singleBonds.end()) {
                sBondDist = singleBonds[bT2];
              }

              double a1a2Dist = pAt1->getCoords()->dist(*pAt2->getCoords()) - sBondDist;
              if (a1a2Dist < mostLikelyAtomDist) {
                mostLikelyAtom = a2;
                mostLikelyAtomDist = a1a2Dist;
              }

              int bondIndex = indexAB(pAt1->getIndex(), pAt2->getIndex(), MAXATOMS);
              int bdPl = bondPlanar[bondIndex];
              if (hybridizations[a2] == 3 or bdPl) { // sp2
                bool gotDoubleBond = false;
                for (unsigned int k = 0; k < bdAtoms[a2].size(); k++) {
                  atom* pAt3 = atomList[bdAtoms[a2][k]];
                  pBond = pParent->getBond(pAt2, pAt3);
                  if (pBond->type == 2) {
                    gotDoubleBond = true;
                  }
                }
                if (!gotDoubleBond) {
                  if (Q[i] > 3 or Q[a2] > 3) continue;
                  if (atSymbols[i] == "C" and atSymbols[a2] == "C") {
                    pBond = pParent->getBond(pAt1, pAt2);
                    pBond->type = 2;
                    hybridizations[i] = 3;
                    hybridizations[a2] = 3;

                    errorMessage += " Setting C-C Bond " + i2s(pBond->atom1->getFileID())
                                    + "-" + i2s(pBond->atom2->getFileID()) +
                                    " to double \n";
                    bOK = true;
                    break;
                  }
                }
              }
            }
          }
          if (!bOK) {
            // Check for 1HFS & 1JH1
            if ((atSymbols[i] == "C") and (atSymbols[mostLikelyAtom] != "C") and (mostLikelyAtom > -1)) {
              atom* pAt2 = atomList[mostLikelyAtom];
              pBond = pParent->getBond(pAt1, pAt2);

              pBond->type = 2;
              hybridizations[i] = 3;
              hybridizations[mostLikelyAtom] = 3;

              errorMessage += " CHECK Setting C-" + pAt2->getElementSymbol() +
                              " Bond " + i2s(pAt1->getFileID()) +
                              "-" + i2s(pAt2->getFileID()) +
                              " to double \n";
              bOK = true;
              break;
            }
          }
        }
      }
    }

    if (errorMessage != "\n") {
      errorLogger.throwError("hybridize::runLabute Step 9.8 ", errorMessage, 4);
    }

    // Set hybridizations, formal charges for each atom
    this->setHybridizations();
    this->setFormalCharges();

    errorMessage = " FINAL :: FORMAL CHARGES | HYBRIDIZATIONS \n";
    for (unsigned int i = 0; i < nAtoms; i++) {
      errorMessage += "   " + i2s(i+1) + " " + i2s(formalCharges[i]) + " " + i2s(hybridizations[i]) + "\n";
    }

    errorMessage += " FINAL :: BOND ORDERS \n";
    errorMessage += " # I-J Type FileID_I FileID_J \n";
    if (!bonds.empty()) {
      int i = 1;
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        errorMessage +=   i2s(i) + " " + i2s(pBond->atom1->getIndex()) + "-" + i2s(pBond->atom2->getIndex())
                  + " " + i2s(pBond->type)
                  + "   " + i2s(pBond->atom1->getFileID()) + "-" + i2s(pBond->atom2->getFileID())
                  + " \n";
        i++;
      }
    }
    errorLogger.throwError("hybridize::runLabute Final Step", errorMessage, 4);
    errorLogger.flush();

    pParent->bAtomHybridizationsAssigned = true;

    return 0;
}

// ============================================================
// Function : COM()
// ------------------------------------------------------------
//
// ============================================================
void hybridize::COM(std::vector<atom*> atoms, vector3d* c)
{
    vector3d* Coord;
    vector3d  Center;

    for (unsigned int i = 0; i < atoms.size(); i++) {
      pAtom = atoms[i];
      Coord = pAtom->getCoords();
      Center = Center+(*Coord);
    }

    Center = Center / atoms.size();
    *c =  Center;
}

// ============================================================
// Function : setHybridizations()
// ------------------------------------------------------------
//
// ============================================================
void hybridize::setHybridizations()
{
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      pAtom->setHybridization(hybridizations[i]);
    }
}

// ============================================================
// Function : setFormalCharges()
// ------------------------------------------------------------
//
// ============================================================
void hybridize::setFormalCharges()
{
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      pAtom->setFormalCharge(formalCharges[i]);
    }
}

// ============================================================
// Function : runMeng()
// ------------------------------------------------------------
// Algorithm: J. Comp. Chem. 12, 891-898, 1991
// ============================================================
int hybridize::runMeng()
{
    std::vector<atom*> atomList = pParent->getAtomList();

    std::vector<int> reliableIndicator;
    for (unsigned int i = 0; i < atomList.size(); i++) {
      reliableIndicator.push_back(1);
    }

    for (unsigned int i = 0; i < atomList.size(); i++) {
      atom* pAtom = atomList[i];
      //std::string atSymbol = pAtom->getElementSymbol();
      int atNumber = pAtom->getAtomicNum();
      std::vector<atom*> bondedHeavyAtoms = pAtom->getBondedHeavyAtoms();

      switch (bondedHeavyAtoms.size()) {
        case 2 : {
          atom* bondedAtom1 = bondedHeavyAtoms[0];
          atom* bondedAtom2 = bondedHeavyAtoms[1];

          //double dist1 = pAtom->getCoords()->dist(*bondedAtom1->getCoords());
          //double dist2 = pAtom->getCoords()->dist(*bondedAtom2->getCoords());
          double ang = angle(*(pAtom->getCoords()), *(bondedAtom1->getCoords()), *(bondedAtom2->getCoords()));

          /*
             Carbon and Nitrogen may be sp3, sp2, or sp hybridized
             Oxygen and sulfur must be sp3 hybridized

               H      R
               |       \
             R-C-R ,    C=R , nitrile R-C%N
               |       /
               H     H

                  H
                /                                               + -             + -              +
             R-N     R=N-H , isocyanate R-N=C=O , isonitrile R-N%C , azide R-N=N=N , diazonium R-N%N
                \
                 R

             R-N=0 nitroso
          */
          switch (atNumber) {
            case 6 : { // sp3, sp2, or sp hybridized Carbon
              reliableIndicator[i] = 0;
              if (ang < AMBIGUOUSANGLE) {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("C3"); // sp3-hybridized Carbon
              }
              if (ang > SPANGLE) {
                pAtom->setHybridization(2); // sp
                pAtom->setMengType("C1"); // sp-hybridized Carbon
              }
              else if (ang > SP2ANGLE) {
                pAtom->setHybridization(3); // sp2
                pAtom->setMengType("C2"); // sp2-hybridized Carbon
              }
              break;
            }
            case 7 : { // sp3, sp2, or sp hybridized N
              reliableIndicator[i] = 0;
              if (ang < AMBIGUOUSANGLE) {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("N3"); // sp3-hybridized Nitrogen
              }
              if (ang > SPANGLE) {
                pAtom->setHybridization(2); // sp
                pAtom->setMengType("N1"); // sp-hybridized Nitrogen
              }
              else if (ang > SP2ANGLE) {
                pAtom->setHybridization(3); // sp2
                pAtom->setMengType("Npl"); // sp2-hybridized Nitrogen
              }
              break;
            }
            case 8 : { // ether, ester, R-O-N=O nitrite, R-O-C%N cyanate O
              pAtom->setHybridization(4);
              pAtom->setMengType("O3"); // sp3-hybridized Oxygen

              break;
            }
            case 16 : { // thioether, thioester, R-S-C%N thiocyanate
              pAtom->setMengType("S3"); // S3   sp3-hybridized Sulfur, neutral
              pAtom->setHybridization(4);
              break;
            }
            default : {
              std::cout << " case 2 unknown " << std::endl;
            }
          }
          break;
        }
        case 3 : {
          /*
            Carbon may be sp3, sp2 hybridized or be a part of a carboxylate group
            Nitrogen may be sp3, sp2 hybridized or part of a nitro group
            Sulfur may be positively charged and sp3 hybridized, or part of a sulfoxide group
            Boron atoms may be in a reduced or oxidized state
            Phosphorous .... phosphine [PR3]

            Atom after R is the current atom:
            carboxylate        amidino
                   O-               N
                 /                /
              O=C             +N=C
                 \                \
                  R                R

                    nitro          nitrate
                   O-        O         O-
                 /         /          /
              O=N+  or  O=N      R-O-N+
                \          \         \\
                 R          R         O

             sulfinic acid   sulfoxide
                   O              O
                 //             //
              R-S            R-S
                 \              \
                  OH             R

             boronic acid
                   OH
                 /
              R-B
                 \
                  OH

            aniline-like
                    Y
                  /
             N--X        X,Y,Z PLANAR.
                 \
                  Z

              P(OR)R2      P(OR)2R      P(OR)3
            Phosphinite  Phosphonite  Phosphite
                  OR         OR          OR
                  |          |           |
                  P          P           P
                 / \        / \         / \
                R   R      R   OR     RO   OR

          */
          atom* bondedAtom1 = bondedHeavyAtoms[0];
          atom* bondedAtom2 = bondedHeavyAtoms[1];
          atom* bondedAtom3 = bondedHeavyAtoms[2];

          double ang1 = angle(*(bondedAtom1->getCoords()), *(pAtom->getCoords()), *(bondedAtom2->getCoords()));
          double ang2 = angle(*(bondedAtom2->getCoords()), *(pAtom->getCoords()), *(bondedAtom3->getCoords()));
          double ang3 = angle(*(bondedAtom1->getCoords()), *(pAtom->getCoords()), *(bondedAtom3->getCoords()));
          double avgAng = (ang1 + ang2 + ang3)/3;
          int nOxygens = pAtom->numBondedOxygens();

          switch (atNumber) {
            case 6 : { // May be sp2, sp3 hybridized or be a part of a carboxylate group
              if (avgAng > SP2ANGLE) {
                pAtom->setHybridization(3); // sp2
                pAtom->setMengType("C2"); // sp2-hybridized Carbon
              }
              else {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("C3"); // sp3-hybridized Carbon
              }
              if (nOxygens == 2) {
                // carboxylate group
                pAtom->setMengType("Cac"); // carboxylate Carbon
              }
              break;
            }
            case 7 : { // May be sp2, sp3 hybridized or part of a nitro group
              if (avgAng > SP2ANGLE) {
                pAtom->setHybridization(3); // sp2
                pAtom->setMengType("Npl"); // sp2-hybridized N
              }
              else {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("N3"); // sp3-hybridized N (neutral)
              }
              if (nOxygens >= 2) {
                pAtom->setMengType("Ntr");  // nitro group
              }
              break;
            }
            case 16 : { // May be positively charged and sp3 hybridized,
                        // or part of a sulfoxide group
              pAtom->setMengType("S"); // default S

              if (avgAng > SP2ANGLE) {
                pAtom->setHybridization(3); // sp2
                pAtom->setMengType("S2"); // sp2-hybridized Sulfur
              }
              else {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("S3+"); // S3+  sp3-hybridized Sulfur, formal +ve charge
              }
              if (nOxygens == 1) {
                pAtom->setHybridization(3); // sulfoxide group
                pAtom->setMengType("Sox");
              }
              //if (nOxygens == 2) {
                // sulfinic acid
              //}
              if (nOxygens == 4) {
                 pAtom->setMengType("Sac"); // sulfate Sulfur
              }
              break;
            }
            case 5 : { // May be in a reduced or oxidized state
              /*
                - Bac  borate Boron
                - Box  other oxidized Boron
              */
              pAtom->setMengType("B"); // default Boron
              break;
            }
            case 15 : { // phosphine [PR3]
              /*
                - P3+  sp3-hybridized Phosphorous, formal +ve charge
              */
              pAtom->setHybridization(5); // sp3d
              pAtom->setMengType("P"); // default Phosphorous

              if (nOxygens > 1) {
                // Phosphinite
                pAtom->setMengType("Pox"); // P-oxide Phosphorous
              }
              break;
            }
            default : {
              std::cout << " case 3 unknown " << std::endl;
            }
          }
          break;
        }
        case 4 : {
          /*
            Carbon must be sp3 hybridized
            Nitrogen must be quaternary or oxidized
            Phosphorous must be part of a phosphate, P-oxide, or quaternary phosphine group
            Sulfur must be apart of a sulfate, sulfone, or sulfoxide group
            Boron may be part of a borate

              N-oxide
                 O-
                 |
               R-N+
                / \
               R   R

             Phosphine  Phosphinate  Phosphonate  Phosphate
               oxide
                  O           O            O            O
                 ||          ||           ||           ||
               R-P-R      RO-P-R       RO-P-OR      RO-P-OR
                 |           |            |            |
                 R           R            R            OR

             sulfonate  sulfone
                O-       O    O
               |         \\ //
             O=S=O        S
               |         / \
               R        R   R
          */
          // Get improper torsion, test if tetrahedral or planar
          //double dImproper = torsion(*(bondedAtoms[0]->getCoords()), *(bondedAtoms[1]->getCoords()),
          //                           *(bondedAtoms[2]->getCoords()), *(bondedAtoms[3]->getCoords()));

          int nOxygens = pAtom->numBondedOxygens();

          switch (atNumber) {
            case 5 : { // borate
              /*
                - Bac  borate Boron
                - Box  other oxidized Boron
              */
              pAtom->setHybridization(4);
              pAtom->setMengType("B"); // default Boron
              break;
            }
            case 6 : { // sp3 hybridized
              pAtom->setHybridization(4);
              pAtom->setMengType("C3"); // sp3-hybridized Carbon
              break;
            }
            case 7 : { // quaternary or oxidized
              pAtom->setHybridization(4); // sp3
              pAtom->setMengType("N3+");
              if (nOxygens == 1) {
                pAtom->setMengType("Nox"); // N-oxide
              }
              break;
            }
            case 15 : { // phosphate, P-oxide or quaternary phosphine group
              pAtom->setMengType("P"); // default Phosphorous
              pAtom->setHybridization(4);

              if (nOxygens == 1) {
                pAtom->setMengType("Pox"); // P-oxide Phosphorous
              }
              else if (nOxygens > 1) {
                pAtom->setMengType("Pac"); // phosphate Phosphorous
              }
              break;
            }
            case 16 : { // sulfate, sulfone, or sulfoxide
              pAtom->setHybridization(4);
              pAtom->setMengType("S"); // default S

              if (nOxygens > 1) {
                pAtom->setMengType("Sac"); // sulfate Sulfur
              }
              break;
            }
            default : {
              pAtom->setHybridization(4);
              std::cout << " case 4 unknown " << std::endl;
            }
          }
          break;
        }
        case 5 : {
          pAtom->setHybridization(1);
          pAtom->setMengType(pAtom->getElementSymbol());
          break;
        }
        case 6 : {
          pAtom->setHybridization(1);
          pAtom->setMengType(pAtom->getElementSymbol());
          break;
        }
        default : {
          std::cout << " > 6 bonded atoms ... unknown " << std::endl;
          pAtom->setHybridization(1);
          pAtom->setMengType(pAtom->getElementSymbol());
        }
      }
    }

    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      //std::string atSymbol = pAtom->getElementSymbol();
      int atNumber = pAtom->getAtomicNum();
      std::vector<atom*> bondedHeavyAtoms = pAtom->getBondedHeavyAtoms();

      atom* bondedAtom = 0;

      switch (bondedHeavyAtoms.size()) {
        case 1 : { // terminal atom, e.g. methyl group carbon; H; halogen.
          bondedAtom = bondedHeavyAtoms[0];
          std::string bdAtSymbol = bondedAtom->getElementSymbol();
          std::string bdAtMengType = bondedAtom->getMengType();
          double dist = pAtom->getCoords()->dist(*bondedAtom->getCoords());

          /*
            Atom after R is the current atom (Hydrogen maybe present):
             R-C-C, R-C=C, R-C%C, R-C-O, R-C=O, R-C-S, R-C=S, R-C-N, R-C=N, R-C%N
             R-O-C, R-OH, R-SH
             R-N-C, R-N=C, R-N-O, R-N=O, R-NH2, R-C=NH, R-N=N+, R-N=N=N
          */
          switch (atNumber) {
            case 1 : {
              if (bdAtSymbol == "C") {
                pAtom->setMengType("HC"); // H bonded to C
              }
              else {
                pAtom->setMengType("HC"); // default H
              }
            }
            case 6 : {
              if (bdAtSymbol == "C") { // C - C bond
                // C1 - C1 1.22 upper limit
                // C2 - C  1.41 upper limit
                // C3 - C  1.45 lower limit
                if (dist < 1.22 and bdAtMengType == "C1") { // C1-C1 bond
                  pAtom->setHybridization(2); // sp
                  pAtom->setMengType("C1");
                }
                else if (dist < 1.41) {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("C2");
                }
                else {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("C3");
                }
              }
              else if (bdAtSymbol == "N") { // C - N bond
                // C2 - N 1.37 upper limit
                if (dist < 1.37) {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("C2");
                }
                else {
                  pAtom->setHybridization(4); // sp2
                  pAtom->setMengType("C3");
                }
              }
              else if (bdAtSymbol == "O") { // C - O bond
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("C3");
              }
              else {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("C3");
              }
            }
            case 7 : {
              if (bdAtSymbol == "C") {
                // N1 - C1 1.20 upper limit
                // N3 - C 1.38 lower limit
                if (dist < 1.20 and bdAtMengType == "C1") {
                  pAtom->setHybridization(2); // sp
                  pAtom->setMengType("N1");
                }
                else if (dist > 1.38) {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("N3");
                }
              }
              if (bdAtSymbol == "N") {
                // N3 - N3 1.43 lower limit
                // N3 - Npl 1.41 lower limit
                if (dist > 1.42 and bdAtMengType == "N3") {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("N3");
                }
                if (dist > 1.41 and bdAtMengType == "Npl") {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("N3");
                }
              }
            }
            case 8 : {
              if (bdAtSymbol == "C") {
                // O2 - C2 1.30 upper limit
                if (dist < 1.30 and bdAtMengType == "C2") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("O2");
                }
              }
              else if (bdAtSymbol == "As") {
                // O2 - As 1.685 upper limit
                if (dist < 1.685) {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("O2");
                }
              }
              else if (bdAtMengType == "Cac" or bdAtMengType == "Ntr" or
                       bdAtMengType == "Son" or bdAtMengType == "Sxd") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("O2-"); // carboxylate, nitro or nitrate O, partial negative charge
              }
              else if (bdAtMengType == "Sac" or bdAtMengType == "Pac") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("O3-"); // resonance-equivalent terminal O on tetrahedral atom: phosphate, sulfate
              }
              else {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("O3"); // sp3 hybridized O
              }
            }
            case 16 : {
              if (bdAtSymbol == "C") {
                // S2 - C2 1.76 upper limit
                if (dist < 1.76 and bdAtMengType == "C2") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("S2");
                }
              }
              else if (bdAtSymbol == "As") {
                // S2 - As 2.11 upper limit
                if (dist < 2.11) {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("O2");
                }
              }
              else {
                pAtom->setHybridization(4); // sp3
                pAtom->setMengType("S3"); // sp3 hybridized S
              }
            }

            default : {
              pAtom->setHybridization(1);
              pAtom->setMengType(pAtom->getElementSymbol());
            }
          }
          break;
        }
      }
    }

    for (unsigned int i = 0; i < atomList.size(); i++) {
      pAtom = atomList[i];
      std::string atSymbol = pAtom->getElementSymbol();
      std::string atMengType = pAtom->getMengType();
      std::vector<atom*> bondedHeavyAtoms = pAtom->getBondedHeavyAtoms();

      switch (bondedHeavyAtoms.size()) {
        case 2 : {
          if (!reliableIndicator[i]) {
            if (atSymbol == "C") {
              for (unsigned int j = 0; j < 2; j++) {
                double cDist = pAtom->getCoords()->dist(*bondedHeavyAtoms[j]->getCoords());
                std::string sy = bondedHeavyAtoms[j]->getElementSymbol();
                if (cDist > 1.53 and sy == "C") {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("C3"); // sp3 hybridized C
                }
                else if (cDist < 1.42 and sy == "C" and atMengType != "C1") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("C2");
                }

                if (cDist > 1.46 and sy == "N") {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("C3"); // sp3 hybridized C
                }
                else if (cDist < 1.41 and sy == "N") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("C2");
                }

                if (cDist > 1.44 and sy == "O") {
                  pAtom->setHybridization(4); // sp3
                  pAtom->setMengType("C3"); // sp3 hybridized C
                }
              }
            }
            else if (atSymbol == "N") {
              for (unsigned int j = 0; j < 2; j++) {
                double cDist = pAtom->getCoords()->dist(*bondedHeavyAtoms[j]->getCoords());
                std::string sy = bondedHeavyAtoms[j]->getElementSymbol();
                if (cDist < 1.38 and sy == "C") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("Npl"); // sp2 hybridized N
                }
                if (cDist < 1.32 and sy == "N") {
                  pAtom->setHybridization(3); // sp2
                  pAtom->setMengType("Npl"); // sp2 hybridized N
                }
              }
            }
          }
        }
      }
    }
    pParent->bAtomHybridizationsAssigned = true;
    return 0;
}

} // MTKpp namespace
