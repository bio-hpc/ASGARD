/*!
   \file ring.cpp
   \brief Determines rings in a molecule
   \author Martin Peters

   $Date: 2010/04/29 18:59:18 $
   $Revision: 1.22 $

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

#include "ring.h"
#include "molecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"
#include "utility.h"

#include "Utils/constants.h"

#include "Utils/vector3d.h"
#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

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
// Function : rings(molecule*)
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
rings::rings(molecule *parent):pParent(parent) {}

// ============================================================
// Function : ~rings()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
rings::~rings() {}

// ============================================================
// Function : determine()
// ------------------------------------------------------------
// Determine all rings in the molecule
// ============================================================
void rings::determine()
{
    molAtomList = pParent->getAtomList();
    if (molAtomList.size() < 3) return;

    molBondMap  = pParent->getBondMap();
    int frerejaque = molBondMap.size() - molAtomList.size() + 1;

#ifdef DEBUG
    std::string errMessage =   " Number of Atoms (a) = " + i2s(molAtomList.size())
                           + "\n Number of Bonds (b) = " + i2s(molBondMap.size())
                           + "\n Frerejaque Number (b - a + 1) = " + i2s(frerejaque) + "\n";
    errorLogger.throwError("rings::determine", errMessage, 4);
#endif

    // if rings are present
    if (frerejaque > 0) {
      this->removeHydrogenAtoms();

      this->removeHydrogenBonds();

      this->removeOpenAcyclic();

      this->removeOpenAcyclicBonds();

      this->findNeighbors();

      this->removeClosedAcyclic();

      this->findNeighbors();

      this->separateBlocks();

      this->findSSSR();
    }
    else {
      this->removeHydrogenAtoms();

      this->removeHydrogenBonds();

      this->removeOpenAcyclic();

      this->removeOpenAcyclicBonds();
    }

    // Set the bond topology
    std::vector<ring*> molRings = pParent->getRings();
    ring* r;
    Bond* pRingBond;
    unsigned int s;
    for (unsigned int i = 0; i < molRings.size(); i++) {
      r = molRings[i];
      s = r->atoms.size();
      for (unsigned int j = 0; j < s; j++) {
        if (j != s-1) {
          pRingBond = pParent->getBond(r->atoms[j], r->atoms[j+1]);
          pRingBond->topology = 1;
        }
        else {
          pRingBond = pParent->getBond(r->atoms[j], r->atoms[0]);
          pRingBond->topology = 1;
        }
      }
    }

    // control
    // The number of the smallest rings in a block can be
    // calculated using the following:
    //      Nr = 1 + n(N3)/2 + n(N4)
    // J. Chem. Inf. Comput. Sci. 2000, 40, 1015-1017
    std::vector<Bond*> blockBonds;
    std::vector<int> nConnections;

    int numberRings = 0;
    int N3 = 0;
    int N4 = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      this->getBonds(blocks[i],cyclicBonds,blockBonds,nConnections);
      for (unsigned int j = 0; j < nConnections.size(); j++) {
        if (nConnections[j] == 3) N3++;
        if (nConnections[j] == 4) N4++;
      }
      numberRings += 1 + N3/2 + N4;
      N3 = 0;
      N4 = 0;
    }

    std::string errMessage2 =   " Number of Rings = " + i2s(numberRings) + "\n";
    errorLogger.throwError("rings::determine", errMessage2, 4);
}

// ============================================================
// Function : removeHydrogenAtoms()
// ------------------------------------------------------------
// Removes all hydrogens
// ============================================================
void rings::removeHydrogenAtoms()
{
    for (unsigned int i = 0; i < molAtomList.size(); i++) {
      pAtom = molAtomList[i];
      if (pAtom->getAtomicNum() != 1) {
        cyclicAtoms.push_back(pAtom);
      }
      else {
        pAtom->setType(1);
        hydrogens.push_back(pAtom);
      }
    }
}

// ============================================================
// Function : removeHydrogenBonds()
// ------------------------------------------------------------
// Removes all hydrogens
// ============================================================
void rings::removeHydrogenBonds()
{
    for (BondMapIterator b = molBondMap.begin(); b != molBondMap.end(); b++) {
      pBond = b->second;
      if ((pBond->atom1->getAtomicNum() == 1) or (pBond->atom2->getAtomicNum() == 1)) {
        continue;
      }
      else {
        pBond->topology = 2;
        cyclicBonds.push_back(pBond);
      }
    }
}

// ============================================================
// Function : removeOpenAcyclic()
// ------------------------------------------------------------
// Removes all OPEN ACYCLIC NODES (terminal atoms)
// ============================================================
void rings::removeOpenAcyclic()
{
    int numBonds = 0;
    bool done = false;
    bool terminal = false;
    std::vector<atom*> bondedAtoms;
    atom* pAtom2;

    while (!done) {
      done = true;
      for (unsigned int i = 0; i < cyclicAtoms.size(); i++ ) {
        terminal = false;
        pAtom = cyclicAtoms[i];
        bondedAtoms = pAtom->getBondedAtoms();
        numBonds = bondedAtoms.size();
        //std::cout << pAtom->getFileID() << " has " << numBonds << " bond/s" << std::endl;
        // Determine if atom is a true terminal atom or a open chain atom
        int testTerminal = 0;
        for (unsigned int j = 0; j < bondedAtoms.size();j++) {
          pAtom2 = bondedAtoms[j];
          if (pAtom2->getAtomicNum() != 1) {
            testTerminal++;
          }
        }
        if (testTerminal == 1 ) {
          terminal = true;
        }

        for (unsigned int j = 0; j < bondedAtoms.size();j++) {
          pAtom2 = bondedAtoms[j];
          atomIterator = std::find(hydrogens.begin(), hydrogens.end(), pAtom2);
          if (atomIterator != hydrogens.end()) {
            numBonds--;
          }
          atomIterator = std::find(acyclicAtoms.begin(), acyclicAtoms.end(), pAtom2);
          if (atomIterator != acyclicAtoms.end()) {
            numBonds--;
          }
        }
        if (numBonds <= 1) {
          atomIterator = std::find(acyclicAtoms.begin(), acyclicAtoms.end(), pAtom);
          if (atomIterator == acyclicAtoms.end()) {

            if (terminal) {
              pAtom->setType(2);
            }
            else {
              pAtom->setType(3);
            }
            acyclicAtoms.push_back(pAtom);
            done = false;
          }
        }
        numBonds = 0;
      }
    }

    for (unsigned int i = 0; i < acyclicAtoms.size(); i++) {
      pAtom = acyclicAtoms[i];
      atomIterator = std::find(cyclicAtoms.begin(), cyclicAtoms.end(), pAtom);
      if (atomIterator != cyclicAtoms.end()) {
        cyclicAtoms.erase(atomIterator);
      }
    }
}

// ============================================================
// Function : removeOpenAcyclicBonds()
// ------------------------------------------------------------
// Removes all OPEN ACYCLIC bonds
// ============================================================
void rings::removeOpenAcyclicBonds()
{
    //std::vector<Bond*>::iterator bondResult;

    for (unsigned int i = 0; i < acyclicAtoms.size(); i++) {
      pAtom = acyclicAtoms[i];
      for (unsigned int j = 0; j < acyclicAtoms.size(); j++) {
        pAtom2 = acyclicAtoms[j];
        pBond = pParent->getBond(pAtom,pAtom2);
        if (pBond) {
          bondIterator = std::find(cyclicBonds.begin(), cyclicBonds.end(), pBond);
          if (bondIterator != cyclicBonds.end()) {
            pBond->topology = 2;
            cyclicBonds.erase(bondIterator);
          }
        }
      }
      for (unsigned int j = 0; j < cyclicAtoms.size(); j++) {
        pAtom2 = cyclicAtoms[j];
        pBond = pParent->getBond(pAtom,pAtom2);
        if (pBond) {
          bondIterator = std::find(cyclicBonds.begin(), cyclicBonds.end(), pBond);
          if (bondIterator != cyclicBonds.end()) {
            pBond->topology = 2;
//#ifdef DEBUG
//        std::cout << " rings:removeOpenAcyclicBonds -> Erasing Terminal Bond = " 
//                  << pBond->atom1->getFileID() << "-" << pBond->atom2->getFileID() << std::endl;
//#endif
            cyclicBonds.erase(bondIterator);
          }
        }
      }
    }
}

// ============================================================
// Function : findNeighbors()
// ------------------------------------------------------------
// Find all neighboring atoms
// ============================================================
void rings::findNeighbors()
{
    neighbors.clear();
    for (unsigned int i = 0; i < cyclicAtoms.size(); i++) {
      pAtom = cyclicAtoms[i];
      std::vector<int> n;
      for (unsigned int j = 0; j < cyclicAtoms.size(); j++) {
        pAtom2 = cyclicAtoms[j];
        pBond = pParent->getBond(pAtom, pAtom2);
        if (pBond) {
          n.push_back(j);
        }
      }
      neighbors.push_back(n);
    }
}

// ============================================================
// Function : findNeighbors(vector<atom*> a, vector<int> b)
// ------------------------------------------------------------
// Find all neighboring atoms
// ============================================================
void rings::findNeighbors(std::vector<atom*> a, std::vector< std::vector<int> > &b)
{
    for (unsigned int i = 0; i < a.size(); i++) {
      pAtom = a[i];
      std::vector<int> n;
      for (unsigned int j = 0; j < a.size(); j++) {
        pAtom2 = a[j];
        pBond = pParent->getBond(pAtom, pAtom2);
        if (pBond) {
          n.push_back(j);
        }
      }
      b.push_back(n);
    }
}

// ============================================================
// Function : getBonds()
// ------------------------------------------------------------
// Find all bonds between atoms in a from the list b
// ============================================================
void rings::getBonds(std::vector<atom*> a, std::vector<Bond*> b,
                     std::vector<Bond*> &c, std::vector<int> &d)
{
    int i = 0;
    c.clear();
    d.clear();

    for (unsigned int q = 0; q < a.size(); q++) {
      atom* pAtomQ = a[q];
      for (unsigned int p = 0; p < a.size(); p++) {
        atom* pAtomP = a[p];
        for (unsigned int f = 0; f < b.size(); f++) {
          Bond* pBondF = b[f];
          if ((pBondF->atom1 == pAtomQ and pBondF->atom2 == pAtomP) or
              (pBondF->atom1 == pAtomP and pBondF->atom2 == pAtomQ) ) {
            bondIterator = std::find(c.begin(), c.end(), pBondF);
            if (bondIterator == c.end()) {
              c.push_back(pBondF);
            }
            i++;
          }
        }
      }
      d.push_back(i);
      i = 0;
    }
}

// ============================================================
// Function : removeClosedAcyclic()
// ------------------------------------------------------------
// Removes all closed acyclic atoms (atoms joining ring systems)
// ============================================================
void rings::removeClosedAcyclic()
{
    // color: 0 == white
    // color: 1 == grey
    // color: 2 == black
    std::vector<int> vertexesColor;
    for (unsigned int i = 0; i < cyclicAtoms.size(); i++) {
      vertexesColor.push_back(0);
    }

    std::vector<int> edgesColor;
    for (unsigned int i = 0; i < cyclicBonds.size(); i++) {
      edgesColor.push_back(0);
    }

    std::vector< std::vector<atom*> > paths;
    std::vector<atom*> path;
    bool loop = false;

    for (unsigned int i = 0; i < cyclicAtoms.size(); i++) {
      std::vector<int> curNeighbors = neighbors[i];
      for (unsigned int j = 0; j < curNeighbors.size(); j++) {
        vertexesColor[i] = 1;
        int curNeighbor = curNeighbors[j];
        int t = 0;
        dfs_visit(curNeighbor,i,cyclicAtoms,cyclicBonds,neighbors,vertexesColor,edgesColor,loop,t,path,paths);
        t = 0;

        for (unsigned int k = 0; k < cyclicAtoms.size(); k++) {
          vertexesColor[k] = 0;
        }
        for (unsigned int l = 0; l < cyclicBonds.size(); l++) {
          edgesColor[l] = 0;
        }
        if (loop) {
          path.clear();
          break;
        }
      }
      if (!loop) {
        cyclicAtoms[i]->setType(4);
        acyclicAtoms.push_back(cyclicAtoms[i]);
        path.clear();
        continue;
      }
      loop = false;
    }

    std::vector<atom*>::iterator result;
    for (unsigned int i = 0; i < acyclicAtoms.size(); i++) {
      pAtom = acyclicAtoms[i];
      result = std::find(cyclicAtoms.begin(), cyclicAtoms.end(), pAtom);
      if (result != cyclicAtoms.end()) {
        cyclicAtoms.erase(result);
      }
    }
    this->removeOpenAcyclicBonds();
}

// ============================================================
// Function : separateBlocks()
// ------------------------------------------------------------
// Separate blocks
// ============================================================
void rings::separateBlocks()
{
    std::vector<int> vertexesColor;
    for (unsigned int i = 0; i < cyclicAtoms.size(); i++) {
      vertexesColor.push_back(0);
    }

    std::vector<int> edgesColor;
    for (unsigned int i = 0; i < cyclicBonds.size(); i++) {
      edgesColor.push_back(0);
    }

    std::vector< std::vector<atom*> > paths;
    std::vector<atom*> path;
    std::vector<atom*> block;
    std::vector<atom*> b;

    int curNeighbor;

    bool loopMatrix[cyclicAtoms.size()][cyclicAtoms.size()];
    for (unsigned int k = 0; k < cyclicAtoms.size(); k++) {
      for (unsigned int h = 0; h < cyclicAtoms.size(); h++) {
        loopMatrix[k][h] = false;
      }
    }

    bool loop = false;
    for (unsigned int i = 0; i < cyclicAtoms.size(); i++) {
      std::vector<int> curNeighbors = neighbors[i];
      for (unsigned int j = 0; j < curNeighbors.size(); j++) {
        vertexesColor[i] = 1;
        curNeighbor = curNeighbors[j];
        int t = 0;
        dfs_visit(curNeighbor,i,cyclicAtoms,cyclicBonds,neighbors,vertexesColor,edgesColor,loop,t,path,paths);
        loopMatrix[i][curNeighbor] = loop;
        t = 0;

        for (unsigned int k = 0; k < cyclicAtoms.size(); k++) {
          vertexesColor[k] = 0;
        }
        for (unsigned int l = 0; l < cyclicBonds.size(); l++) {
          edgesColor[l] = 0;
        }
        path.clear();
        loop = false;
      }
    }

    block.push_back(cyclicAtoms[0]);
    blocks.push_back(block);
    int curBlock = 0;

    unsigned int sumInBlocks = 1;
    bool someAdded = false;
    int k = 0;

    while (sumInBlocks < cyclicAtoms.size()) {
      for (unsigned int i = 0; i < blocks[curBlock].size(); i++) {
        for (unsigned int j = 0; j < cyclicAtoms.size(); j++) {
          if (blocks[curBlock][i] == cyclicAtoms[j]) k = j;
        }
        for (unsigned int j = 0; j < cyclicAtoms.size(); j++) {
          if (loopMatrix[k][j] == true) {
            atomIterator = std::find(blocks[curBlock].begin(), blocks[curBlock].end(), cyclicAtoms[j]);
            if (atomIterator == blocks[curBlock].end()) {
              blocks[curBlock].push_back(cyclicAtoms[j]);
              someAdded = true;
            }
          }
        }
        k = 0;
      }

      sumInBlocks = 0;
      for (unsigned int x = 0; x < blocks.size(); x++) {
        sumInBlocks += blocks[x].size();
      }

      if (sumInBlocks == cyclicAtoms.size()) break;

      if (someAdded) {
        someAdded = false;
      }
      else {
        atom* a = 0;
        bool gotIt = false;
        for (unsigned int i = 1; i < cyclicAtoms.size(); i++) {
          //std::cout << cyclicAtoms[i]->getFileID() << std::endl;
          for (unsigned int j = 0; j < blocks.size(); j++) {
            for (unsigned int j2 = 0; j2 < blocks[j].size(); j2++) {
              if (cyclicAtoms[i] == blocks[j][j2]) {
                gotIt = true;
                //std::cout << " gotIT " << std::endl;
              }
            }
          }
          if (!gotIt) {
            //std::cout << " !gotIt " << cyclicAtoms[i]->getFileID() << std::endl;
            a = cyclicAtoms[i];
            break;
          }
          gotIt = false;
        }
        curBlock++;
        std::vector<atom*> c;
        c.push_back(a);
        blocks.push_back(c);
      }
    }

    // REMOVE INTER BLOCK BONDS
    std::vector<atom*> block1;
    std::vector<atom*> block2;
    atom* block1Atom;
    atom* block2Atom;
/*
#ifdef DEBUG
   std::cout << blocks.size() << std::endl;
   for (unsigned int i = 0; i < blocks.size(); i++) {
     std::cout << " ring::block: " << std::endl;
     for (unsigned int j = 0; j < blocks[i].size(); j++) {
       std::cout << blocks[i][j]->getFileID() << std::endl;
     }
   }
#endif
*/
    if (blocks.size() > 1) {
      for (unsigned int i = 0; i < blocks.size(); i++) {
        block1 = blocks[i];
        for (unsigned int j = i+1; j < blocks.size(); j++) {
          block2 = blocks[j];
          for (unsigned int d = 0; d < block1.size(); d++) { // k to d
            block1Atom = block1[d];
            for (unsigned int l = 0; l < block2.size(); l++) {
              block2Atom = block2[l];
              pBond = pParent->getBond(block1Atom,block2Atom);
              if (pBond) {
                bondIterator = std::find(cyclicBonds.begin(), cyclicBonds.end(), pBond);
                if (bondIterator != cyclicBonds.end()) {
                  cyclicBonds.erase(bondIterator);
                  //std::cout << block1Atom->getName() << " " << block2Atom->getName() << std::endl;
                }
              }
            }
          }
        }
      }
    }
}

// ============================================================
// Function : pickRootAtom()
// ------------------------------------------------------------
// Pick Root Atom
// ============================================================
int rings::pickRootAtom(std::vector<atom*> a)
{
    int root = 0;
    int nBonds = 0;

    typedef std::vector<int>::iterator iVectorIt ;
    iVectorIt start, end;

    std::vector<int> possibleRoots;

    for (unsigned int i = 0; i < a.size(); i++) {
      pAtom = a[i];
      for (unsigned int j = 0; j < a.size(); j++) {
        pAtom2 = a[j];
        pBond = pParent->getBond(pAtom, pAtom2);
        if (pBond) {
          nBonds++;
        }
      }
      if (nBonds == 2) {
        possibleRoots.push_back(i);
      }
      nBonds = 0;
    }
    start = possibleRoots.begin();
    end = possibleRoots.end();
    std::random_shuffle(start, end);
    root = possibleRoots[0];
    return root;
}
// ============================================================
// Function : findSSSR()
// ------------------------------------------------------------
// Find the Smallest Set of Smallest Rings (SSSR)
// ============================================================
void rings::findSSSR()
{
    for (unsigned int i = 0; i < blocks.size(); i++) {
      this->decomposeBlock(blocks[i]);
    }
}

// ============================================================
// Function : decomposeBlock()
// ------------------------------------------------------------
//
// ============================================================
void rings::decomposeBlock(std::vector<atom*> block)
{
    if (this->numberRingsInBlock(block) > 1) {
      for (unsigned int i = 0; i < block.size(); i++) {
        //int root = this->pickRootAtom(block);
        std::vector<atom*> path;
        this->getIrreducibleClosedPath(block, i, path);

        if (path.size() > 2) {
/*
#ifdef DEBUG
          std::cout << " Try To Add Ring:";
          for (unsigned int cv = 0; cv < path.size(); cv++) {
            std::cout << " " << path[cv]->getFileID();
          }
          std::cout << " \n" << std::endl;
#endif
*/
          pParent->addRing(path);
        }
      }
    }
    else {
      std::vector<atom*> path;
      this->getIrreducibleClosedPath(block, 0, path);
/*
#ifdef DEBUG
      std::cout << " Try To Add Ring:";
      for (unsigned int cv = 0; cv < path.size(); cv++) {
        std::cout << " " << path[cv]->getFileID();
      }
      std::cout << " \n" << std::endl;
#endif
*/
      pParent->addRing(path);
    }
}

// ============================================================
// Function : getIrreducibleClosedPath()
// ------------------------------------------------------------
//
// ============================================================
void rings::getIrreducibleClosedPath(std::vector<atom*> block, int root, std::vector<atom*>& path)
{
    bool loop = false;
    std::vector<Bond*> blockBonds;
    std::vector< std::vector<int> > blockNeighbors;
    std::vector<int> vertexesColor;
    std::vector<int> edgesColor;
    std::vector<atom*> curPath;
    std::vector<int> nConnections;

    blockBonds.clear();
    nConnections.clear();
    blockNeighbors.clear();
    vertexesColor.clear();
    edgesColor.clear();
    this->getBonds(block,cyclicBonds,blockBonds,nConnections);
    this->findNeighbors(block, blockNeighbors);

    for (unsigned int tI = 0; tI < block.size(); tI++) {
      vertexesColor.push_back(0);
    }
    for (unsigned int tI = 0; tI < blockBonds.size(); tI++) {
      edgesColor.push_back(0);
    }
    int curAtom = blockNeighbors[root][0];

    bool done = false;
    bool firstTime = true;
    int nN3 = 0;
    std::vector<int> nN3s;
    while (!done) {

      if (firstTime) {
        bool ok = false;
        while (!ok) {
          curPath.push_back(block[root]);
          curPath.push_back(block[curAtom]);
          vertexesColor[root] = 1; // black
          Bond* pBond;
          for (unsigned int j = 0; j < blockBonds.size(); j++) {
            pBond = blockBonds[j];
            if ((pBond->atom1 == block[root] and pBond->atom2 == block[curAtom]) or
                (pBond->atom2 == block[root] and pBond->atom1 == block[curAtom]) ) {
              edgesColor[j] = 1; // black
            }
          }
          this->dfs_visitNEW(curAtom,root,1,loop,block,blockBonds,
                           blockNeighbors,vertexesColor,edgesColor,curPath);

          this->checkPath(curPath, ok);

          if (!ok) {
            //for (unsigned int cv = 0; cv < curPath.size(); cv++) {
            //  std::cout << " '" << curPath[cv]->getFileID() << "'";
            //}
            //std::cout << " " << std::endl;

            loop = false;
            curPath.clear();
            for (unsigned int tI = 0; tI < vertexesColor.size(); tI++) {
              vertexesColor[tI] = 0;
            }
            for (unsigned int tI = 0; tI < edgesColor.size(); tI++) {
              edgesColor[tI] = 0;
            }
          }
        }

        path = curPath;

        //std::cout << " First curPath:";
        //for (unsigned int cv = 0; cv < curPath.size(); cv++) {
        //  std::cout << " '" << curPath[cv]->getFileID() << "'";
        //}
        //std::cout << " " << std::endl;

        firstTime = false;
        int index;
        std::vector<atom*> tempPath = curPath;
        curPath.clear();
        for (unsigned int q = 0; q < tempPath.size(); q++) {
          for (unsigned int e = 0; e < block.size(); e++) {
            if (tempPath[q] == block[e]) {
              index = e;
            }
          }
          if (nConnections[index] > 2) {
            if (q != 0) {
              curAtom = index;
              nN3 = index;
              nN3s.push_back(index);
              curPath.push_back(block[index]);
              break;
            }
            else {
              curPath.push_back(block[index]);
            }
          }
          else {
            curPath.push_back(block[index]);
          }
        }
        if (curPath.size() == path.size()) done = true;
      }
      else {
        loop = false;

        if (block[curAtom] != path[path.size()-1]) {

          //std::cout << " curPath beforeX:";
          //for (unsigned int cv = 0; cv < curPath.size(); cv++) {
          //  std::cout << " " << curPath[cv]->getFileID();
          //}
          //std::cout << " " << std::endl;

          this->dfs_visitNEW(curAtom,root,0,loop,block,blockBonds,
                             blockNeighbors,vertexesColor,edgesColor,curPath);

          if (curPath.size() == 0) {
            curPath = path;
          }
          //std::cout << " curPath afterX:";
          //for (unsigned int cv = 0; cv < curPath.size(); cv++) {
          //  std::cout << " " << curPath[cv]->getFileID();
          //}
          //std::cout << " " << std::endl;

          if (curPath.size() < path.size()) {
            this->getNewPath(path,curPath,done);
            curPath.clear();
            curPath = path;
          }
          else {
            curPath.clear();
            curPath = path;
            //std::cout << " curPath reset:";
            //for (unsigned int cv = 0; cv < curPath.size(); cv++) {
            //  std::cout << " " << curPath[cv]->getFileID();
            //}
            //std::cout << " " << std::endl;
          }
        }
        else {
          done = true;
        }

        if (!done) {
          std::vector<atom*> tempPath = curPath;
          curPath.clear();
          int index;
          //int pastPreviousN3 = 0;
          bool newN3 = false;
          for (unsigned int q = 0; q < tempPath.size(); q++) {
            for (unsigned int e = 0; e < block.size(); e++) {
              if (tempPath[q] == block[e]) {
                index = e;
              }
            }

            if (nConnections[index] > 2) {
              if (q != 0) {
                std::vector<int>::iterator iVectorIt;
                iVectorIt = std::find(nN3s.begin(), nN3s.end(), index);
                if (iVectorIt == nN3s.end()) {
                  curAtom = index;
                  nN3 = index;
                  newN3 = true;
                  nN3s.push_back(index);
//std::cout << " adding to curPath: "<< block[index]->getFileID() << std::endl;
                  curPath.push_back(block[index]);
                  break;
                }
//std::cout << " adding to curPath: "<< block[index]->getFileID() << std::endl;
                curPath.push_back(block[index]);
              }
              else {
//std::cout << " adding to curPath: "<< block[index]->getFileID() << std::endl;
                curPath.push_back(block[index]);
              }
            }
            else {
//std::cout << " adding to curPath: "<< block[index]->getFileID() << std::endl;
              curPath.push_back(block[index]);
            }
          }
          if (!newN3) done = true;
        }
      }
    }
}

// ============================================================
// Function : numberRingsInBlock()
// ------------------------------------------------------------
//
// ============================================================
int rings::numberRingsInBlock(std::vector<atom*> block)
{
    std::vector<Bond*> blockBonds;
    std::vector<int> nConnections;

    int numberRings = 0;
    int N3 = 0;
    int N4 = 0;

    this->getBonds(block,cyclicBonds,blockBonds,nConnections); // dangerous use of cyclicBonds
    for (unsigned int j = 0; j < nConnections.size(); j++) {
       if (nConnections[j] == 3) N3++;
       if (nConnections[j] == 4) N4++;
    }
    numberRings = 1 + N3/2 + N4;
    return numberRings;
}

// ============================================================
// Function : getNewPath()
// ------------------------------------------------------------
//
// ============================================================
void rings::getNewPath(std::vector<atom*> &path, std::vector<atom*> curPath, bool &done)
{
    std::vector<atom*> newPath;
    int lastPoint = 0;
    bool error = false;
    for (unsigned int i = 0; i < curPath.size(); i++) {
      newPath.push_back(curPath[i]);
      lastPoint = i;
    }
    bool record = false;
    for (unsigned int i = 0; i < path.size(); i++) {
      if (record) {
        if (pParent->hasBond(path[i], newPath[newPath.size()-1])) {
          newPath.push_back(path[i]);
        }
        else {
          // not in sequence
          error = true;
          break;
        }
      }
      if (path[i] == newPath[lastPoint]) {
        record = true;
      }
    }
    if (error) return;

    if (newPath.size() < path.size()) {
      path.clear();
      path = newPath;
    }
    //else {
    //  done = true;
    //}

    //std::cout << " getNewPath: path:";
    //for (unsigned int cv = 0; cv < path.size(); cv++) {
    //  std::cout << " " << path[cv]->getFileID();
    //}
    //std::cout << " ....done = " << done << std::endl;
}

// ============================================================
// Function : checkPath()
// ------------------------------------------------------------
//
// ============================================================
void rings::checkPath(std::vector<atom*> &path, bool &ok)
{
    ok = true;
    Bond* pRingBond;
    for (unsigned int i = 0; i < path.size(); i++) {
//#ifdef DEBUG
//    std::cout << " rings:checkPath ->" << path[i] << " "<< path[i]->getFileID() << std::endl;
//#endif
      if (i != path.size()-1) {
        pRingBond = pParent->getBond(path[i], path[i+1]);
        if (!pRingBond) {
          ok = false;
        }
      }
      else {
        pRingBond = pParent->getBond(path[i], path[0]);
        if (!pRingBond) {
          ok = false;
        }
      }
    }
}

// ============================================================
// Function : eliminateReducibleAtoms()
// ------------------------------------------------------------
//
// ============================================================
void rings::eliminateReducibleAtoms(std::vector<atom*> block, std::vector<atom*> path)
{

}

// ============================================================
// Function : dfs_visit()
// ------------------------------------------------------------
//
// ============================================================
void rings::dfs_visit(int curAtom, int rootAtom,
//                      std::vector<atom*> vertexes, std::vector<Bond*> edges,
//                      std::vector< std::vector<int> > neighbors,
                      std::vector<atom*> &vertexes, std::vector<Bond*> edges,
                      std::vector< std::vector<int> > &neighbors,
                      std::vector<int>& vertexesColor, std::vector<int>& edgesColor,
                      bool& loop, int& t,
                      std::vector<atom*>& curPath,
                      std::vector< std::vector<atom*> >& paths)
{
    vertexesColor[curAtom] = 1; // grey
    std::vector<int> curNeighbors = neighbors[curAtom];
    int curNeighbor;
    int ec = 0;
    int curEdge;

    t++;
    if (!loop) {

      typedef std::vector<int>::iterator iVectorIt ;
      iVectorIt start, end;
      start = curNeighbors.begin();
      end = curNeighbors.end();
      std::random_shuffle(start, end);

      for (unsigned int i = 0; i < curNeighbors.size(); i++) {
        curNeighbor = curNeighbors[i];

        for (unsigned int j = 0; j < edges.size(); j++) {
          pBond = edges[j];
          if ((pBond->atom1 == vertexes[curAtom] and pBond->atom2 == vertexes[curNeighbor]) or
              (pBond->atom2 == vertexes[curAtom] and pBond->atom1 == vertexes[curNeighbor]) ) {
            curEdge = j;
            ec = edgesColor[j];
            edgesColor[j] = 2; // black
          }
        }
        if (vertexesColor[curNeighbor] == 0) {
          curPath.push_back(vertexes[curNeighbor]);
          this->dfs_visit(curNeighbor, rootAtom, vertexes, edges, neighbors,
                          vertexesColor,edgesColor,loop,t,curPath,paths);
        }

        if (t > 1) {
          if (ec == 0 and vertexesColor[curNeighbor] == 2) {
            if (pParent->hasBond(vertexes[rootAtom],vertexes[curNeighbor])) {
              loop = true;
              paths.push_back(curPath);
            }
          }
        }
      }
      vertexesColor[curAtom] = 2;
    }
}

// ============================================================
// Function : dfs_visitNEW()
// ------------------------------------------------------------
//
// ============================================================
void rings::dfs_visitNEW(int curAtom, int rootAtom, bool first, bool &loop,
//                         std::vector<atom*> vertexes, std::vector<Bond*> edges,
//                         std::vector< std::vector<int> > neighbors,
                         std::vector<atom*> &vertexes, std::vector<Bond*> edges,
                         std::vector< std::vector<int> > &neighbors,
                         std::vector<int>& vertexesColor, std::vector<int>& edgesColor,
                         std::vector<atom*>& path)
{
    vertexesColor[curAtom] = 1; // black
    std::vector<int> curNeighbors = neighbors[curAtom];
    int curNeighbor;
    int curEdge = 0;
    int bondIndex = 0;

    if (loop) return;
    //std::cout << " \ncurAtom = " << vertexes[curAtom]->getFileID() << std::endl;

    typedef std::vector<int>::iterator iVectorIt ;
    iVectorIt start, end;
    start = curNeighbors.begin();
    end = curNeighbors.end();
    std::random_shuffle(start, end);

    for (unsigned int i = 0; i < curNeighbors.size(); i++) {
      curNeighbor = curNeighbors[i];
      //std::cout << " curNeighbor = " << vertexes[curNeighbor]->getFileID() << std::endl;
      int ec = 0;
      for (unsigned int j = 0; j < edges.size(); j++) {
        pBond = edges[j];
        if (((pBond->atom1 == vertexes[curAtom] and pBond->atom2 == vertexes[curNeighbor]) or
             (pBond->atom2 == vertexes[curAtom] and pBond->atom1 == vertexes[curNeighbor])) and !loop ) {
          curEdge = j;
          //std::cout << " edgesColor = " << edgesColor[j] << std::endl;
          ec = edgesColor[j];
          edgesColor[j] = 1; // black
        }
        if ((pBond->atom1 == vertexes[rootAtom] and pBond->atom2 == vertexes[curNeighbor]) or
            (pBond->atom2 == vertexes[rootAtom] and pBond->atom1 == vertexes[curNeighbor]) ) {
           bondIndex = j;
        }
      }
      if (vertexesColor[curNeighbor] == 0 and !loop) {
        path.push_back(vertexes[curNeighbor]);
        //std::cout << " pushing back1 = "
        //          << vertexes[rootAtom]->getFileID() << " "
        //          << vertexes[curNeighbor]->getFileID() << " "
        //          << pParent->hasBond(vertexes[rootAtom],vertexes[curNeighbor])
        //          << std::endl;

          if (!pParent->hasBond(vertexes[rootAtom],vertexes[curNeighbor])) {
            this->dfs_visitNEW(curNeighbor,rootAtom,first,loop,vertexes,edges,neighbors,vertexesColor,edgesColor,path);
          }
          else {
            edgesColor[bondIndex] = 1;
            vertexesColor[curNeighbor] = 1; // black
          }
      }

      if (ec == 0 and vertexesColor[curNeighbor] == 1 and !loop) {
        if (first) {
          //std::cout << " rootAtom = " << vertexes[rootAtom]->getFileID() << std::endl;
          //std::cout << " curAtom = " << vertexes[curAtom]->getFileID() << std::endl;
          //std::cout << " curNeighbor = " << vertexes[curNeighbor]->getFileID() << std::endl;
          if (pParent->hasBond(vertexes[rootAtom],vertexes[curNeighbor]) and edgesColor[bondIndex] == 1) {
            edgesColor[bondIndex] = 1;
            loop = true;
          }
          else {
            edgesColor[curEdge] = 0;
          }
        }
        else {
          loop = true;
          atomIterator = std::find(path.begin(), path.end(), vertexes[curNeighbor]);
          if (atomIterator == path.end()) {// and vertexesColor[curNeighbor] == 0) {
            path.push_back(vertexes[curNeighbor]);
            //std::cout << " pushing back2 = " << vertexes[curNeighbor]->getFileID() << std::endl;
          }
        }
      }
    }
}

// ============================================================
// Function : kekulize()
// ------------------------------------------------------------
// Determine if ring is aromatic or not
// ============================================================
void rings::kekulize(ring* r)
{
    //std::cout << "\n\n rings::kekulize " << std::endl;
    int ringSize = r->atoms.size();
    std::vector<atom*> bondedAtoms;
    double dTorsion = 0.0;
    int planar = 0;
    double largestTorsion = 0.0;
    int huckel = 0;
    std::string symbol = "";
    int formalCharge = 0;

    // Check if planar
    for (int i = 0; i < ringSize-3; i++) {
      dTorsion = torsion(*(r->atoms[i]->getCoords()),
                         *(r->atoms[i+1]->getCoords()),
                         *(r->atoms[i+2]->getCoords()),
                         *(r->atoms[i+3]->getCoords()));

      // make all torsion positive (from 0 to 360 degrees)
      if (dTorsion < 0) dTorsion = dTorsion + 2*PI;
      if (dTorsion > largestTorsion) largestTorsion = dTorsion;

      // 15 degree cutoff (was 10 degrees, change 12 to 18 below)
      if ((std::abs(dTorsion - 0)    < PI/12) or
          (std::abs(dTorsion - PI)   < PI/12) or
          (std::abs(dTorsion - 2*PI) < PI/12) ) {
/*
        std::cout << r->atoms[i]->getFileID()   << "-"
                  << r->atoms[i+1]->getFileID() << "-"
                  << r->atoms[i+2]->getFileID() << "-"
                  << r->atoms[i+3]->getFileID() << " :"
                  << dTorsion * RAD2DEG << std::endl;
*/
        planar = 1;
        r->planar = 1;
      }
      else {
        planar = 0;
        r->planar = 0;
        break;
      }
    }
    //std::cout << "   ring::kekulize, largest ring torsion [15 degree cut-off] =  "
    //          << largestTorsion * RAD2DEG << std::endl;

    // check to see if ring is nonaromatic
    //  1) no double bonds
    int nonaromatic1 = 1;
    for (int i = 0; i < ringSize; i++) {
      if (i != ringSize-1) {
        pBond = pParent->getBond(r->atoms[i], r->atoms[i+1]);
      }
      else {
        pBond = pParent->getBond(r->atoms[i], r->atoms[0]);
      }
      if (pBond) {
        // std::cout << "bond " << pBond->atom1->getHybridization() << "-" << pBond->atom2->getHybridization() << std::endl;
        if (pBond->atom1->getHybridization() == 3 and pBond->atom2->getHybridization() == 3) {
//        if ((pBond->type == 2) or (pBond->type == 4)) {
          nonaromatic1 = 0;
        }
      }
    }

    //  2) contains quaternary atom (bonded to 4 other heavy atoms)
    int nonaromatic2 = 0;
    int nSaturatedCarbons = 0;

    for (int i = 0; i < ringSize; i++) {
      bondedAtoms = r->atoms[i]->getBondedAtoms();
      if (bondedAtoms.size() == 4) {
        bool hFound = false;
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          if (bondedAtoms[j]->getElement()->symbol == "H") hFound = true;
          if (bondedAtoms[j]->getElement()->symbol == "C") nSaturatedCarbons++;
        }
          if (!hFound) nonaromatic2 = 1;
      }
    }

    //  3) contains more than one saturated carbon
    int nonaromatic3 = 0;
    if (nSaturatedCarbons > 1) {
      nonaromatic3 = 1;
    }

    //  4) contains a monoradical
    int nonaromatic4 = 0;

///////////////// THIS CODE DOESN'T WORK ////////////////////

    //  5) contains sulfoxide or sulfone
    /*
         check for:
              sulfoxide  or  sulfone
                  O          O   O
                 ||          \\ //
                 S             S
                / \           / \
         if true then ring cannot be aromatic
    */
    int nonaromatic5 = 0;
    bool ringMember = false;
    for (int i = 0; i < ringSize; i++) {
      if (r->atoms[i]->getElement()->symbol == "S") {
        bondedAtoms = r->atoms[i]->getBondedAtoms();
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          for (int x = 0; x < ringSize; x++) {
            if (r->atoms[x] == bondedAtoms[j]) {
              ringMember = true;
            }
          }
          if (!ringMember) {
            pBond = pParent->getBond(r->atoms[i], bondedAtoms[j]);
          }
          ringMember = false;
        }
      }
    }
//////////////////

    int nonaromatic = nonaromatic1 + nonaromatic2 + nonaromatic3 + nonaromatic4 + nonaromatic5;
/*
    std::cout << "   ring::kekulize: nonaromatic1 = " << nonaromatic1 << std::endl;
    std::cout << "   ring::kekulize: nonaromatic2 = " << nonaromatic2 << std::endl;
    std::cout << "   ring::kekulize: nonaromatic3 = " << nonaromatic3 << std::endl;
    std::cout << "   ring::kekulize: nonaromatic5 = " << nonaromatic5 << std::endl;
    std::cout << "   ring::kekulize: nonaromatic : " << nonaromatic << std::endl;
    std::cout << "   ring::kekulize: planar : " << planar << std::endl;
    std::cout << "   ring::kekulize: ringSize : " << ringSize << std::endl;
*/

    // If its already aromatic, then there is no need to test for aromaticity
    bool doHuckelTest = true;
    int numAro = 0;
    for (int i = 0; i < ringSize; i++) {
      if (i != ringSize-1) {
        pBond = pParent->getBond(r->atoms[i], r->atoms[i+1]);
        if (pBond) {
          if (pBond->type == 4) numAro++;
        }
      }
      else {
        pBond = pParent->getBond(r->atoms[i], r->atoms[0]);
        if (pBond) {
          if (pBond->type == 4) numAro++;
        }
      }
    }
    if (numAro == ringSize) {
      doHuckelTest = false;
      huckel = 1;
    }


    // if planar, count the number of pi electrons
    if ( (planar) and (!nonaromatic) and (ringSize > 4) and (doHuckelTest)) {
      int huckelNumbers[8] = {2,6,10,14,18,22,26,30};
      //int nAromatic = 2 + 4 * ( (int) ( (double)(ringSize-2) * 0.25 + 0.5) );
      int numberPiElectrons = 0;
      for (int i = 0; i < ringSize; i++) {
        symbol = r->atoms[i]->getElement()->symbol;
        formalCharge = r->atoms[i]->getFormalCharge();

        if (symbol == "C") {
          if (formalCharge == 1) { // cationic
            // add zero to numberPiElectrons
          }
          else if (formalCharge == -1) { // anionic
            numberPiElectrons += 2;
          }
          else { // neutral
            numberPiElectrons++;
          }
        }

        if (symbol == "O") {
          numberPiElectrons += 2;
        }

        if (symbol == "S") {
          numberPiElectrons += 2;
        }

        if (symbol == "N") {
          Bond* pBond1 = 0;
          Bond* pBond2 = 0;
          if (i == 0) {
            pBond1 = pParent->getBond(r->atoms[i], r->atoms[ringSize-1]);
          }
          else {
            pBond1 = pParent->getBond(r->atoms[i], r->atoms[i-1]);
          }
          if (i == ringSize-1) {
            pBond2 = pParent->getBond(r->atoms[i], r->atoms[0]);
          }
          else {
            pBond2 = pParent->getBond(r->atoms[i], r->atoms[i+1]);
          }
          if (pBond1 and pBond2) {
            // non-basic, lone pair of electrons delocalized into the pi system, e.g. pyrrole
            if ( (pBond1->type == 1) and (pBond2->type == 1) ) {
              numberPiElectrons += 2;
            }
            // basic, lone pair doesn't contribute to pi system, e.g. imidazole
            else {
              numberPiElectrons++;
            }
          }
        }
//#ifdef DEBUG
//        std::cout << "   ring::kekulize: numberPiElectrons = " << numberPiElectrons << std::endl;
//#endif
      }

////// TEST CODE
      // Check for exocyclic pi bonds
      for (int i = 0; i < ringSize; i++) {
        symbol = r->atoms[i]->getElement()->symbol;
        bondedAtoms = r->atoms[i]->getBondedAtoms();
        for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
          for (int x = 0; x < ringSize; x++) {
            if (r->atoms[x] == bondedAtoms[j]) {
              ringMember = true;
            }
          }
          if (!ringMember) {
            pBond = pParent->getBond(r->atoms[i], bondedAtoms[j]);
            if ((pBond->type == 2) and (symbol == "C")) {
              std::string exoAtom = bondedAtoms[j]->getElement()->symbol;
              if ((exoAtom == "O") or (exoAtom == "S")) {
                numberPiElectrons--;
              }
            }
          }
          ringMember = false;
        }
      }
////////
      for (int j = 0; j < 8; j++) {
        if (numberPiElectrons == huckelNumbers[j]) {
          huckel = 1;
        }
      }
    }

//
//    std::cout << "   huckel = " << huckel << std::endl;
//

    // if number of pi electrons = 4n+2 => 2,6,10,14,... then its aromatic
    if (planar and huckel) {
      r->aromatic = 1;
      for (int i = 0; i < ringSize; i++) {
        if (i != ringSize-1) {
          pBond = pParent->getBond(r->atoms[i], r->atoms[i+1]);
        }
        else {
          pBond = pParent->getBond(r->atoms[i], r->atoms[0]);
        }

        if (pBond) {
          if (pBond->type == 1) {
            pBond->type = 6; // Single or Aromatic
          }
          else if (pBond->type == 2) {
            pBond->type = 7; // Double or Aromatic
          }
          pBond->topology = 1; // ring
        }
        else {
          errorLogger.throwError("rings::kekulize", " Can't find bond", 1);
          //exit(0);
          std::stringstream ss;
          ss<<"rings::kekulize"<< " Can't find bond";
                  throw MTKException(ss.str());
        }
        r->atoms[i]->setType(6); // Aromatic Ring Heavy Atom
        r->atoms[i]->setHybridization(3); // sp2
      }
    }

    if (planar and !huckel) {
      // if ring contains hetro atoms then make that atom sp2
      for (int i = 0; i < ringSize; i++) {
        symbol = r->atoms[i]->getElement()->symbol;
        if (symbol != "C") {
          r->atoms[i]->setHybridization(3);
        }
        if (r->atoms[i]->getType() != 6) {
          r->atoms[i]->setType(5);
        }
      }
      r->aromatic = 0;
    }

    if (!planar and !huckel) {
      for (int i = 0; i < ringSize; i++) {
        if (r->atoms[i]->getType() < 5) {
          r->atoms[i]->setType(5);
        }
      }
      r->aromatic = 0;
    }
    pParent->bBondTypes2Assigned = true;
}

// ============================================================
// Function : calcCentroid()
// ------------------------------------------------------------
// Determine the center of the ring
// ============================================================
void rings::calcCentroid(ring* r)
{
    unsigned int ringSize = r->atoms.size();
    r->centroid.resize(3);

    for (unsigned int i = 0; i < r->centroid.size(); i++) {
      r->centroid(i) = 0.0;
    }

    for (unsigned int i = 0; i < ringSize; i++) {
      atom* at = r->atoms[i];
      vector3d* coord = at->getCoords();
      for (unsigned int j = 0; j < 3; j++) {
        r->centroid(j) += (*coord)[j];
      }
    }
    for (unsigned int i = 0; i < r->centroid.size(); i++) {
      r->centroid(i) /= double(ringSize);
    }
}

// ============================================================
// Function : getPlaneNormal()
// ------------------------------------------------------------
// Determine the plane and normal of the ring
// ============================================================
int rings::getPlaneNormal(ring* r)
{
    this->calcCentroid(r);

    unsigned int ringSize = r->atoms.size();

    // Boost
    // ublas::matrix<double, ublas::column_major> coordMatrix(ringSize,3); // rows, columns
    // ublas::vector<double> R(3);

    // Eigen 
    Eigen::MatrixXd coordMatrix(ringSize, 3);

    r->planeNormal.resize(3,3);

    // size1 == rows
    // size2 == columns
    //int s1 = coordMatrix.size1();
    //int s2 = coordMatrix.size2();
    unsigned int s1 = coordMatrix.rows();
    unsigned int s2 = coordMatrix.cols();

    for (unsigned int i = 0; i < s1; i++) {
      for (unsigned int j = 0; j < s2; j++) {
        coordMatrix(i,j) = 0.0;
      }
    }

    for (unsigned int i = 0; i < 3; i++) {
      //R(i) = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        r->planeNormal(i,j) = 0.0;
      }
    }

    for (unsigned int i = 0; i < s1; i++) { // rows
      atom* at = r->atoms[i];
      vector3d* coord = at->getCoords();
      for (unsigned int j = 0; j < s2; j++) { // columns
        coordMatrix(i,j) = (*coord)[j] - r->centroid(j);;
      }
    }

    r->planeNormal = coordMatrix.transpose() * coordMatrix;

    SelfAdjointEigenSolver<MatrixXd> eigensolver(r->planeNormal);
    MatrixXd evectors = eigensolver.eigenvectors();
    VectorXd R = eigensolver.eigenvalues();

    eigenValueSort(evectors, R, 1);

/*
    r->planeNormal = ublas::prod(ublas::trans(coordMatrix),coordMatrix);
    int result = diagonalize(r->planeNormal,R);
    if (result != 0) return 1;

    eigenValueSort(r->planeNormal,R,1);
*/

#ifdef DEBUG
    std::string errMessage = "\n RING: ";
    for (unsigned int i = 0; i < coordMatrix.size1(); i++) { // rows
      atom* at = r->atoms[i];
      errMessage += i2s(at->getIndex()) + " ";
    }
    errMessage += "\n CENTROID Centroid_" + pParent->getName() + "_" + i2s(r->index);
    for (unsigned int i = 0; i < r->centroid.size(); i++) {
      errMessage += " " + d2s(r->centroid(i));
    }
    errMessage +=" \n";

/*
   printing using mat(j,i):
   eigenvectors:
    EV1_x EV1_y EV1_z
    EV2_x EV2_y EV2_z
    EV3_x EV3_z EV3_z

   eigenvalues:
    V1 V2 V3
*/
    for (unsigned int i = 0; i < 2; i++) {
      errMessage += " PLANE P_" + pParent->getName() + "_" + i2s(r->index) + "_" + i2s(i+1);
      for (unsigned int j = 0; j < 3; j++) {
        errMessage += " " + d2s(r->planeNormal(j,i) + r->centroid(j));
      }
      errMessage += " \n";
    }

    errMessage += " PLANE N_" + pParent->getName() + "_" + i2s(r->index);
    for (unsigned int j = 0; j < 3; j++) {
      errMessage += " " + d2s(r->planeNormal(j,2) + r->centroid(j));
    }
    errMessage += " ";
    errorLogger.throwError("rings::getPlaneNormal", errMessage, 4);
#endif

    return 0;
}

} // MTKpp namespace

