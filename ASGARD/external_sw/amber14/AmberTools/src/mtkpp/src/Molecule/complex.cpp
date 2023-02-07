/*!
   \file complex.cpp
   \brief Determine the receptor, solvent and ligands
   \author Martin B. Peters

   $Date: 2010/03/29 20:42:28 $
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
#include "complex.h"
#include "collection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "functionalize.h"

#include "stdLibrary.h"
#include "stdGroup.h"
#include "stdFrag.h"

#include "Utils/vector3d.h"
#include "element.h"

// Graph
#include "Graph/graph.h"
#include "Graph/edge.h"
#include "Graph/vertex.h"

#include "Log/errorHandler.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : complex()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
complex::complex(collection* pCollection, bool l)
{
    this->cofactors.push_back("SO4");

    this->pCollection = pCollection;

    this->rec = 0;
    this->lig = 0;
    this->cof = 0;
    this->pLigand = 0;
    this->sol = 0;
    this->recSol = 0;
    this->atomIndices = 0;

    this->resRec = 0;
    this->resLig = 0;
    this->resSol = 0;
    this->resRecSol = 0;
    this->resIndices = 0;
    this->chain = 0;
    this->chainCounter = 0;

    this->bLigOnly = l;

    this->setup();
}

// ============================================================
// Function : ~complex()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
complex::~complex()
{
    delete [] this->rec;
    delete [] this->lig;
    delete [] this->sol;
    delete [] this->recSol;
    delete [] this->atomIndices;
    delete [] this->chain;

    delete [] this->resRec;
    delete [] this->resLig;
    delete [] this->resSol;
    delete [] this->resRecSol;
    delete [] this->resIndices;
}

// ============================================================
// Function : setup()
// ------------------------------------------------------------
// Define receptor, ligand, and solvent molecules
// ============================================================
void complex::setup()
{
    this->nAtoms = pCollection->getNumAtoms();
    this->nResidues = pCollection->getNumberSubMolecules();

    this->nAtoms = nAtoms;
    try {
      this->rec         = new int [this->nAtoms];
      this->lig         = new int [this->nAtoms];
      this->cof         = new int [this->nAtoms];
      this->sol         = new int [this->nAtoms];
      this->recSol      = new int [this->nAtoms];
      this->atomIndices = new int [this->nAtoms];
      this->chain       = new int [this->nAtoms];

      this->resRec      = new int [this->nResidues];
      this->resLig      = new int [this->nResidues];
      this->resSol      = new int [this->nResidues];
      this->resRecSol   = new int [this->nResidues];
      this->resIndices  = new int [this->nResidues];
    }
    catch (std::bad_alloc) {
      MTKpp::errorLogger.throwError("complex::setup", " Memory Allocation Failure ", 1);
      //exit(1);
      std::stringstream ss;
      ss << "complex::setup"<< " Memory Allocation Failure ";
      throw MTKException(ss.str());
    }

    for (int a = 0; a < this->nAtoms; a++) {
      this->rec[a] = 0;
      this->lig[a] = 0;
      this->cof[a] = 0;
      this->sol[a] = 0;
      this->recSol[a] = 0;
      this->atomIndices[a] = 0;
      this->chain[a] = 0;
    }

    for (int a = 0; a < this->nResidues; a++) {
      this->resRec[a] = 0;
      this->resLig[a] = 0;
      this->resSol[a] = 0;
      this->resRecSol[a] = 0;
      this->resIndices[a] = 0;
    }

    std::vector<molecule*> molList = pCollection->getMoleculeList();
    std::vector<atom*> colAtomList = pCollection->getAtomList();

    std::vector<double> molWeights;
    for (unsigned int m = 0; m < molList.size(); m++) {
      double mw = molList[m]->getMolecularWeight();
      if (mw < 0.5) {
        MTKpp::errorLogger.throwError("complex::setup", " Please assign elements ", 1);
      }
      molWeights.push_back(mw);
    }

    this->nRecAtoms = 0;
    this->nLigAtoms = 0;
    int atomIndex = 0;
    int resIndex = 0;

    std::vector<int> lChains;

    for (unsigned int m = 0; m < molList.size(); m++) {
      std::vector<submolecule*> sMolList = molList[m]->getSubMoleculeList();
      if (molWeights[m] > 1000.0) {
        this->chainCounter++;
        lChains.push_back(0);
        for (unsigned int s = 0; s < sMolList.size(); s++) {
          receptor.push_back(sMolList[s]);
          this->nRecAtoms += static_cast<int>(sMolList[s]->getNumAtoms());

          std::vector<atom*> atomList = sMolList[s]->getAtomList();
          for (unsigned int a = 0; a < atomList.size(); a++) {
            this->rec[atomIndex] = 1;
            this->chain[atomIndex] = this->chainCounter;
            atomIndex++;
          }
          this->resRec[resIndex] = 1;
          resIndex++;
        }
      }
      else {
        for (unsigned int s = 0; s < sMolList.size(); s++) {
          std::string resName = sMolList[s]->getName();

          bool bCofactor = false;
          for (unsigned r = 0; r < cofactors.size(); r++) {
            if (resName == cofactors[r]) {
              bCofactor = true;
            }
          }

          if (bCofactor) {
            this->cofactor.push_back(sMolList[s]);
            this->nRecAtoms += static_cast<int>(sMolList[s]->getNumAtoms());

            std::vector<atom*> atomList = sMolList[s]->getAtomList();
            for (unsigned int a = 0; a < atomList.size(); a++) {
              this->cof[atomIndex] = 1;
              atomIndex++;
            }
            this->resSol[resIndex] = 1;
            resIndex++;
          }
          else if (resName == "HOH" or resName == "WAT") {
            solvent.push_back(sMolList[s]);
            this->nRecAtoms += static_cast<int>(sMolList[s]->getNumAtoms());

            std::vector<atom*> atomList = sMolList[s]->getAtomList();
            for (unsigned int a = 0; a < atomList.size(); a++) {
              this->sol[atomIndex] = 1;
              atomIndex++;
            }
            this->resSol[resIndex] = 1;
            resIndex++;
          }
          else {
            std::vector<atom*> atomList = sMolList[s]->getAtomList();
            if (atomList.size() == 1) {
              if (atomList[0]->getElement()->number > 18) {
                receptor.push_back(sMolList[s]);
                this->nRecAtoms += static_cast<int>(sMolList[s]->getNumAtoms());

                for (unsigned int a = 0; a < atomList.size(); a++) {
                  this->rec[atomIndex] = 1;
                  atomIndex++;
                }
                this->resRec[resIndex] = 1;
                resIndex++;
              }
            }
            else {
              ligand.push_back(sMolList[s]);
              pLigand = sMolList[s]->getParent();
              this->nLigAtoms += static_cast<int>(sMolList[s]->getNumAtoms());

              for (unsigned int a = 0; a < atomList.size(); a++) {
                this->lig[atomIndex] = 1;
                atomIndex++;
              }
              this->resLig[resIndex] = 1;
              resIndex++;
            }
          }
        }
      }
    }

    // Define chains
    atomIndex = 0;
    int atomIndex2 = 0;
    atom* pAtom1 = 0;
    atom* pAtom2 = 0;
    vector3d Coord1;
    vector3d Coord2;
    double lDist = 0.0;

    for (unsigned int m = 0; m < molList.size(); m++) {
      atomIndex2 = atomIndex;
      std::vector<submolecule*> sMolList = molList[m]->getSubMoleculeList();
      if (molWeights[m] < 1000.0) {

        for (int c = 0; c < this->chainCounter; c++) {
          lChains[c] = 0;
        }
        int nlAts = 0;

        for (unsigned int s = 0; s < sMolList.size(); s++) {
          std::string resName = sMolList[s]->getName();
          if (resName != "HOH" and resName != "WAT") {
            std::vector<atom*> atomList = sMolList[s]->getAtomList();
            for (unsigned int a = 0; a < atomList.size(); a++) {
              nlAts++;
              pAtom1 = atomList[a];
              Coord1 = (*pAtom1->getCoords());
              double minDist = 100.0;
              int closestAtom = 0;
              for (unsigned int c = 0; c < colAtomList.size(); c++) {
                if (this->chain[c] > 0) {
                  pAtom2 = colAtomList[c];
                  Coord2 = (*pAtom2->getCoords());
                  lDist = Coord1.dist(Coord2);
                  if (lDist < minDist) {
                    minDist = lDist;
                    closestAtom = c;
                  }
                }
              }
              this->chain[atomIndex] = -this->chain[closestAtom];
              lChains[this->chain[closestAtom]-1]++;
              atomIndex++;
            }
          }
          else { // waters
            std::vector<atom*> atomList = sMolList[s]->getAtomList();
            for (unsigned int a = 0; a < atomList.size(); a++) {
              atomIndex++;
            }
          }
        }

        // majority
        if (nlAts > 0) {
          for (int c = 0; c < this->chainCounter; c++) {
            double p = (double(lChains[c])/double(nlAts)) * 100.0;
            if (p > 99.9) {
              break;
            }
            else if (p > 50.0) {
              for (unsigned int s = 0; s < sMolList.size(); s++) {
                std::string resName = sMolList[s]->getName();
                if (resName != "HOH" and resName != "WAT") {
                  std::vector<atom*> atomList = sMolList[s]->getAtomList();
                  for (unsigned int a = 0; a < atomList.size(); a++) {
                    this->chain[atomIndex2] = -(c+1);
                    atomIndex2++;
                  }
                }
              }
            }
          }
        }

      }
      else { // protein
        for (unsigned int s = 0; s < sMolList.size(); s++) {
          std::vector<atom*> atomList = sMolList[s]->getAtomList();
          for (unsigned int a = 0; a < atomList.size(); a++) {
            atomIndex++;
          }
        }
      }
    }

    // set recSol and atomIndices
    for (int a = 0; a < this->nAtoms; a++) {
      if (rec[a] or sol[a]) {
        this->recSol[a] = 1;
      }
    }

    int p = 0;
    for (int a = 0; a < this->nAtoms; a++) {
      if (recSol[a] == 1) {
        this->atomIndices[a] = p;
        p++;
      }
    }

    for (int a = 0; a < this->nAtoms; a++) {
      if (recSol[a] == 0) {
        this->atomIndices[a] = p;
        p++;
      }
    }

    // Set resRecSol and resIndices
    for (int a = 0; a < this->nResidues; a++) {
      if (resRec[a] == 1 or resSol[a] == 1) {
        this->resRecSol[a] = 1;
      }
    }

    p = 0;
    for (int a = 0; a < this->nResidues; a++) {
      if (resRecSol[a] == 1) {
        this->resIndices[a] = p;
        p++;
      }
    }

    for (int a = 0; a < this->nResidues; a++) {
      if (resRecSol[a] == 0) {
        this->resIndices[a] = p;
        p++;
      }
    }

    std::string recStr = "";
    // Checks
    if (receptor.size() == 0) {
      MTKpp::errorLogger.throwError("complex::setup", " No Receptor Found ", 1);
      //exit(1);
      std::stringstream ss;
      ss << "complex::setup"<< " No Receptor Found ";
      throw MTKException(ss.str());
    }
    else {
      recStr += " RECEPTOR: ";
      for (unsigned int s = 0; s < receptor.size(); s++) {
        recStr += receptor[s]->getName() + " ";
      }
      MTKpp::errorLogger.throwError("complex::setup", recStr, 4);
    }

    std::string ligStr = "";
    if (ligand.size() == 0) {
      MTKpp::errorLogger.throwError("complex::setup", " No Ligand Found ", 1);
      //exit(0);
      std::stringstream ss;
      ss << "complex::setup"<< " No Ligand Found ";
      throw MTKException(ss.str());
    }
    else {
      ligStr += " LIGAND: ";
      for (unsigned int s = 0; s < ligand.size(); s++) {
        ligStr += ligand[s]->getName() + " ";
      }
      MTKpp::errorLogger.throwError("complex::setup", ligStr, 4);
    }

    std::string solStr = "";
    if (solvent.size() == 0) {
      MTKpp::errorLogger.throwError("complex::setup", " No Solvent Found ", 4);
    }
    else {
      solStr += " SOLVENT: ";
      for (unsigned int s = 0; s < solvent.size(); s++) {
        solStr += solvent[s]->getName() + " ";
      }
      MTKpp::errorLogger.throwError("complex::setup", solStr, 4);
    }

    this->nRecResidues = static_cast<int>(receptor.size() + solvent.size());

    this->nRecRecAtoms = (this->nRecAtoms * (this->nRecAtoms+1))/2;
    this->nRecRecResidues  = (this->nRecResidues * (this->nRecResidues+1))/2;

    if (this->bLigOnly) {
      this->atomMatrixSize = (this->nAtoms * (this->nAtoms+1))/2 - this->nRecRecAtoms;
      this->resMatrixSize  = (this->nResidues * (this->nResidues+1))/2 - this->nRecRecResidues;
    }
    else {
      this->atomMatrixSize = (this->nAtoms * (this->nAtoms+1))/2;
      this->resMatrixSize  = (this->nResidues * (this->nResidues+1))/2;
    }
}

// ============================================================
// Function : getNumAtoms()
// ------------------------------------------------------------
// Get the number of atoms
// ============================================================
int complex::getNumAtoms()
{
    return this->nAtoms;
}

// ============================================================
// Function : getNumResidues()
// ------------------------------------------------------------
// Get the number of residues
// ============================================================
int complex::getNumResidues()
{
    return this->nResidues;
}

// ============================================================
// Function : getNumRecAtoms()
// ------------------------------------------------------------
// Get the number of receptor atoms
// ============================================================
int complex::getNumRecAtoms()
{
    return this->nRecAtoms;
}

// ============================================================
// Function : getNumRecResidues()
// ------------------------------------------------------------
// Get the number of receptor residues
// ============================================================
int complex::getNumRecResidues()
{
    return this->nRecResidues;
}

// ============================================================
// Function : getNumRecResidues()
// ------------------------------------------------------------
// Get the number of ligand atoms
// ============================================================
int complex::getNumLigAtoms()
{
    return this->nLigAtoms;
}

// ============================================================
// Function : getNumLigResidues()
// ------------------------------------------------------------
// Get the number of ligand residues
// ============================================================
int complex::getNumLigResidues()
{
    return this->nLigResidues;
}

// ============================================================
// Function : getAtomsMatrixSize()
// ------------------------------------------------------------
//
// ============================================================
int complex::getAtomsMatrixSize()
{
    return this->atomMatrixSize;
}

// ============================================================
// Function : getResMatrixSize()
// ------------------------------------------------------------
//
// ============================================================
int complex::getResMatrixSize()
{
    return this->resMatrixSize;
}

// ============================================================
// Function : getLigOnly()
// ------------------------------------------------------------
//
// ============================================================
bool complex::getLigOnly()
{
    return this->bLigOnly;
}

// ============================================================
// Function : getRecRecAtom()
// ------------------------------------------------------------
//
// ============================================================
int complex::getRecRecAtoms()
{
    return this->nRecRecAtoms;
}

// ============================================================
// Function : getRecRecResidues()
// ------------------------------------------------------------
//
// ============================================================
int complex::getRecRecResidues()
{
    return this->nRecRecResidues;
}

// ============================================================
// Function : getRecFlags()
// ------------------------------------------------------------
//
// ============================================================
int complex::getRecFlags(int recFlags[])
{
    try {
      for (int i = 0; i < nAtoms; i++) {
        recFlags[i] = rec[i];
      }
    }
    catch (std::bad_alloc) {
      MTKpp::errorLogger.throwError("complex::getRecFlags", " Memory Out of bounds Failure ", 1);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getLigFlags()
// ------------------------------------------------------------
//
// ============================================================
int complex::getLigFlags(int ligFlags[])
{
    try {
      for (int i = 0; i < nAtoms; i++) {
        ligFlags[i] = lig[i];
      }
    }
    catch (std::bad_alloc) {
      MTKpp::errorLogger.throwError("complex::getLigFlags", " Memory Out of bounds Failure ", 1);
      return 1;
    }
    return 0;
}

// ============================================================
// Function : isAtomLig()
// ------------------------------------------------------------
//
// ============================================================
bool complex::isAtomLig(int a)
{
    if (a > this->nAtoms-1) {
      return false;
    }

    if (this->lig[a]) {
      return true;
    }
    return false;
}

// ============================================================
// Function : getAtomIndex()
// ------------------------------------------------------------
//
// ============================================================
int complex::getAtomIndex(int a)
{
    if (a > this->nAtoms-1) {
      return -1;
    }
    return atomIndices[a];
}

// ============================================================
// Function : isResLig()
// ------------------------------------------------------------
//
// ============================================================
bool complex::isResLig(int r)
{
    if (r > this->nResidues-1) {
      return false;
    }

    if (this->resLig[r]) {
      return true;
    }

    return false;
}

// ============================================================
// Function : isChain()
// ------------------------------------------------------------
//
// ============================================================
bool complex::isChain(int c, int a)
{
    if (c > this->nAtoms-1) {
      return false;
    }

    if (std::abs(this->chain[a]) == c) {
      return true;
    }

    return false;
}

// ============================================================
// Function : getResIndex()
// ------------------------------------------------------------
//
// ============================================================
int complex::getResIndex(int r)
{
    if (r > this->nResidues-1) {
      return -1;
    }
    return resIndices[r];
}

// ============================================================
// Function : setCoreFragment()
// ------------------------------------------------------------
//
// ============================================================
int complex::setCoreFragment(std::string core)
{
    if (pLigand == 0) {
      return 1;
    }

    int resIndex = 1;

    // Add receptor residues
    for (unsigned int i = 0; i < receptor.size(); i++) {
      std::vector<atom*> a = receptor[i]->getAtomList();
      std::vector<int> va;
      for (unsigned j = 0; j < a.size(); j++) {
        va.push_back(a[j]->getColIndex());
      }
      residues.push_back(va);

      // Set residue map
      std::string rName = receptor[i]->getName();
      std::stringstream rId;
      rId << receptor[i]->getIndex();
      rName+=rId.str();
      this->addMapping(rName, rName);
      residueNAMEIDs.push_back(rName);
      resIndex++;
    }

    // Add solvent residues
    for (unsigned int i = 0; i < solvent.size(); i++) {
      std::vector<atom*> a = solvent[i]->getAtomList();
      std::vector<int> va;
      for (unsigned j = 0; j < a.size(); j++) {
        va.push_back(a[j]->getColIndex());
      }
      residues.push_back(va);

      // Set residue map
      std::string rName = solvent[i]->getName();
      std::stringstream rId;
      rId << solvent[i]->getIndex();
      rName+=rId.str();
      this->addMapping(rName, rName);
      residueNAMEIDs.push_back(rName);

      resIndex++;
    }

    // Add ligand residues

    // Check to see if user doesn't want to fragment the ligand
    if (core == " ") {
      for (unsigned int i = 0; i < ligand.size(); i++) {
        std::vector<atom*> a = ligand[i]->getAtomList();
        std::vector<int> va;
        for (unsigned j = 0; j < a.size(); j++) {
          va.push_back(a[j]->getColIndex());
        }
        residues.push_back(va);

        // Set residue map
        std::string rName = ligand[i]->getName();
        std::stringstream rId;
        rId << ligand[i]->getIndex();
        rName+=rId.str();
        this->addMapping(rName, rName);
        residueNAMEIDs.push_back(rName);
        ligandIndices.push_back(resIndex);
        resIndex++;
      }
      return 0;
    }

    funcGroup* coreFuncGroup = 0;
    stdFrag* coreStdFrag = 0;
    std::vector<atom*> coreAtoms;
    typedef std::map<stdAtom*, atom*>::iterator funcGroupMapIterator;
    std::vector<funcGroup*> funcGroups = pLigand->getFunctionalGroups();
    std::vector<atom*> ligAtoms = pLigand->getAtomList();
    for (unsigned int w = 0; w < funcGroups.size(); w++) {
      if (funcGroups[w]->pStdFrag->getSymbol() == core) {
        std::vector<int> va;
        coreFuncGroup = funcGroups[w];
        coreStdFrag = funcGroups[w]->pStdFrag;

        for (funcGroupMapIterator f = funcGroups[w]->atomMap.begin();
                                  f != funcGroups[w]->atomMap.end(); f++) {
          atom* pFGAtom = f->second;
          coreAtoms.push_back(pFGAtom);
          va.push_back(pFGAtom->getColIndex());
        }
        residues.push_back(va);

        // Set residue map
        std::string rName = core;
        std::stringstream rId;
        rId << resIndex;
        rName+=rId.str();
        this->addMapping(rName, rName);
        residueNAMEIDs.push_back(rName);
        ligandIndices.push_back(resIndex);
        resIndex++;
      }
    }

    // Create other fragments, create graph, and find subgraphs
    graph* molGraph = new graph();
    vertex* pVertexI = 0;
    vertex* pVertexJ = 0;

    for (unsigned int i = 0; i < ligAtoms.size(); i++) {
      atom* pA1 = ligAtoms[i];
      bool inCore = false;
      for (unsigned int j = 0; j < coreAtoms.size(); j++) {
        atom* pA2 = coreAtoms[j];
        if (pA1 == pA2) {
          inCore = true;
        }
      }
      if (!inCore) {
        int atI = pA1->getIndex();
        std::stringstream stAtI;
        stAtI << atI;
        pVertexI = molGraph->addVertex(atI);
        pVertexI->setName(stAtI.str().c_str());
      }
    }

    std::vector<vertex*> vertices = molGraph->getVertices();
    for (unsigned int i = 0; i < vertices.size(); i++) {
      pVertexI = vertices[i];
      atom* pA1 = pLigand->getAtom(pVertexI->getIndex());
      for (unsigned int j = i+1; j < vertices.size(); j++) {
        pVertexJ = vertices[j];
        atom* pA2 = pLigand->getAtom(pVertexJ->getIndex());
        if (pLigand->hasBond(pA1, pA2)) {
          molGraph->addEdge(pVertexI, pVertexJ);
        }
      }
    }

    // Find subgraphs
    std::vector<graph*> subGraphs;
    std::vector<vertex*> blockVertices;

    for (unsigned int v = 0; v < vertices.size(); v++) {
      if (!vertices[v]->isVisited()) {
        graph* subGraph = new graph();
        molGraph->dfs(vertices[v]);

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

    if (subGraphs.size() == 0) {
      std::cout << " No subgraphs in complex " << std::endl;
      //exit(0);
      throw MTKException(" No subgraphs in complex ");
    }

    std::vector<int> rightOrder;
    std::vector<int> connPts = coreStdFrag->getStdConnPtsList();
    for (unsigned int c = 0; c < connPts.size(); c++) {
      stdAtom* pSA = coreStdFrag->getStdAtom(connPts[c]);
      atom* pA = coreFuncGroup->atomMap[pSA];
      molecule* mol = pA->getParent()->getParent();

      for (unsigned int i = 0; i < subGraphs.size(); i++) {
        std::vector<vertex*> verticesI = subGraphs[i]->getVertices();
        for (unsigned int j = 0; j < verticesI.size(); j++) {
          pVertexI = verticesI[j];
          atom* pA2 = pLigand->getAtom(pVertexI->getIndex());

          if (pA->getParent() == pA2->getParent()) {
            if (mol->hasBond(pA, pA2)) {
              rightOrder.push_back(i);
            }
          }

        }
      }
    }
/*
    for (unsigned int k = 0; k < rightOrder.size(); k++) {
      std::cout << rightOrder[k] << " ";
    }
    std::cout << "  " << std::endl;
*/
    for (unsigned int i = 0; i < subGraphs.size(); i++) {
      int r = rightOrder[i];
      std::vector<vertex*> verticesI = subGraphs[r]->getVertices();
      std::vector<int> va;
      for (unsigned int j = 0; j < verticesI.size(); j++) {
        pVertexI = verticesI[j];
        atom* pA = pLigand->getAtom(pVertexI->getIndex());
        va.push_back(pA->getColIndex());
      }
      residues.push_back(va);

      // Set residue map
      std::string rName = "FGX";
      std::stringstream rId;
      rId << i+1;
      rName+=rId.str();
      this->addMapping(rName, rName);
      residueNAMEIDs.push_back(rName);
      ligandIndices.push_back(resIndex);

      resIndex++;
    }

    return 0;
}

// ============================================================
// Function : addMapping()
// ------------------------------------------------------------
//
// ============================================================
void complex::addMapping(std::string m1, std::string m2)
{
    this->residueMap[m1] = m2;
}

// ============================================================
// Function : getMapping()
// ------------------------------------------------------------
//
// ============================================================
std::string complex::getMapping(std::string m1)
{
    stringMapIterator b = this->residueMap.find(m1);

    if (b != this->residueMap.end()){
      return this->residueMap[m1];
    }
    return "";
}

// ============================================================
// Function : getResNameID()
// ------------------------------------------------------------
//
// ============================================================
std::string complex::getResNameID(int i)
{
    return this->residueNAMEIDs[i];
}

// ============================================================
// Function : getResidues()
// ------------------------------------------------------------
//
// ============================================================
std::vector<std::vector<int> > complex::getResidues()
{
    return this->residues;
}

// ============================================================
// Function : getLigandFragments()
// ------------------------------------------------------------
//
// ============================================================
std::vector<int> complex::getLigandFragments()
{
    return this->ligandIndices;
}

// ============================================================
// Function : getReceptorSubMols()
// ------------------------------------------------------------
//
// ============================================================
std::vector<submolecule*> complex::getReceptorSubMols()
{
    return this->receptor;
}

// ============================================================
// Function : getLigandSubMols()
// ------------------------------------------------------------
//
// ============================================================
std::vector<submolecule*> complex::getLigandSubMols()
{
    return this->ligand;
}

} // MTKpp namespace
