/*!
   \file hydrophobize.cpp
   \brief Determines the hydrophobic groups in a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:43:22 $
   $Revision: 1.4 $

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

#include "hydrophobize.h"

#include "molecule.h"
#include "atom.h"
#include "element.h"
#include "bond.h"
#include "Utils/vector3d.h"

// Graph
#include "Graph/graph.h"
#include "Graph/edge.h"
#include "Graph/vertex.h"

namespace MTKpp
{
// ============================================================
// Function : hydrophobize()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
hydrophobize::hydrophobize(molecule *parent)
{
    this->pParent = parent;

    atomList = pParent->getAtomList();
    nAtoms = atomList.size();

    bonds = pParent->getBondMap();

    atHydrophobicities.assign(nAtoms, -1);
    bondOrders.assign(nAtoms, 0);

    for (unsigned int i = 0; i < nAtoms; i++) {
      atNumbers.push_back(atomList[i]->getAtomicNum());
      atGroups.push_back(atomList[i]->getElement()->group);
      atPeriods.push_back(atomList[i]->getElement()->period);
      atSymbols.push_back(atomList[i]->getElementSymbol());
      formalCharges.push_back(atomList[i]->getFormalCharge());
      atHybridizations.push_back(atomList[i]->getHybridization());
    }

    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      std::vector<int> b;
      std::vector<atom*> bondedAtoms = pAtom->getBondedAtoms();
      int bondOrder = 0;
      for (unsigned int j = 0; j < bondedAtoms.size(); j++) {
        b.push_back(bondedAtoms[j]->getIndex()-1);
        pBond = pParent->getBond(pAtom, bondedAtoms[j]);
        if (pBond) {
          if (pBond->type > bondOrder) {
            bondOrder = pBond->type;
          }
        }
      }
      bdAtoms.push_back(b);
      bondOrders[i] = bondOrder;
    }

    OHs.assign(nAtoms, 0);
    NHs.assign(nAtoms, 0);
    SHs.assign(nAtoms, 0);

    // OH groups
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "O") {
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (atSymbols[bdAtoms[i][j]] == "H") {
            OHs[i] = 1;
          }
        }
      }
    }

    // NH groups
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "N") {
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (atSymbols[bdAtoms[i][j]] == "H") {
            OHs[i] = 1;
          }
        }
      }
    }

    // SH groups
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "S") {
        for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
          if (atSymbols[bdAtoms[i][j]] == "H") {
            OHs[i] = 1;
          }
        }
      }
    }
}

// ============================================================
// Function : ~hydrophobize()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
hydrophobize::~hydrophobize() {}

// ============================================================
// Function : run()
// ------------------------------------------------------------
//
// ============================================================
int hydrophobize::run()
{
#ifdef DEBUG
      std::cout << " hydrophobize::run " << std::endl;
#endif
    std::vector<atom*> atomList = this->pParent->getAtomList();

    // Step 1
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "O" or atSymbols[i] == "N") {
        atHydrophobicities[i] = 1;
      }
    }

    // Step 2
    for (unsigned int i = 0; i < nAtoms; i++) {
      if (atSymbols[i] == "S" and SHs[i] == 1) {
        atHydrophobicities[i] = 2;
      }
    }

    // Step 3
    for (unsigned int i = 0; i < nAtoms; i++) {
      if ((atSymbols[i] == "S") and (bondOrders[i] > 1)) {
        atHydrophobicities[i] = 3;
      }
    }

    // Step 4
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        if (formalCharges[bdAtoms[i][j]] != 0) {
          atHydrophobicities[i] = 4;
        }
      }
    }

    // Step 5
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        if (NHs[bdAtoms[i][j]] or OHs[bdAtoms[i][j]]) {
          if (atHybridizations[i] == 4) { // sp3 carbon???
            atHydrophobicities[i] = 5;
          }
        }
      }
    }

    // Step 6
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        if (SHs[bdAtoms[i][j]]) {
          if (atHybridizations[i] == 4) { // sp3 carbon???
            atHydrophobicities[i] = 6;
          }
        }
      }
    }

    // Step 7
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        if (atSymbols[bdAtoms[i][j]] == "O" and bondOrders[bdAtoms[i][j]] == 2) {
          atHydrophobicities[i] = 7;
        }
      }
    }

    // Step 8
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        if (atSymbols[bdAtoms[i][j]] == "S" and bdAtoms[bdAtoms[i][j]].size() > 2) {
          atHydrophobicities[i] = 8;
        }
      }
    }

    // Step 9
    for (unsigned int i = 0; i < nAtoms; i++) {
      int nHetero = 0;
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        if (atSymbols[bdAtoms[i][j]] == "N" or atSymbols[bdAtoms[i][j]] == "O") {
          nHetero++;
        }
      }
      if (nHetero and (atHybridizations[i] == 4)) {
        atHydrophobicities[i] = 9;
      }
    }

    // Step 10
    for (unsigned int i = 0; i < nAtoms; i++) {
      int As = 0;
      int Bs = 0;
      for (unsigned int j = 0; j < bdAtoms[i].size(); j++) {
        for (unsigned int k = 0; k < bdAtoms[  bdAtoms[i][j]   ].size(); k++) {
          if (bdAtoms[  bdAtoms[i][j]   ][k] != static_cast<int>(i)) {
            if (atSymbols[ bdAtoms[  bdAtoms[i][j]   ][k]   ] == "O") {
              if (bondOrders[  bdAtoms[  bdAtoms[i][j]   ][k]  ] > 1) {
                As++;
              }
            }
          }
        }

        if (atSymbols[bdAtoms[i][j]] == "S" and bdAtoms[ bdAtoms[i][j] ].size() > 2) {
          Bs++;
        }
      }
      if ((As + Bs) >= 2) {
        atHydrophobicities[i] = 10;
      }
    }

#ifdef DEBUG
    std::cout << "    ID   NUM GRP PRD SYM BO FC HYBD HYDRP " << std::endl;
    for (unsigned int i = 0; i < nAtoms; i++) {
      std::cout << std::setw(6) << i+1
                << std::setw(6) << atNumbers[i] << " "
                << std::setw(3) << atGroups[i] << " "
                << std::setw(3) << atPeriods[i] << " "
                << std::setw(3) << atSymbols[i] << " "
                << std::setw(2) << bondOrders[i] << " "
                << std::setw(2) << formalCharges[i] << "   "
                << std::setw(2) << atHybridizations[i] << "    "
                << std::setw(2) << atHydrophobicities[i] << std::endl;
    }
#endif

    //
    graph* molGraph = new graph();
    vertex* pVertexI = 0;
    vertex* pVertexJ = 0;

    if (!bonds.empty()) {
      for (BondMapIterator b = bonds.begin(); b != bonds.end(); b++) {
        pBond = b->second;
        int atI = pBond->atom1->getIndex();
        int atJ = pBond->atom2->getIndex();
        pVertexI = 0;
        pVertexJ = 0;

        if (atHydrophobicities[atI] < 1 and atHydrophobicities[atJ] < 1) {
          std::stringstream stAtI;
          stAtI << atI;
          std::stringstream stAtJ;
          stAtJ << atJ;

          pVertexI = molGraph->getVertex(atI);
          pVertexJ = molGraph->getVertex(atJ);

          if (!pVertexI) {
            pVertexI = molGraph->addVertex(atI);
            pVertexI->setName(stAtI.str().c_str());
          }
          if (!pVertexJ) {
            pVertexJ = molGraph->addVertex(atJ);
            pVertexJ->setName(stAtJ.str().c_str());
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

//        std::cout << "\n " << vertices[i]->getIndex() << " ";
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
      std::vector<atom*> hydph;
      for (unsigned int l = 0; l < subGV.size(); l++) {
        hydph.push_back(atomList[subGV[l]->getIndex() - 1]);
        for (unsigned int p = l; p < subGV.size(); p++) {
          if (molGraph->hasEdge(subGV[l]->getIndex(), subGV[p]->getIndex())) {
            subGraphs[k]->addEdge(subGV[l], subGV[p]);
          }
        }
      }
      pParent->addHydrophobicGroup(hydph);
    }

#ifdef DEBUG
    std::cout << "\n HYDROPHOBIC REGIONS " << std::endl;
    typedef std::map<int, edge*>::iterator edgeIterator;
    for (unsigned int k = 0; k < subGraphs.size(); k++) {
      std::vector<vertex*> subGV = subGraphs[k]->getVertices();
      std::map<int, edge*>  subGE = subGraphs[k]->getEdges();
      if (subGV.size() >= 3) {
        std::cout << " \n subGraph: " << k << " ::";
        for (unsigned int l = 0; l < subGV.size(); l++) {
          std::cout << " " << subGV[l]->getIndex();
        }
        std::cout << " \n";

        if (!subGE.empty()) {
          for (edgeIterator e = subGE.begin(); e != subGE.end(); e++) {
            edge* pEdge = e->second;
            std::cout << pEdge->v1->getIndex() << "-" << pEdge->v2->getIndex() << "; ";
          }
        }
      }
    }
    std::cout << " \n";
#endif

    // Assign hydrophobicities
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAtom = atomList[i];
      pAtom->setHydrophilicity(atHydrophobicities[i]);
    }
    return 0;
}

}

