/*!
   \file graph.cpp
   \brief Container for graph
   \author Martin Peters

   $Date: 2010/03/29 20:25:55 $
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

#include "Utils/constants.h"

#include "graph.h"
#include "vertex.h"
#include "edge.h"

#include "Utils/index.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : graph()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
graph::graph()
{
    vertexIndex = 1;
    maxEdges = MAXATOMS;
    counter = 0;
}

// ============================================================
// Function : graph()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
graph::graph(graph* g)
{
    vertexIndex = g->vertexIndex;
    counter = g->counter;
    maxEdges = MAXATOMS;

    std::vector<vertex*> vs = g->getVertices();
    for (unsigned int i = 0; i < vs.size(); i++) {
      this->addVertex(vs[i]);
    }

    for (unsigned int i = 0; i < vs.size(); i++) {
      int i_index = vs[i]->getIndex();
      for (unsigned int j = i+1; j < vs.size(); j++) {
        int j_index = vs[j]->getIndex();
        if (g->hasEdge(i_index, j_index)) {
          this->addEdge(this->getVertex(i_index), this->getVertex(j_index));
        }
      }
    }

}

// ============================================================
// Function : ~graph()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
graph::~graph()
{
    itsVertices.erase(itsVertices.begin(), itsVertices.end());
    itsVertices.clear();

    itsEdges.erase(itsEdges.begin(), itsEdges.end());
    itsEdges.clear();
}

// ============================================================
// Function : addVertex()
// ------------------------------------------------------------
// Add vertex
// ============================================================
vertex* graph::addVertex()
{
    pVertex = new vertex();
    pVertex->setIndex(this->vertexIndex);
    this->vertexIndex++;
    this->itsVertices.push_back(pVertex);
    return pVertex;
}

// ============================================================
// Function : addVertex()
// ------------------------------------------------------------
// Add vertex
// ============================================================
vertex* graph::addVertex(vertex* v)
{
    pVertex = new vertex();
    pVertex->setIndex(v->getIndex());
    if (v->getIndex() > this->vertexIndex) {
      this->vertexIndex = v->getIndex()+1;
    }
    this->itsVertices.push_back(pVertex);
    return pVertex;
}

// ============================================================
// Function : addVertex(int)
// ------------------------------------------------------------
// Add vertex
// ============================================================
vertex* graph::addVertex(int index)
{
    pVertex = new vertex();
    pVertex->setIndex(index);
    this->itsVertices.push_back(pVertex);
    return pVertex;
}

// ============================================================
// Function : getVertex()
// ------------------------------------------------------------
// Get vertex
// ============================================================
vertex* graph::getVertex(std::string name)
{
    for (unsigned int i = 0; i < this->itsVertices.size(); i++) {
      if (this->itsVertices[i]->getName() == name) {
        return this->itsVertices[i];
      }
    }
    return 0;
}

// ============================================================
// Function : getVertex()
// ------------------------------------------------------------
// Get vertex
// ============================================================
vertex* graph::getVertex(int index)
{
    for (unsigned int i = 0; i < this->itsVertices.size(); i++) {
      if (this->itsVertices[i]->getIndex() == index) {
        return this->itsVertices[i];
      }
    }
    return 0;
}

// ============================================================
// Function : getVertices()
// ------------------------------------------------------------
// Get vertices
// ============================================================
std::vector<vertex*> graph::getVertices()
{
    return this->itsVertices;
}

// ============================================================
// Function : addEdge()
// ------------------------------------------------------------
// Add edge
// ============================================================
void graph::addEdge(vertex* vt1, vertex* vt2)
{
    pEdge= new edge();
    pEdge->v1 = vt1;
    pEdge->v2 = vt2;
    pEdge->visited = false;
    vt1->addNeighbor(vt2);
    vt2->addNeighbor(vt1);

    int edgeIndex = indexAB(vt1->getIndex(), vt2->getIndex(), maxEdges);
    itsEdges[edgeIndex] = pEdge;
}

// ============================================================
// Function : getEdge()
// ------------------------------------------------------------
// Get edge
// ============================================================
edge* graph::getEdge(vertex* vt1, vertex* vt2)
{
    if (vt1 == vt2) return 0;

    int edgeIndex = indexAB(vt1->getIndex(), vt2->getIndex(), maxEdges);
    edgeIterator e = this->itsEdges.find(edgeIndex);

    if (e != this->itsEdges.end()){
      return this->itsEdges[edgeIndex];
    }
    return 0;
}

// ============================================================
// Function : hasEdge()
// ------------------------------------------------------------
// has edge
// ============================================================
bool graph::hasEdge(int i, int j)
{
    if (i == j) return false;

    int edgeIndex = indexAB(i, j, maxEdges);
    edgeIterator e = this->itsEdges.find(edgeIndex);

    if (e != this->itsEdges.end()){
      return true;
    }
    return false;
}

// ============================================================
// Function : delEdge()
// ------------------------------------------------------------
// Delete edge
// ============================================================
void graph::delEdge(int i, int j)
{
    if (i != j) {
      int edgeIndex = indexAB(i, j, maxEdges);
      edgeIterator e = this->itsEdges.find(edgeIndex);

      if (e != this->itsEdges.end()) {
        pEdge = e->second;
        pEdge->v1->delNeighbor(pEdge->v2);
        pEdge->v2->delNeighbor(pEdge->v1);
        this->itsEdges.erase(e);
      }
    }
}

// ============================================================
// Function : getEdges()
// ------------------------------------------------------------
// Get edges
// ============================================================
std::map<int, edge*> graph::getEdges()
{
    return this->itsEdges;
}

// ============================================================
// Function : dfs()
// ------------------------------------------------------------
// Depth-first search
// ============================================================
void graph::dfs(vertex* v)
{
    int vLayer = v->getLayer();
    v->setVisited();

    if (v->isLeaf()) {
      counter++;
    }

    std::vector<vertex*> n = v->getNeighbors();
    for (unsigned int i = 0; i < n.size(); i++) {
      if (n[i]->getLayer() >= vLayer and !n[i]->isVisited()) {
        pEdge = this->getEdge(v, n[i]);
        if (pEdge) {
          if (!pEdge->visited) { // new
            pEdge->visited = true;
            this->dfs(n[i]);
          }
        }
        //else {
        //  std::cout << " Can not find edge in graph::dfs ... stopping " << std::endl;
        //  throw MTKException(" Can not find edge in graph::dfs ... stopping ");
        //}
      }
    }
}

// ============================================================
// Function : reset()
// ------------------------------------------------------------
// Reset graph
// ============================================================
void graph::reset()
{
    for (unsigned int i = 0; i < this->itsVertices.size(); i++) {
      this->itsVertices[i]->setVisited(false);
    }
    if (!itsEdges.empty()) {
      for (edgeIterator e = this->itsEdges.begin(); e != this->itsEdges.end(); e++) {
        pEdge = e->second;
        pEdge->visited = false;
      }
    }
}

// ============================================================
// Function : getNumVertices()
// ------------------------------------------------------------
// Return number of vertices
// ============================================================
unsigned int graph::getNumVertices()
{
    return this->itsVertices.size();
}

} // MTK++ namespace
