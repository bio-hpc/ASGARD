/*!
   \file vertex.cpp
   \brief Container for vertex information
   \author Martin Peters

   $Date: 2010/08/11 21:15:03 $
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

#include "vertex.h"

namespace MTKpp
{

// ============================================================
// Function : vertex()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
vertex::vertex()
{
    this->itsName = "";
    this->itsParent = 0;
    this->itsLayer = 0;
    this->itsIndex = 0;
    this->itsValue = 0;
    this->visited = false;
    this->leaf = false;
}

// ============================================================
// Function : ~vertex()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
vertex::~vertex() {}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// Set name
// ============================================================
void vertex::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// Get name
// ============================================================
std::string vertex::getName()
{
    return this->itsName;
}

// ============================================================
// Function : setIndex()
// ------------------------------------------------------------
// Set index
// ============================================================
void vertex::setIndex(int i)
{
    this->itsIndex = i;
}

// ============================================================
// Function : getIndex()
// ------------------------------------------------------------
// Get index
// ============================================================
int vertex::getIndex()
{
    return this->itsIndex;
}

// ============================================================
// Function : setParent()
// ------------------------------------------------------------
// Set parent vertex
// ============================================================
void vertex::setParent(vertex* v)
{
    this->itsParent = v;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
// Get parent vertex
// ============================================================
vertex* vertex::getParent()
{
    return this->itsParent;
}

// ============================================================
// Function : setVisited()
// ------------------------------------------------------------
// Sets visited
// ============================================================
void vertex::setVisited()
{
    this->visited = true;
}

// ============================================================
// Function : setVisited(bool)
// ------------------------------------------------------------
// Sets visited
// ============================================================
void vertex::setVisited(bool v)
{
    this->visited = v;
}

// ============================================================
// Function : isVisited()
// ------------------------------------------------------------
// Gets visited
// ============================================================
bool vertex::isVisited()
{
    return this->visited;
}

// ============================================================
// Function : setLeaf()
// ------------------------------------------------------------
// Sets leaf
// ============================================================
void vertex::setLeaf()
{
    this->leaf = true;
}

// ============================================================
// Function : isLeaf()
// ------------------------------------------------------------
// Gets leaf
// ============================================================
bool vertex::isLeaf()
{
    return this->leaf;
}

// ============================================================
// Function : setLayer()
// ------------------------------------------------------------
// Sets layer value
// ============================================================
void vertex::setLayer(int l)
{
    this->itsLayer = l;
}

// ============================================================
// Function : getLayer()
// ------------------------------------------------------------
// Gets layer value
// ============================================================
int vertex::getLayer()
{
    return this->itsLayer;
}

// ============================================================
// Function : setValue()
// ------------------------------------------------------------
// Sets value
// ============================================================
void vertex::setValue(double d)
{
    this->itsValue = d;
}

// ============================================================
// Function : getValue()
// ------------------------------------------------------------
// Gets value
// ============================================================
double vertex::getValue()
{
    return this->itsValue;
}

// ============================================================
// Function : addNeighbor()
// ------------------------------------------------------------
// Add neighbor
// ============================================================
void vertex::addNeighbor(vertex* v)
{
    this->itsNeighbors.push_back(v);
}

// ============================================================
// Function : delNeighbor()
// ------------------------------------------------------------
// Delete neighbor
// ============================================================
void vertex::delNeighbor(vertex* v)
{
    std::vector<vertex*>::iterator nb = itsNeighbors.begin();
    std::vector<vertex*>::iterator ne = itsNeighbors.end();

    VertexIterator c = std::find(nb, ne, v);

    if (c != ne) {
      itsNeighbors.erase(c);
    }
}

// ============================================================
// Function : getNeighbors()
// ------------------------------------------------------------
// Get bonded vertices
// ============================================================
std::vector<vertex*> vertex::getNeighbors()
{
    return this->itsNeighbors;
}

} // MTK++ namespace
