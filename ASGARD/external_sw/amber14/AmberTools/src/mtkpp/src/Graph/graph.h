/*!
   \file graph.h
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

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

namespace MTKpp
{

class vertex;
struct edge;

// ============================================================
// Class : graph()
// ------------------------------------------------------------
/*!
   \class graph
   \brief Container for graphs
   \author Martin Peters
   \version 0.1
   \date 2006
*/
// ============================================================

class graph
{
public:

    /*!
       \brief graph Constructor
    */
    graph();

    /*!
       \brief Copy Constructor
       \param g graph pointer
    */
    graph(graph* g);

    //! graph Destructor
    virtual ~graph();

    /*!
       \brief add vertex
       \return vertex pointer
    */
    vertex*                  addVertex();

    /*!
       \brief add vertex
       \param index vertex index
       \return vertex pointer
    */
    vertex*                  addVertex(int index);

    /*!
       \brief add vertex
       \param v vertex pointer
       \return vertex pointer
    */
    vertex*                  addVertex(vertex* v);

    /*!
       \brief get vertex
       \param n name
       \return vertex pointer
    */
    vertex*                  getVertex(std::string n);

    /*!
       \brief get vertex
       \param n index
       \return vertex pointer
    */
    vertex*                  getVertex(int n);

    /*!
       \brief get vertices
       \return vertex list
    */
    std::vector<vertex*>     getVertices();

    /*!
       \brief Add an edge between a and b
       \param a vertex pointer 1
       \param b vertex pointer 2
    */
    void                     addEdge(vertex* a, vertex* b);

    /*!
       \brief Get an edge between a and b
       \param a vertex pointer 1
       \param b vertex pointer 2
       \return edge pointer
    */
    edge*                    getEdge(vertex* a, vertex* b);

    /*!
       \brief Get an edge between a and b
       \param i vertex index 1
       \param j vertex index 2
       \return bool
    */
    bool                     hasEdge(int i, int j);

    /*!
       \brief Delete an edge between a and b
       \param i vertex index 1
       \param j vertex index 2
       \return bool
    */
    void                     delEdge(int i, int j);

    /*!
       \brief get edeges
       \return edge list
    */
    std::map<int, edge*>     getEdges();

    /*!
       \brief Depth-first search
       \param v node vertex pointer
    */
    void                     dfs(vertex* v);

    /*!
       \brief Rest
    */
    void                     reset();

    //!
    int                      counter;

    //! 
    unsigned int             getNumVertices();

protected:

    //! vertex pointer
    vertex* pVertex;
 
    //! edge pointer
    edge* pEdge;

    //! list of vertices
    std::vector<vertex*>     itsVertices;

    //! edge map
    std::map<int, edge*>     itsEdges;

    //! vertex index
    int                      vertexIndex;

    //! edge index
    int                      edgeIndex;

    //! max edges
    unsigned int             maxEdges;

    //! vertices iterator
    //typedef std::vector<vertex*>::iterator vertexIterator;

    //! edge map iterator
    typedef std::map<int, edge*>::iterator edgeIterator;
};

} // MTK++ namespace

#endif // GRAPH_H

