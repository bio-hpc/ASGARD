/*!
   \file vertex.h
   \brief Container for vertex information
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

#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>

namespace MTKpp
{

// ============================================================
// Class : vertex()
// ------------------------------------------------------------
/*!
   \class vertex
   \brief Container for vertex info
   \author Martin Peters
   \version 0.1
   \date 2006
*/
// ============================================================

class vertex
{
public:

    /*!
       \brief vertex Constructor
    */
    vertex();

    //! vertex Destructor.
    virtual ~vertex();

    /*!
       \brief Set name of vertex
       \param name vertex name
    */
    void                     setName(std::string name);

    /*!
       \brief Get name of vertex
       \return vertex name
    */
    std::string              getName();

    /*!
       \brief Set index of vertex
       \param i internal index
    */
    void                     setIndex(int i);

    /*!
       \brief Get index of vertex
       \return vertex index
    */
    int                      getIndex();

    /*!
       \brief Set Parent vertex
       \param v parent vertex
    */
    void                     setParent(vertex* v);

    /*!
       \brief Get Parent vertex
       \return parent vertex
    */
    vertex*                  getParent();

    /*!
       \brief Set visited boolean
    */
    void                     setVisited();

    /*!
       \brief Set visited boolean
       \param v boolean
    */
    void                     setVisited(bool v);

    /*!
       \brief Get visited boolean
    */
    bool                     isVisited();

    /*!
       \brief Set leaf boolean
    */
    void                     setLeaf();

    /*!
       \brief Get leaf boolean
    */
    bool                     isLeaf();

    /*!
       \brief Set layer value
       \param l layer value
    */
    void                     setLayer(int l);

    /*!
       \brief Get layer value
    */
    int                      getLayer();

    /*!
       \brief Set  value
       \param d value
    */
    void                     setValue(double d);

    /*!
       \brief Get value
    */
    double                   getValue();

    /*!
       \brief Add bonded vertex
       \param v neighboring vertex
    */
    void                     addNeighbor(vertex* v);

    /*!
       \brief Delete bonded vertex
       \param v neighboring vertex
    */
    void                     delNeighbor(vertex* v);

    /*!
       \brief Get bonded vertices
       \return neighboring vertices
    */
    std::vector<vertex*>     getNeighbors();

protected:

    //! vertex name
    std::string              itsName;

    //! parent vertex
    vertex*                  itsParent;

    //! bonded vertices
    std::vector<vertex*>     itsNeighbors;

    //! layer value
    int                      itsLayer;

    //! value
    double                   itsValue;

    //! internal index
    int                      itsIndex;

    //! Was visited or not
    bool                     visited;

    //! leaf or not
    bool                     leaf;

    //! vertex iterator
    typedef std::vector<vertex*>::iterator VertexIterator;
};

} // MTK++ namespace

#endif // VERTEX_H


