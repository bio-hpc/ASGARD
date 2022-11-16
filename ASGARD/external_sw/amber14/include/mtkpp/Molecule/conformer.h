/*! 
   \file conformer.h
   \brief Rotatable bond based conformer class
   \author Martin Peters

    - For each torsion we will only rotate on the rhs, X1, 4, X2, ... etc.
    \code
                        X1
                      /
                2 -- 3     X2
              /       \  /
            1          4
    \endcode
    - We need to make a list of these atoms for each torsion

    All torsions are stored in radians:
     360 deg = 2 pi radian
     180 deg = pi radian
      90 deg = pi/2 radian
       1 deg = pi/180 radian = 0.0174533 radian
      10 deg = 10 * 0.0174533 = 0.174533
    1 radian = 180/pi deg = 57.29578 deg

   $Date: 2010/03/29 20:42:28 $
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

#ifndef CONFORMER_H
#define CONFORMER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

namespace MTKpp
{

class molecule;
class atom;
class Bond;
class vector3d;
struct Torsion;

// ============================================================
// Struct : conformer
// ------------------------------------------------------------
/*! 
   \struct conformer
   \brief Container for conformer info
   \author Martin Peters
   \date 2005
*/
// ============================================================
struct conformer
{
    //! conformer name
    std::string name;

    //! torsion values
    std::vector<double> torsions;
};

// ============================================================
// Class : conformers()
// ------------------------------------------------------------
/*! 
   \class conformers
   \brief This class creates and operators on conformers, however, the conformers are stored in molecule
   \author Martin Peters
*/
// ============================================================
class conformers
{
public:
    /*!
       \brief conformers Constructor
       \param parent molecule pointer
    */
    conformers(molecule *parent = 0);

    //! conformers Destructor.
    virtual ~conformers();

    /*!
       \brief Set rotatable bonds
       \param a atom1 index
       \param b atom2 index
       \param c atom3 index
       \param d atom4 index
    */
    void addRotatableTorsion(int a, int b, int c, int d);

    /*!
       \brief Determine atom list
    */
    void determineAtomList();

    /*!
       \brief Get the coordinates of the conformer
       \param c conformer to be considered
       \param newCoords conformer coordinates to be returned
    */
    void getConformerCoordinates(conformer* c, std::vector< vector3d > &newCoords);

    /*!
       \brief Get all atomic coordinates
       \param c conformer pointer
       \param coords double array of coordinates
    */
    void getConformerCoordinates(conformer* c, double coords[][3]);

    /*!
       \brief Get the center of mass
       \param coords conformer coordinates to be considered
       \param center conformer center of mass
    */
    void centerOfMass(std::vector< vector3d > coords, double center[3]);

protected:

    //! molecule pointer
    molecule*                     pParent;

    //! atom iterator
    std::vector<atom*>::iterator  atomIterator;

    // Bond iterator
    std::vector<Bond*>::iterator  bondIterator;

    // int iterator
    std::vector<int>::iterator    intIterator;

    //! atom pointer
    atom*                         pAtom1;

    //! atom pointer
    atom*                         pAtom2;

    //! Bond pointer
    Bond*                         pBond;

    //! parent molecule coordinates
    std::vector< vector3d >       itsCoords;

    //! rotatable torsions
    std::vector< std::vector<int> >    itsRotatableTorsions;

    /*!
       \brief vector of atom vectors

       - For each torsion we will only rotate on the rhs, 4, X, etc
       \code
                2 -- 3     X
              /       \  /
            1          4
       \endcode
       - We need to make a list of these atoms for each torsion
    */
    std::vector< std::vector<unsigned int> > rhsAtomList;
};

} // MTKpp namespace

#endif // CONFORMER_H
