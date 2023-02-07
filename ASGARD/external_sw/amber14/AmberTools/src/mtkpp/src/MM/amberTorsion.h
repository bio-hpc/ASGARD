/*! 
   \file amberTorsion.h
   \brief AMBER torsion energy and gradient
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
   $Revision: 1.7 $

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

#ifndef AMBERTORSION_H
#define AMBERTORSION_H

namespace MTKpp
{

class amber;

// ============================================================
// Class : amberTorsion()
// ------------------------------------------------------------
/*! 
   \class amberTorsion

   \brief AMBER torsion energy and gradient

   \author Martin Peters

   \version 0.1

   \date 2006

   \section AMBERTorsionEnG AMBER Torsion Energy And Gradient

*/
// ============================================================
class amberTorsion
{
public:
    /*!
       \brief amberTorsion Constructor
    */
    amberTorsion();

    /*!
       \brief amberTorsion Constructor
       \param pAmber parent potential
    */
    amberTorsion(amber* pAmber);

    /*!
       \brief amberTorsion Destructor
    */
    virtual ~amberTorsion();

    /*!
       \brief Calculate torsion Energy
       \return torsion energy
    */
    double calculateE(int i);

    /*!
       \brief Calculate torsion Energy and Gradient
       \return torsion energy
    */
    double calculateG(int i);

protected: // FUNCTIONS
    /*!
       \brief Determine the torsion between atoms 1, 2, 3, and 4
       \param x1 atom 1 x coordinate
       \param y1 atom 1 y coordinate
       \param z1 atom 1 z coordinate
       \param x2 atom 2 x coordinate
       \param y2 atom 2 y coordinate
       \param z2 atom 2 z coordinate
       \param x3 atom 3 x coordinate
       \param y3 atom 3 y coordinate
       \param z3 atom 3 z coordinate
       \param x4 atom 4 x coordinate
       \param y4 atom 4 y coordinate
       \param z4 atom 4 z coordinate
       \return torsion between atoms 1, 2, 3, and 4
    */
    inline double torsion(double x1, double y1, double z1,
                          double x2, double y2, double z2,
                          double x3, double y3, double z3,
                          double x4, double y4, double z4);

protected:
    //! torsion energy
    double         energy;

    //! parent potential
    amber*         pAmber;
};

} // MTKpp namespace

#endif // AMBERTORSION_H

