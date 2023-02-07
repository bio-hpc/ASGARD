/*! 
   \file amberBond.h
   \brief AMBER bond energy and gradient
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
   $Revision: 1.8 $

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

#ifndef AMBERBOND_H
#define AMBERBOND_H

#include "Utils/constants.h"

namespace MTKpp
{

class amber;

// ============================================================
// Class : amberBond()
// ------------------------------------------------------------
/*! 
   \class amberBond

   \brief AMBER bond energy and gradient

   \author Martin Peters

   \version 0.1

   \date 2006

   \section AMBERBondEnG AMBER Bond Energy And Gradient

   \f{eqnarray}
     E_{\rm bond} = K_r\left(r_{ij} - r_{eq}\right)^2
   \label{AMBER:bond}
   \f}
   where \f$K_r\f$ is the bond force constant, \f$r_{ij}\f$ is the bond length, \f$r_{eq}\f$ is the standard bond length.

   \f{eqnarray}
   \nabla_iE & = & {\partial E \over \partial r} \cdot \nabla_ir \\
             & = & 2K_r (r_{ij} - r_{eq}) \cdot {r_{ij} \over |r_{ij}|}
   \label{AMBER:dBond}
   \f}

   \f{eqnarray}
   {\partial r \over \partial x_i} & = & {1 \over 2} ((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j))^{-1/2} \cdot 2(x_i - x_j) \\
                & = & {1 \over r} (x_i - x_j) \\
   \nabla_ir       & = & {r_{ij} \over |r_{ij}|} \\
   \nabla_jr       & = & - \nabla_ir
   \label{dr:dx}
   \f}
*/
// ============================================================
class amberBond
{
public:

    /*!
       \brief amberBond Constructor
    */
    amberBond();

    /*!
       \brief amberBond Constructor
       \param pAmber parent potential
    */
    amberBond(amber* pAmber);

    /*!
       \brief amberBond Destructor
    */
    virtual ~amberBond();

    /*!
       \brief Calculate bond Energy
       \return bond energy
    */
    double calculateE();

    /*!
       \brief Calculate bond Energy and Gradient
       \return bond energy
    */
    double calculateG();

protected:
    //! parent potential
    amber*         pAmber;

    // bond energy
    double         energy;

    //! distance
    double         deltaDistance;
};

} // MTKpp namespace

#endif // AMBERBOND_H
