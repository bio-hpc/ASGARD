/*! 
   \file amberBond.cpp
   \brief AMBER bond energy and gradient
   \author Martin Peters

   $Date: 2007/09/14 10:29:56 $
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

//=============================================================
// description:
//
//           ___          2
// E      =  \  ( K (r-r )  )
//  bonds    /__   r    0 
//          bonds
//
//=============================================================

#include "amberBond.h"
#include "amber.h"
#include "math.h"
#include <iostream>

namespace MTKpp
{

// ============================================================
// Function : amberBond()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberBond::amberBond() {}

// ============================================================
// Function : amberBond()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberBond::amberBond(amber *parent):pAmber(parent) {}

// ============================================================
// Function : ~amber()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
amberBond::~amberBond() {}

// ============================================================
// Function : calculateE()
// ------------------------------------------------------------
//
// ============================================================
double amberBond::calculateE()
{
    energy        = 0.0;
    deltaDistance = 0.0;

    int nBonds = pAmber->getNumBonds();
    double *xyz = pAmber->getCoords();

    int *bonds = pAmber->getBonds();
    double *bondParams = pAmber->getBondParams();

    for (int i = 0; i <= nBonds*2-2; i+=2) {
      deltaDistance = (sqrt (pow((xyz[bonds[i]*3  ] - xyz[bonds[i+1]*3  ]),2) +
                             pow((xyz[bonds[i]*3+1] - xyz[bonds[i+1]*3+1]),2) +
                             pow((xyz[bonds[i]*3+2] - xyz[bonds[i+1]*3+2]),2) 
                                   ) - bondParams[i]);
      energy += bondParams[i+1] * deltaDistance * deltaDistance;
    }
    return energy;
}

// ============================================================
// Function : calculateG()
// ------------------------------------------------------------
//
// ============================================================
double amberBond::calculateG()
{
    energy        = 0.0;
    deltaDistance = 0.0;

    int nBonds = pAmber->getNumBonds();
    double *xyz = pAmber->getCoords();
    double *grad = pAmber->getGradients();

    int *bonds = pAmber->getBonds();
    double *bondParams = pAmber->getBondParams();

    double x12 = 0.0;
    double y12 = 0.0;
    double z12 = 0.0;
    double dist12 = 0.0;
    double dEdrr = 0.0;
    double deltaDist = 0.0;
    double dX = 0.0;
    double dY = 0.0;
    double dZ = 0.0;

    for (int i = 0; i <= nBonds*2-2; i+=2) {
      x12 = xyz[bonds[i]*3  ] - xyz[bonds[i+1]*3  ];
      y12 = xyz[bonds[i]*3+1] - xyz[bonds[i+1]*3+1];
      z12 = xyz[bonds[i]*3+2] - xyz[bonds[i+1]*3+2];

      dist12 = sqrt(x12*x12 + y12*y12 + z12*z12);
      deltaDist = dist12 - bondParams[i];

      // Energy
      energy += bondParams[i+1] * deltaDist * deltaDist;

      // Gradient
      dEdrr = 2 * bondParams[i+1] * deltaDist/dist12;

      dX = dEdrr * x12;
      dY = dEdrr * y12;
      dZ = dEdrr * z12;

      grad[bonds[i]*3    ] += dX;
      grad[bonds[i]*3+1  ] += dY;
      grad[bonds[i]*3+2  ] += dZ;
      grad[bonds[i+1]*3  ] -= dX;
      grad[bonds[i+1]*3+1] -= dY;
      grad[bonds[i+1]*3+2] -= dZ;

    }
    return energy;
}

} // MTKpp namespace
