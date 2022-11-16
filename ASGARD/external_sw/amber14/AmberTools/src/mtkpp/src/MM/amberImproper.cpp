/*! 
   \file amberImproper.cpp
   \brief AMBER improper energy and gradient
   \author Martin Peters

   $Date: 2007/09/14 10:29:56 $
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

//=============================================================
// description:
//                      ___
// E              (1/2) \  ( V [1+cos(n phi - gamma)])
//  dihedrals  =        /__   n
//                   dihedrals
//
//  V     = barrier height
//   n
//
//  n     = multiplicity or the number of minimum points in the function
//
//  gamma = phase factor
//
//  phi   = torsion angle
//
//=============================================================

#include "amberImproper.h"
#include "amber.h"
#include "math.h"
#include <iostream>

namespace MTKpp
{

// ============================================================
// Function : amberImproper()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberImproper::amberImproper() {}

// ============================================================
// Function : amberImproper()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberImproper::amberImproper(amber *parent):pAmber(parent) {}

// ============================================================
// Function : ~amber()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
amberImproper::~amberImproper(){}

// ============================================================
// Function : calculateE()
// ------------------------------------------------------------ 
//
// ============================================================
double amberImproper::calculateE()
{
#ifdef DEBUG
    std::cout << " amberImproper::calculateE " << std::endl;
#endif

    energy       = 0.0;
    double       dImproper = 0.0;

    int nImpropers = pAmber->getNumImpropers();
    double *xyz = pAmber->getCoords();
    int *impropers = pAmber->getImpropers();
    double *improperParams = pAmber->getImproperParams();

    int j = 0;
    for (int i = 0; i <= nImpropers*4-4; i+=4) {
      if (improperParams[j] > 0) {
        dImproper = improper(xyz[impropers[i]*3  ], xyz[impropers[i]*3+1  ], xyz[impropers[i]*3+2  ],
                             xyz[impropers[i+1]*3], xyz[impropers[i+1]*3+1], xyz[impropers[i+1]*3+2],
                             xyz[impropers[i+2]*3], xyz[impropers[i+2]*3+1], xyz[impropers[i+2]*3+2], 
                             xyz[impropers[i+3]*3], xyz[impropers[i+3]*3+1], xyz[impropers[i+3]*3+2]);

        energy += ( improperParams[j] * (1 + cos(improperParams[j+1] * dImproper - improperParams[j+2]) ));

        //std::cout << " Improper2 " << impropers[i] << "-" << impropers[i+1] << "-" << impropers[i+2] << "-" << impropers[i+3] 
        //          << " " << improperParams[j] << " " << improperParams[j+1] << improperParams[j+2] << " " 
        //          << " " << dImproper << " " << energy << std::endl;
      }
      j+=3;
    }
    return energy;
}

// ============================================================
// Function : calculateG()
// ------------------------------------------------------------ 
//
// ============================================================
double amberImproper::calculateG()
{
#ifdef DEBUG
    std::cout << " amberImproper::calculateG " << std::endl;
#endif
    return 0.0;
}

// ============================================================
// Function : improper()
// ------------------------------------------------------------ 
//
// ============================================================
double amberImproper::improper(double x1, double y1, double z1, double x2, double y2, double z2,
                               double x3, double y3, double z3, double x4, double y4, double z4)
{
    // a-b
    double abX = x1 - x2;
    double abY = y1 - y2;
    double abZ = z1 - z2;

    // c-b
    double cbX = x3 - x2;
    double cbY = y3 - y2;
    double cbZ = z3 - z2;

    // c-d
    double cdX = x3 - x4;
    double cdY = y3 - y4;
    double cdZ = z3 - z4;

    // p = ab cross cb
    double pX = (abY*cbZ - cbY*abZ);
    double pY = (abZ*cbX - cbZ*abX);
    double pZ = (abX*cbY - cbX*abY);

    // normalize p
    double pLength = sqrt(pX*pX + pY*pY + pZ*pZ);
    pX/=pLength;
    pY/=pLength;
    pZ/=pLength;

    // q = cb cross cd
    double qX = (cbY*cdZ - cdY*cbZ);
    double qY = (cbZ*cdX - cdZ*cbX);
    double qZ = (cbX*cdY - cdX*cbY);

    // normalize q
    double qLength = sqrt(qX*qX + qY*qY + qZ*qZ);
    qX/=qLength;
    qY/=qLength;
    qZ/=qLength;

    double ang = acos(pX*qX + pY*qY + pZ*qZ);

    // pXq = p cross q
    double pcqX = (pY*qZ - qY*pZ);
    double pcqY = (pZ*qX - qZ*pX);
    double pcqZ = (pX*qY - qX*pY);

    double s = cbX*pcqX + cbY*pcqY + cbZ*pcqZ;

    if (s < 0) ang = -ang;

    return ang;
}

} // MTKpp namespace

