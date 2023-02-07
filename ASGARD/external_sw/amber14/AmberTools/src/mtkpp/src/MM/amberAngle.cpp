/*! 
   \file amberAngle.cpp
   \brief AMBER angle energy and gradient
   \author Martin Peters

   $Date: 2010/04/29 19:46:03 $
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

//=============================================================
// description:
//             ___                   2
// E        =  \  (K (theta - theta )  )
//  angles     /__  theta          0
//            angles
//
//=============================================================

#include "amberAngle.h"
#include "amber.h"
#include <math.h>
#include <cmath>
#include <iostream>

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : amberAngle()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberAngle::amberAngle() {}

// ============================================================
// Function : amberAngle()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberAngle::amberAngle(amber *parent):pAmber(parent) {}

// ============================================================
// Function : ~amber()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
amberAngle::~amberAngle(){}

// ============================================================
// Function : calculateE()
// ------------------------------------------------------------ 
//
// ============================================================
double amberAngle::calculateE()
{
//#ifdef DEBUG
//        std::cout << " amberAngle::calculateE " << std::endl;
//#endif

    energy     = 0.0;
    deltaAngle = 0.0;
    bError = false;

    int nAngles = pAmber->getNumAngles();
    double *xyz = pAmber->getCoords();
    int *angles = pAmber->getAngles();
    double *angleParams = pAmber->getAngleParams();

    int j = 0;
    for (int i = 0; i <= nAngles*3-3; i+=3) {
      if (angleParams[j] > 0) {
        //std::cout << "amberAngle::calculateE angle = " << angles[i] << "-" << angles[i+1] << "-" << angles[i+2] << std::endl;
        deltaAngle = angle(
         xyz[angles[i]*3  ], xyz[angles[i]*3+1  ], xyz[angles[i]*3+2  ],
         xyz[angles[i+1]*3], xyz[angles[i+1]*3+1], xyz[angles[i+1]*3+2],
         xyz[angles[i+2]*3], xyz[angles[i+2]*3+1], xyz[angles[i+2]*3+2], bError) -
         angleParams[j];

        energy += angleParams[j+1] * deltaAngle * deltaAngle;
      }
      j+=2;
    }
    return energy;
}

// ============================================================
// Function : calculateG()
// ------------------------------------------------------------ 
//
// ============================================================
double amberAngle::calculateG()
{
//#ifdef DEBUG
//        std::cout << " amberAngle::calculateG " << std::endl;
//#endif
    energy     = 0.0;
    deltaAngle = 0.0;

    int nAngles = pAmber->getNumAngles();
    double *xyz = pAmber->getCoords();
    double *grad = pAmber->getGradients();
    int *angles = pAmber->getAngles();
    double *angleParams = pAmber->getAngleParams();

    double x12 = 0.0;
    double y12 = 0.0;
    double z12 = 0.0;

    double x32 = 0.0;
    double y32 = 0.0;
    double z32 = 0.0;

    double r12 = 0.0;
    double r32 = 0.0;
    double r13 = 0.0;
    double cosAng = 0.0;
    double ang = 0.0;
    double dt = 0.0;
    double dt1 = 0.0;
    double d13 = 0.0;
    double d11 = 0.0;
    double d33 = 0.0;

    double dx1 = 0.0;
    double dy1 = 0.0;
    double dz1 = 0.0;
    double dx2 = 0.0;
    double dy2 = 0.0;
    double dz2 = 0.0;
    double dx3 = 0.0;
    double dy3 = 0.0;
    double dz3 = 0.0;

    int j = 0;

    for (int i = 0; i <= nAngles*3-3; i+=3) {
      if (angleParams[j] > 0) {
/*
std::cout << xyz[angles[i]*3  ] << " " 
          << xyz[angles[i]*3+1] << " " 
          << xyz[angles[i]*3+2] << std::endl;
std::cout << xyz[angles[i+1]*3  ] << " " 
          << xyz[angles[i+1]*3+1] << " " 
          << xyz[angles[i+1]*3+2] << std::endl;
std::cout << xyz[angles[i+2]*3  ] << " " 
          << xyz[angles[i+2]*3+1] << " " 
          << xyz[angles[i+2]*3+2] << std::endl;
std::cout << angleParams[j] << " " << angleParams[j+1] << std::endl;
*/
        // vectors 1-2 and 3-2
        x12 = xyz[angles[i]*3  ] - xyz[angles[i+1]*3  ];
        y12 = xyz[angles[i]*3+1] - xyz[angles[i+1]*3+1];
        z12 = xyz[angles[i]*3+2] - xyz[angles[i+1]*3+2];

        x32 = xyz[angles[i+2]*3  ] - xyz[angles[i+1]*3  ];
        y32 = xyz[angles[i+2]*3+1] - xyz[angles[i+1]*3+1];
        z32 = xyz[angles[i+2]*3+2] - xyz[angles[i+1]*3+2];

        r12 = x12*x12 + y12*y12 + z12*z12;
        r32 = x32*x32 + y32*y32 + z32*z32;

        r13 = sqrt(r12*r32);
        cosAng = (x12*x32 + y12*y32 + z12*z32)/r13;
        ang = acos(cosAng);

        if (ang != ang) {
          std::cout << xyz[angles[i]*3]  << " "<< xyz[angles[i]*3+1]<< " " << xyz[angles[i]*3+2] << std::endl;
          std::cout << xyz[angles[i+1]*3]<<" "<< xyz[angles[i+1]*3+1]<< " " << xyz[angles[i+1]*3+2] << std::endl;
          std::cout << xyz[angles[i+2]*3]<< " " << xyz[angles[i+2]*3+1]<< " " << xyz[angles[i+2]*3+2] << std::endl;
/*
std::cout << "amberAngle::calculateG (DEBUGGING ANGLE IS NAN ... exiting) " << std::endl;
exit(0);
*/
          throw MTKException(" ANGLE IS NAN ... exiting ");
        }

        //deltaAngle = ang - angleParams[j];
        deltaAngle = angleParams[j] - ang;

        // energy
        energy += angleParams[j+1] * deltaAngle * deltaAngle;

        // gradient
        dt = (-2 * angleParams[j+1] * deltaAngle)/sin(ang);

        //dt1 = dt/cosAng;
        dt1 = dt*cosAng;
        d13 = dt/r13;
        d11 = dt1/r12;
        d33 = dt1/r32;

        dx1 = d13*x32 - d11*x12;
        dy1 = d13*y32 - d11*y12;
        dz1 = d13*z32 - d11*z12;

        dx3 = d13*x12 - d33*x32;
        dy3 = d13*y12 - d33*y32;
        dz3 = d13*z12 - d33*z32;

        dx2 = -dx1-dx3;
        dy2 = -dy1-dy3;
        dz2 = -dz1-dz3;

        grad[angles[i]*3    ] -= dx1;
        grad[angles[i]*3+1  ] -= dy1;
        grad[angles[i]*3+2  ] -= dz1;

        grad[angles[i+1]*3  ] -= dx2;
        grad[angles[i+1]*3+1] -= dy2;
        grad[angles[i+1]*3+2] -= dz2;

        grad[angles[i+2]*3  ] -= dx3;
        grad[angles[i+2]*3+1] -= dy3;
        grad[angles[i+2]*3+2] -= dz3;
      }
      j+=2;
    }

    return energy;
}

// ============================================================
// Function : angle()
// ------------------------------------------------------------
// see vector3d::angle for details
// ============================================================
double amberAngle::angle(double x1, double y1, double z1,
                         double x2, double y2, double z2,
                         double x3, double y3, double z3,
                         bool &bE)
{
    double abX = x1 - x2;
    double abY = y1 - y2;
    double abZ = z1 - z2;
    double abLength = sqrt(abX*abX + abY*abY + abZ*abZ);

    double cbX = x3 - x2;
    double cbY = y3 - y2;
    double cbZ = z3 - z2;
    double cbLength = sqrt(cbX*cbX + cbY*cbY + cbZ*cbZ);

    abX = abX/abLength;
    abY = abY/abLength;
    abZ = abZ/abLength;

    cbX = cbX/cbLength;
    cbY = cbY/cbLength;
    cbZ = cbZ/cbLength;

    double ang = acos(abX*cbX + abY*cbY + abZ*cbZ);

    if (ang != ang) {
/*
      std::cout << "amberAngle::angle (DEBUGGING "<< x1 << ", " << y1 << ", " << z1 << ")" << std::endl;
      std::cout << "amberAngle::angle (DEBUGGING "<< x2 << ", " << y2 << ", " << z2 << ")" << std::endl;
      std::cout << "amberAngle::angle (DEBUGGING "<< x3 << ", " << y3 << ", " << z3 << ")" << std::endl;

      std::cout << "amberAngle::angle (DEBUGGING ANGLE IS NAN) ... exiting ";
      exit(1);
*/
      bE = true;
    }
    return ang;
}

} // MTKpp namespace

