/*!
   \file amberTorsion.cpp
   \brief AMBER torsion energy and gradient
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
   $Revision: 1.9 $

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

#include "amberTorsion.h"
#include "amber.h"
#include "math.h"
#include "Utils/constants.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : amberTorsion()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberTorsion::amberTorsion() {}

// ============================================================
// Function : amberTorsion()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberTorsion::amberTorsion(amber *parent):pAmber(parent) {}

// ============================================================
// Function : ~amber()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
amberTorsion::~amberTorsion() {}

// ============================================================
// Function : calculateE()
// ------------------------------------------------------------
//
// ============================================================
double amberTorsion::calculateE(int k)
{
    energy = 0.0;

    double *xyz = pAmber->getCoords();
    int nTorsions = 0;
    int *torsions = 0;
    double *torsionParams = 0;

    if (k == 0) {
      nTorsions = pAmber->getNumTorsions();
      torsions = pAmber->getTorsions();
      torsionParams = pAmber->getTorsionParams();
    }
    else {
      nTorsions = pAmber->getNumImpropers();
      torsions = pAmber->getImpropers();
      torsionParams = pAmber->getImproperParams();
    }

/*
   torsionParams[j] = pTorsionParam->Vn / pTorsionParam->npth;
   torsionParams[j+1] = pTorsionParam->Nt;
   torsionParams[j+2] = pTorsionParam->gamma;
*/
    double tor = 0.0;
    int j = 0;
    for (int i = 0; i <= nTorsions*4-4; i+=4) {

      tor = torsion(xyz[torsions[i  ]*3], xyz[torsions[i  ]*3+1], xyz[torsions[i  ]*3+2],
                    xyz[torsions[i+1]*3], xyz[torsions[i+1]*3+1], xyz[torsions[i+1]*3+2],
                    xyz[torsions[i+2]*3], xyz[torsions[i+2]*3+1], xyz[torsions[i+2]*3+2],
                    xyz[torsions[i+3]*3], xyz[torsions[i+3]*3+1], xyz[torsions[i+3]*3+2]);

      energy += ( torsionParams[j] * (1 + cos(torsionParams[j+1] * tor - torsionParams[j+2]) ));

      j+=3;
    }
    return energy;
}

// ============================================================
// Function : calculateG()
// ------------------------------------------------------------
//
// ============================================================
double amberTorsion::calculateG(int k)
{
//#ifdef DEBUG
//    std::cout << " amberTorsion::calculateG " << std::endl;
//#endif
    energy       = 0.0;

    double *xyz = pAmber->getCoords();
    double *grad = pAmber->getGradients();
    int nTorsions = 0;
    int *torsions = 0;
    double *torsionParams = 0;

    if (k == 0) {
      nTorsions = pAmber->getNumTorsions();
      torsions = pAmber->getTorsions();
      torsionParams = pAmber->getTorsionParams();
    }
    else {
      nTorsions = pAmber->getNumImpropers();
      torsions = pAmber->getImpropers();
      torsionParams = pAmber->getImproperParams();
    }


    double x21 = 0.0;
    double y21 = 0.0;
    double z21 = 0.0;

    double x31 = 0.0;
    double y31 = 0.0;
    double z31 = 0.0;

    double x32 = 0.0;
    double y32 = 0.0;
    double z32 = 0.0;

    double x42 = 0.0;
    double y42 = 0.0;
    double z42 = 0.0;

    double x43 = 0.0;
    double y43 = 0.0;
    double z43 = 0.0;

    double xt = 0.0;
    double yt = 0.0;
    double zt = 0.0;

    double xu = 0.0;
    double yu = 0.0;
    double zu = 0.0;

    //double xtu = 0.0;
    //double ytu = 0.0;
    //double ztu = 0.0;

    double rt2 = 0.0;
    double ru2 = 0.0;
    double rtru = 0.0;
    double r32 = 0.0;

    double ang = 0.0;
    double cosAng = 0.0;
    //double sinAng = 0.0;
    double s = 0.0;
    double nPhi = 0.0;
    double cosGamma = 0.0;
    double sinGamma = 0.0;
    double cosNPhi = 0.0;
    double sinNPhi = 0.0;
    double e1 = 0.0;
    double de1 = 0.0;
    double dedphi = 0.0;

    double dPhi_dtx = 0.0;
    double dPhi_dty = 0.0;
    double dPhi_dtz = 0.0;
    double dPhi_dux = 0.0;
    double dPhi_duy = 0.0;
    double dPhi_duz = 0.0;
    double dedx1 = 0.0;
    double dedy1 = 0.0;
    double dedz1 = 0.0;
    double dedx2 = 0.0;
    double dedy2 = 0.0;
    double dedz2 = 0.0;
    double dedx3 = 0.0;
    double dedy3 = 0.0;
    double dedz3 = 0.0;
    double dedx4 = 0.0;
    double dedy4 = 0.0;
    double dedz4 = 0.0;

    int p = -3;
    for (int i = 0; i <= nTorsions*4-4; i+=4) {
      p+=3;
      x21 = xyz[torsions[i+1]*3  ] - xyz[torsions[i]*3  ];
      y21 = xyz[torsions[i+1]*3+1] - xyz[torsions[i]*3+1];
      z21 = xyz[torsions[i+1]*3+2] - xyz[torsions[i]*3+2];

      x32 = xyz[torsions[i+2]*3  ] - xyz[torsions[i+1]*3  ];
      y32 = xyz[torsions[i+2]*3+1] - xyz[torsions[i+1]*3+1];
      z32 = xyz[torsions[i+2]*3+2] - xyz[torsions[i+1]*3+2];

      x43 = xyz[torsions[i+3]*3  ] - xyz[torsions[i+2]*3  ];
      y43 = xyz[torsions[i+3]*3+1] - xyz[torsions[i+2]*3+1];
      z43 = xyz[torsions[i+3]*3+2] - xyz[torsions[i+2]*3+2];

      xt = y21*z32 - y32*z21;
      yt = z21*x32 - z32*x21;
      zt = x21*y32 - x32*y21;

      xu = y32*z43 - y43*z32;
      yu = z32*x43 - z43*x32;
      zu = x32*y43 - x43*y32;

      //xtu = yt*zu - yu*zt;
      //ytu = zt*xu - zu*xt;
      //ztu = xt*yu - xu*yt;

      rt2 = xt*xt + yt*yt + zt*zt;
      ru2 = xu*xu + yu*yu + zu*zu;
      rtru = sqrt(rt2 * ru2);

      if (rtru != 0.0) {
/*
std::cout << xyz[torsions[i]*3  ] << " "
          << xyz[torsions[i]*3+1] << " "
          << xyz[torsions[i]*3+2] << std::endl;
std::cout << xyz[torsions[i+1]*3  ] << " "
          << xyz[torsions[i+1]*3+1] << " "
          << xyz[torsions[i+1]*3+2] << std::endl;
std::cout << xyz[torsions[i+2]*3  ] << " "
          << xyz[torsions[i+2]*3+1] << " "
          << xyz[torsions[i+2]*3+2] << std::endl;
std::cout << xyz[torsions[i+3]*3  ] << " "
          << xyz[torsions[i+3]*3+1] << " "
          << xyz[torsions[i+3]*3+2] << std::endl;
std::cout << torsionParams[p] << " " << torsionParams[p+1] << " " << torsionParams[p+2] <<std::endl;
*/
        r32 = sqrt(x32*x32 + y32*y32 + z32*z32);
        cosAng = -1.0 * (xt*xu + yt*yu + zt*zu)/rtru;
        ang = acos(cosAng);

        if (ang != ang ) {
          std::cout << xyz[torsions[i]*3]  << " "<< xyz[torsions[i]*3+1]<< " " << xyz[torsions[i]*3+2] << std::endl;
          std::cout << xyz[torsions[i+1]*3]<<" "<< xyz[torsions[i+1]*3+1]<< " " << xyz[torsions[i+1]*3+2] << std::endl;
          std::cout << xyz[torsions[i+2]*3]<< " " << xyz[torsions[i+2]*3+1]<< " " << xyz[torsions[i+2]*3+2] << std::endl;
          std::cout << xyz[torsions[i+3]*3]<< " " << xyz[torsions[i+3]*3+1]<< " " << xyz[torsions[i+3]*3+2] << std::endl;

          std::cout << " TORSION/IMPROPER IS NAN ... exiting " << std::endl;
          throw MTKException(" TORSION/IMPROPER IS NAN ... exiting ");
          //exit(0);
        }

        s = x32*(zt*yu - yt*zu) + y32*(xt*zu - zt*xu) + z32*(yt*xu - xt*yu);
        if (s > 0.0) ang = -ang;
        ang = PI - ang;
        cosAng = cos(ang);
        //sinAng = sin(ang);

/*
torsionParams[j]   = pTorsionParam->Vn / pTorsionParam->npth;
torsionParams[j+1] = pTorsionParam->Nt;
torsionParams[j+2] = pTorsionParam->gamma;
*/

        nPhi = torsionParams[p+1] * ang;
        cosGamma = cos(torsionParams[p+2]);
        sinGamma = sin(torsionParams[p+2]);
        cosNPhi = cos(nPhi);
        sinNPhi = sin(nPhi);

        e1 = 1.0 + (cosNPhi*cosGamma + sinNPhi*sinGamma);

        // energy
        energy += (torsionParams[p] * e1);

        // gradient
        de1 = torsionParams[p+1] * (cosNPhi*sinGamma - sinNPhi*cosGamma);

        dedphi = (torsionParams[p] * de1);

        x31 = xyz[torsions[i+2]*3  ] - xyz[torsions[i]*3  ];
        y31 = xyz[torsions[i+2]*3+1] - xyz[torsions[i]*3+1];
        z31 = xyz[torsions[i+2]*3+2] - xyz[torsions[i]*3+2];

        x42 = xyz[torsions[i+3]*3  ] - xyz[torsions[i+1]*3  ];
        y42 = xyz[torsions[i+3]*3+1] - xyz[torsions[i+1]*3+1];
        z42 = xyz[torsions[i+3]*3+2] - xyz[torsions[i+1]*3+2];

        dPhi_dtx =  dedphi * (yt * z32 - y32 * zt) / (rt2 * r32);
        dPhi_dty =  dedphi * (zt * x32 - z32 * xt) / (rt2 * r32);
        dPhi_dtz =  dedphi * (xt * y32 - x32 * yt) / (rt2 * r32);

        dPhi_dux = -dedphi * (yu * z32 - y32 * zu) / (ru2 * r32);
        dPhi_duy = -dedphi * (zu * x32 - z32 * xu) / (ru2 * r32);
        dPhi_duz = -dedphi * (xu * y32 - x32 * yu) / (ru2 * r32);

        dedx1 = dPhi_dty * z32 - dPhi_dtz * y32;
        dedy1 = dPhi_dtz * x32 - dPhi_dtx * z32;
        dedz1 = dPhi_dtx * y32 - dPhi_dty * x32;

        dedx2 = dPhi_dtz * y31 - dPhi_dty * z31 + dPhi_duy * z43 - dPhi_duz * y43;
        dedy2 = dPhi_dtx * z31 - dPhi_dtz * x31 + dPhi_duz * x43 - dPhi_dux * z43;
        dedz2 = dPhi_dty * x31 - dPhi_dtx * y31 + dPhi_dux * y43 - dPhi_duy * x43;

        dedx3 = dPhi_dty * z21 - dPhi_dtz * y21 + dPhi_duz * y42 - dPhi_duy * z42;
        dedy3 = dPhi_dtz * x21 - dPhi_dtx * z21 + dPhi_dux * z42 - dPhi_duz * x42;
        dedz3 = dPhi_dtx * y21 - dPhi_dty * x21 + dPhi_duy * x42 - dPhi_dux * y42;

        dedx4 = dPhi_duy * z32 - dPhi_duz * y32;
        dedy4 = dPhi_duz * x32 - dPhi_dux * z32;
        dedz4 = dPhi_dux * y32 - dPhi_duy * x32;

        grad[torsions[i]*3    ] += dedx1;
        grad[torsions[i]*3+1  ] += dedy1;
        grad[torsions[i]*3+2  ] += dedz1;

        grad[torsions[i+1]*3  ] += dedx2;
        grad[torsions[i+1]*3+1] += dedy2;
        grad[torsions[i+1]*3+2] += dedz2;

        grad[torsions[i+2]*3  ] += dedx3;
        grad[torsions[i+2]*3+1] += dedy3;
        grad[torsions[i+2]*3+2] += dedz3;

        grad[torsions[i+3]*3  ] += dedx4;
        grad[torsions[i+3]*3+1] += dedy4;
        grad[torsions[i+3]*3+2] += dedz4;
      }
    }
    return energy;
}

// ============================================================
// Function : torsion()
// ------------------------------------------------------------
//
// ============================================================
inline double amberTorsion::torsion(double x1, double y1, double z1,
                                    double x2, double y2, double z2,
                                    double x3, double y3, double z3,
                                    double x4, double y4, double z4)
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

    double pXq = pX*qX + pY*qY + pZ*qZ;
    if (std::abs(pXq) > 1) {
      pXq = 1;
    }
    else if (pXq < -1) {
      pXq = -1;
    }

    double ang = acos(pXq);

    // pXq = p cross q
    double pcqX = (pY*qZ - qY*pZ);
    double pcqY = (pZ*qX - qZ*pX);
    double pcqZ = (pX*qY - qX*pY);

    double s = cbX*pcqX + cbY*pcqY + cbZ*pcqZ;

    if (s < 0) ang = -ang;

    return ang;
}

} // MTKpp namespace

