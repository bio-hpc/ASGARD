/*!
   \file amber.cpp
   \brief AMBER MM potential
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

#include "amber.h"
#include "amberBond.h"
#include "amberAngle.h"
#include "amberTorsion.h"
#include "amberImproper.h"
#include "amberNonBonded.h"

namespace MTKpp
{

// ============================================================
// Function : amber()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amber::amber():mmPotential()
{
    pAmberBond = 0;
    pAmberAngle = 0;
    pAmberTorsion = 0;
    pAmberNonBonded = 0;
    this->amberInitialize();
}

// ============================================================
// Function : amberInitialize()
// ------------------------------------------------------------
// 
// ============================================================
void amber::amberInitialize()
{
    pAmberBond      = new amberBond(this);
    pAmberAngle     = new amberAngle(this);
    pAmberTorsion   = new amberTorsion(this);
    pAmberNonBonded = new amberNonBonded(this);
}

// ============================================================
// Function : ~amber()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
amber::~amber()
{
    delete pAmberBond;
    delete pAmberAngle;
    delete pAmberTorsion;
    delete pAmberNonBonded;
}

// ============================================================
// Function : calcBondEnergy()
// ------------------------------------------------------------
// 
// ============================================================
double amber::calcBondEnergy()
{
    if (bGradient) {
      this->bondEnergy = pAmberBond->calculateG();
    }
    else {
      this->bondEnergy = pAmberBond->calculateE();
    }
    return this->bondEnergy;
}

// ============================================================
// Function : calcAngleEnergy()
// ------------------------------------------------------------
// 
// ============================================================
double amber::calcAngleEnergy()
{
    if (bGradient) {
      this->angleEnergy = pAmberAngle->calculateG();
    }
    else {
      this->angleEnergy = pAmberAngle->calculateE();
    }
    return this->angleEnergy;
}

// ============================================================
// Function : calcTorsionEnergy()
// ------------------------------------------------------------
// 
// ============================================================
double amber::calcTorsionEnergy()
{
    this->torsionEnergy = 0.0;
    if (bGradient) {
      this->torsionEnergy = pAmberTorsion->calculateG(0);
    }
    else {
      this->torsionEnergy = pAmberTorsion->calculateE(0);
    }
    return this->torsionEnergy;
}

// ============================================================
// Function : calcImproperEnergy()
// ------------------------------------------------------------
// 
// ============================================================
double amber::calcImproperEnergy()
{
    this->improperEnergy = 0.0;
    if (bGradient) {
      this->improperEnergy = pAmberTorsion->calculateG(1);
    }
    else {
      this->improperEnergy = pAmberTorsion->calculateE(1);
    }
    return this->improperEnergy;
}

// ============================================================
// Function : calcNonBondedEnergy()
// ------------------------------------------------------------
// 
// ============================================================
double amber::calcNonBondedEnergy()
{
    if (bGradient) {
      this->nonBondedEnergy = pAmberNonBonded->calculateG();
    }
    else {
      if (this->bNonBondedDecompose) {
        this->nonBondedEnergy = pAmberNonBonded->decompose();
      }
      else {
        this->nonBondedEnergy = pAmberNonBonded->calculateE();
      }
    }
    return this->nonBondedEnergy;
}

// ============================================================
// Function : calcHBondEnergy()
// ------------------------------------------------------------
// 
// ============================================================
double amber::calcHBondEnergy()
{
    this->hBondEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : printPWD()
// ------------------------------------------------------------
// 
// ============================================================
void amber::printPWD(int i, int j, double r6, double r12, double vdw_energy, double electrostatic, double r)
{
#ifdef HAVE_ZLIB
    gzprintf(pwdGZFileStream,
             "%7d%7d %10.7f %10.7f %10.7f %12.7f %10.4f\n",
             i, j, r6, r12, vdw_energy, electrostatic, r);
#endif
}

// ============================================================
// Functions : nonBonded set functions
// ------------------------------------------------------------
// 
// ============================================================
void amber::setR6(double d){this->R6 = d;}
void amber::setR12(double d){this->R12 = d;}
void amber::setVDW14(double d){this->vdW14Energy = d;}
void amber::setVDW(double d){this->vdWEnergy = d;}
void amber::setEle14(double d){this->ele14Energy = d;}
void amber::setEle(double d){this->eleEnergy = d;}

} // MTKpp namespace
