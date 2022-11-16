/*! 
   \file gaGaussian.cpp
   \brief Contains gaussian information used by gaSelection
   \author Martin Peters

   $Date: 2007/09/14 10:28:10 $
   $Revision: 1.4 $

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

#include "gaGaussian.h"
#include "math.h"

namespace MTKpp
{

// ============================================================
// Function : gaGaussian()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaGaussian::gaGaussian()
{
    this->start = 0;
    this->end   = 0;
    this->numIndividuals = 0;
    this->stdDev = 0.0;
    this->invStdDevPi2 = 0.0;
}

// ============================================================
// Function : ~gaGaussian()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaGaussian::~gaGaussian() {}

// ============================================================
// Function : setup()
// ------------------------------------------------------------
// 
// ============================================================
void gaGaussian::setup(const int &st, const int &en, const double &stDe, const int &numInd)
{
    this->start = st;
    this->end = en;
    this->numIndividuals = numInd;
    this->stdDev = stDe * numInd;
    this->invStdDevPi2 = 1 / (this->stdDev * 2.506628);
    this->makeGaussian();
}

// ============================================================
// Function : makeGaussian()
// ------------------------------------------------------------
// 
// ============================================================
void gaGaussian::makeGaussian()
{
    double gaussianSum = 0.0;
    int index = 0;
    for (int i = this->start; i < this->end; i++) {
      this->gaussian.push_back( this->invStdDevPi2 * exp( -(i*i) / (2*this->stdDev*this->stdDev) ) );
      gaussianSum += this->gaussian[index];
      index++;
    }
    index = 0;
    for (int i = this->start; i < this->end; i++) {
      this->gaussian[index] /= gaussianSum;
      index++;
    }
    index = 0;
    gaussianSum = 0.0;
    for (int i = this->start; i < this->end; i++) {
      gaussianSum += this->gaussian[index];
      this->gaussian[index] = 1.0 - gaussianSum;
      index++;
    }
}

// ============================================================
// Function : pick()
// ------------------------------------------------------------
// 
// ============================================================
int gaGaussian::pick(double randomNum)
{
    int index = 0;
    double g = 1.0;
    for (int i = this->start; i < this->end; i++) {
      if ( (randomNum < g) and (randomNum > this->gaussian[index])) {
        return index;
      }
      g = this->gaussian[index];
      index++;
    }
    return index;
}

} // MTKpp namespace

