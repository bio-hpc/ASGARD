/*! 
   \file gaSelection.cpp
   \brief Selects gaIndividuals from gaPopulation based on fitness
   \author Martin Peters

   $Date: 2007/09/14 10:28:10 $
   $Revision: 1.3 $

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


#include "gaSelection.h"
#include "gaOperators.h"
#include "gaGaussian.h"
#include "gaRanNumGen.h"

namespace MTKpp
{

// ============================================================
// Function : gaSelection()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaSelection::gaSelection(gaOperators *parent, std::string selection, int maxInds) {

    this->pParent = parent;
    this->selection = selection;
    this->maxInds = maxInds;

    pGaGaussian = new gaGaussian();
}

// ============================================================
// Function : ~gaSelection()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaSelection::~gaSelection() {}

// ============================================================
// Function : select()
// ------------------------------------------------------------
// Main function of gaSelection
// ============================================================
int gaSelection::select()
{
    if (selection == "random") return this->__random();
    if (selection == "semi-random") return this->__semi_random();
    if (selection == "rouletteWheel") return this->__rouletteWheel();
    if (selection == "tournament") return this->__tournament();
    if (selection == "truncation") return this->__truncation();
    return 0;
}

// ============================================================
// Function : setGaussian()
// ------------------------------------------------------------
// If semi-random is used
// ============================================================
void gaSelection::setGaussian(int start, int end, double selPress)
{
    pGaGaussian->setup(start, end, selPress, this->maxInds);
}

// ============================================================
// Function : __random()
// ------------------------------------------------------------
// Returns uniform random individuals
// ============================================================
int gaSelection::__random()
{
    return int(ranNumBetweenZeroAndX(this->maxInds));
}

// ============================================================
// Function : __semi_random()
// ------------------------------------------------------------
// returns semi-random individuals weighted by a normal gaussian distribution
// ============================================================
int gaSelection::__semi_random()
{
    double r = ranNumBetweenZeroAndOne();
    return pGaGaussian->pick(r);
}

// ============================================================
// Function : __rouletteWheel()
// ------------------------------------------------------------
// returns individuals using proportional fitness/roulette wheel
// ============================================================
int gaSelection::__rouletteWheel()
{
    return 1;
}

// ============================================================
// Function : __tournament()
// ------------------------------------------------------------
// returns individuals by perform a tournament (size N) of fitness
// ============================================================
int gaSelection::__tournament()
{
    return 1;
}

// ============================================================
// Function : __truncation()
// ------------------------------------------------------------
// Used by the keep operator
// a certain percentage of the generation is allowed to live
// without competition or selection
// ============================================================
int gaSelection::__truncation()
{
    return 1;
}

} // MTKpp namespace
