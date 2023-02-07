/*! 
   \file gaMutate.cpp
   \brief Performs mutation of gaChromosomes
   \author Martin Peters

   $Date: 2007/09/14 10:28:10 $
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

#include "gaMutate.h"
#include "gaOperators.h"
#include "gaIndividual.h"
#include "gaRanNumGen.h"

namespace MTKpp
{

// ============================================================
// Function : gaMutate()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaMutate::gaMutate(gaOperators *parent, std::vector<int> genePerChr, 
                   std::string mutate, int chrPerInd)
{
    pParent = parent;
    this->mutate = mutate;
    this->chrPerInd = chrPerInd;
    this->genePerChr = genePerChr;
}

// ============================================================
// Function : ~gaMutate()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaMutate::~gaMutate() {}

// ============================================================
// Function : mutate()
// ------------------------------------------------------------
// 
// ============================================================
bool gaMutate::mutation(gaIndividual* newInd)
{
    bool successful = false;
    if (!newInd) return successful;
    if (mutate == "single-gene")    successful = this->__singleGene(newInd);
    return successful;
}

// ============================================================
// Function : __singleGene()
// ------------------------------------------------------------
// Random gene is selected per chromosome
// followed by point mutation
// ============================================================
bool gaMutate::__singleGene(gaIndividual* newInd)
{
    for (int i = 0; i < this->chrPerInd; i++) {
      int mutateGene = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
      newInd->mutateGene(i, mutateGene);
    }
    return true;
}

// ============================================================
// Function : __multipleGene()
// ------------------------------------------------------------
// Multiple Random genes are selected per chromosome 
// followed by point mutation
// ============================================================
bool gaMutate::__multipleGene(gaIndividual* newInd)
{
    for (int i = 0; i < this->chrPerInd; i++) {
      int nGenes = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
      for (int j = 0; j < nGenes; j++) {
        int mutateGene = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
        newInd->mutateGene(i, mutateGene);
      }
    }
    return true;
}

} // MTKpp namespace

