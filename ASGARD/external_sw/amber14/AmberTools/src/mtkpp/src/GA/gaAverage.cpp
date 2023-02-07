/*!
   \file gaAverage.cpp
   \brief Performs chromosome averaging
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

#include "gaAverage.h"
#include "gaOperators.h"
#include "gaPopulation.h"
#include "gaIndividual.h"
#include "gaRanNumGen.h"

namespace MTKpp
{

// ============================================================
// Function : gaAverage()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaAverage::gaAverage(gaOperators *parent, std::vector<int> genePerChr,
                     int nChild, std::string aver, int chrPerInd)
{
    pParent = parent;
    this->nChild = nChild;
    this->aver = aver;
    this->chrPerInd = chrPerInd;
    this->genePerChr = genePerChr;
}

// ============================================================
// Function : ~gaAverage()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaAverage::~gaAverage() {}

// ============================================================
// Function : reproduce()
// ------------------------------------------------------------
// Main function of gaAverage
// ============================================================
bool gaAverage::average(gaIndividual* ind1, gaIndividual* ind2,
                        gaPopulation* pop, std::string name)
{
    bool successful = false;
    gaIndividual* newInd;
    gaIndividual* parent;

    // CHECK TO SEE IF THEY'VE REPRODUCED BEFORE
    // CHECK TO SEE IF THEY'VE REPRODUCED BEFORE
    int numKids = ind1->mated(ind2->getName());

    if (numKids >= this->nChild) {
      return successful;
    }

    ind1->addPartner(ind2->getName());
    ind2->addPartner(ind1->getName());

    if (this->chrPerInd == 1) {
      for (unsigned int i = 0; i < this->genePerChr.size(); i++) {
        if (this->genePerChr[i] == 1) {
          std::cout << "gaAverage:__sGene: number of chromosomes and genes is one" << std::endl;
          // exit here
        }
      }
    }

    int mainParent = int(ranNumBetweenZeroAndOne() * 2);
    if (mainParent == 0) {
       newInd = pop->addIndividual(ind1);
       parent = ind2;
    }
    else {
       newInd = pop->addIndividual(ind2);
       parent = ind1;
    }

    if (aver == "single-gene")    successful = this->__singleGene(newInd, parent);
    //if (aver == "multiple-gene")  successful = this->__multipleGene(newInd, parent);
    //if (aver == "random-locus")   successful = this->__randomLocus(newInd, parent);

    newInd->setName(name);
    newInd->setEvaluate(true);

    return successful;
}

// ============================================================
// Function : __singleGene()
// ------------------------------------------------------------
// single gene average per chromosome
//
//    {} == Individual
//    || == chromosome
//    [] == gene
//
// ind1: {|[aa][cc]||[ff][tt]|}
//                             --> {|[(a+b)/2 (a+b)/2][cc]||[(h+t)/2 (h+t)/2][tt]|}
// ind2: {|[bb][dd]||[hh][jj]|}
//
// ============================================================
bool gaAverage::__singleGene(gaIndividual* newInd, gaIndividual* parent)
{
    for (int i = 0; i < this->chrPerInd; i++) {
      int swapGene = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
      newInd->setGene(i, swapGene, parent);
    }
    return true;
}

} // MTKpp namespace

