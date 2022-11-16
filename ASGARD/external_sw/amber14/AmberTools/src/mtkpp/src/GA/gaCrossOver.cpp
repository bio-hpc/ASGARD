/*! 
   \file gaCrossOver.cpp
   \brief Performs recombination of gaChromosomes
   \author Martin Peters

   CrossOver Functions Available:

   1) random single gene crossover per chromosome

   2) random multiple gene crossover per chromosome

   3) random loci crossover over all chromosomes

   $Date: 2007/09/14 10:28:10 $
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

#include "gaCrossOver.h"
#include "gaPopulation.h"
#include "gaIndividual.h"
#include "gaOutput.h"
#include "gaRanNumGen.h"

namespace MTKpp
{

// ============================================================
// Function : gaCrossOver()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaCrossOver::gaCrossOver(gaOutput* out, std::vector<int> genePerChr, 
                         int nChild, std::string crossover,int chrPerInd)
{
    pGaOutput = out;
    this->nChild = nChild;
    this->crossover = crossover;
    this->chrPerInd = chrPerInd;
    this->genePerChr = genePerChr;
}

// ============================================================
// Function : ~gaCrossOver()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaCrossOver::~gaCrossOver() {}

// ============================================================
// Function : reproduce()
// ------------------------------------------------------------
// Main function of gaCrossOver
// ============================================================
bool gaCrossOver::reproduce(gaIndividual* ind1, gaIndividual* ind2,
                            gaPopulation* pop, std::string name)
{
#ifdef DEBUG
    std::cout << "gaCrossOver::reproduce Ind1 = " << ind1->getName() << " Ind2 = "
              << ind2->getName() << " popName = " << pop->getName() << " newInd = "
              << name << std::endl;
#endif

    bool successful = false;
    gaIndividual* newInd;
    gaIndividual* parent;

    // CHECK TO SEE IF THEY'VE REPRODUCED BEFORE
    int numKids = ind1->mated(ind2->getName());

    if (numKids >= this->nChild) {
      return successful;
    }

    ind1->addPartner(ind2->getName());
    ind2->addPartner(ind1->getName());

    int mainParent = int(ranNumBetweenZeroAndOne() * 2);
    if (mainParent == 0) {
       newInd = pop->addIndividual(ind1);
       parent = ind2;
    }
    else {
       newInd = pop->addIndividual(ind2);
       parent = ind1;
    }
    if (crossover == "single-gene")    successful = this->__singleGene(newInd, parent);
    if (crossover == "multiple-gene")  successful = this->__multipleGene(newInd, parent);
    if (crossover == "random-locus")   successful = this->__randomLocus(newInd, parent);

    newInd->setName(name);
    newInd->setEvaluate(true);
    newInd->addParent(ind1->getName());
    newInd->addParent(ind2->getName());

    return successful;
}

// ============================================================
// Function : __singleGene()
// ------------------------------------------------------------
// single gene crossover per chromosome
//
//    {} == Individual
//    || == chromosome
//    [] == gene
//
// ind1: {|[aaaa][cccc]||[ffff][tttt]|}
//                                     --> {|[aaaa][dddd]||[hhhh][tttt]|}
// ind2: {|[bbbb][dddd]||[hhhh][jjjj]|}
//
// ============================================================
bool gaCrossOver::__singleGene(gaIndividual* newInd, gaIndividual* parent)
{
    for (int i = 0; i < this->chrPerInd; i++) {
      int swapGene = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
      newInd->setGene(i, swapGene, parent);
    }
    return true;
}

// ============================================================
// Function : __multipleGene()
// ------------------------------------------------------------
// multiple gene crossover per chromosome
//
//    {} == Individual
//    || == chromosome
//    [] == gene
//
// ind1: {|[aa][bb][cc]||[dd][ee][ff]|}
//                                     --> {|[aa][bb][xx]||[gg][ee][ii]|}
// ind2: {|[zz][yy][xx]||[gg][hh][ii]|}
//
// ============================================================
bool gaCrossOver::__multipleGene(gaIndividual* newInd, gaIndividual* parent)
{
    int numberCrossOvers = 0;

    for (int i = 0; i < this->chrPerInd; i++) {
      numberCrossOvers = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
      if (numberCrossOvers == 0) numberCrossOvers = 1;
      for (int j = 0; j < numberCrossOvers; j++) {
        int swapGene = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
        newInd->setGene(i, swapGene, parent);
      }
    }
    return true;
}

// ============================================================
// Function : __randomLocus()
// ------------------------------------------------------------
// random locus crossover per chromosome
//
//    {} == Individual
//    || == chromosome
//    [] == gene
//
//             |                  |
// ind1: {|[aa]|[bb][cc]||[dd][ee]|[ff]|}
//             |                  |      --> {|[aa][yy][xx]||[dd][ee][ii]}
// ind2: {|[zz]|[yy][xx]||[gg][hh]|[ii]|}
//             |                  |
//
// ============================================================
bool gaCrossOver::__randomLocus(gaIndividual* newInd, gaIndividual* parent) 
{
#ifdef DEBUG
    std::cout << "gaCrossOver::__randomLocus " << this->genePerChr.size() << " " << this->chrPerInd << std::endl;
#endif

    int locus = 0;

    for (int i = 0; i < this->chrPerInd; i++) {
      locus = int(ranNumBetweenZeroAndOne() * this->genePerChr[i]);
      for (int j = locus; j < this->genePerChr[i]; j++) {
        newInd->setGene(i, j, parent);
      }
    }
    return true;
}

} // MTKpp namespace

