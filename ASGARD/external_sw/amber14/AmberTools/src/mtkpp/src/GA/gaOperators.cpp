/*! 
   \file gaOperators.cpp
   \brief This class performs the GA operations including crossover
          mutation, averaging, etc.
   \author Martin Peters

   $Date: 2010/03/29 20:24:52 $
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


#include "gaOperators.h"

#include "gaWorld.h"
#include "gaRegion.h"
#include "gaPopulation.h"
#include "gaIndividual.h"
#include "gaSelection.h"
#include "gaCrossOver.h"
#include "gaMutate.h"
#include "gaAverage.h"
#include "gaRanNumGen.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : gaOperators()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaOperators::gaOperators(gaRegion *parent) {
    pGaOutput = parent->getParent()->getGaOutput();
    pParent = parent;

    this->nKeep = 0;
    this->nCrossover = 0;
    this->nMutate = 0;
    this->nAverage = 0;
    this->chreDiff = 0.1;
    mutantName = "_m";
}

// ============================================================
// Function : ~gaOperators()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaOperators::~gaOperators() {}

// ============================================================
// Function : setup()
// ------------------------------------------------------------
// Calculate nKeep, nCrossover, nMutate and nAverage
// ============================================================
void gaOperators::setup()
{
#ifdef DEBUG
    std::cout << "gaOperators::setup" << std::endl;
#endif

    pGaSelection = new gaSelection(this, pParent->selection, pParent->maxInds);

    pGaCrossOver = new gaCrossOver(pGaOutput, pParent->getParent()->getGenePerChr(),
                                   pParent->nChild, pParent->crossover,
                                   pParent->getParent()->getChrPerInd());

    pGaMutate = new gaMutate(this, pParent->getParent()->getGenePerChr(), 
                             pParent->mutate,
                             pParent->getParent()->getChrPerInd());

    pGaAverage = new gaAverage(this, pParent->getParent()->getGenePerChr(), 
                               pParent->nChild, pParent->average,
                               pParent->getParent()->getChrPerInd());

    int maxInds = pParent->maxInds;

    this->nKeep      = int(pParent->pKeep * maxInds) + 1;
    this->nCrossover = int(pParent->pCrossover * maxInds);
    this->nMutate    = int(pParent->pMutate * maxInds);
    this->nAverage   = maxInds - (this->nKeep + this->nCrossover + this->nMutate);

    this->chreDiff = pParent->chreDiff;

    if (pParent->selection == "semi-random") {
      pGaSelection->setGaussian(0, maxInds, pParent->selectionPressure);
    }
}

// ============================================================
// Function : keep()
// ------------------------------------------------------------
// nKeep individuals of the new generation are individuals from 
// the previous generation.
// ============================================================
void gaOperators::keep(gaPopulation* curPop, gaPopulation* prePop)
{
    gaIndividual* ind;
    gaIndividual* newInd;
    int i = 0;
    int kept = 0;
    std::vector<int> indKept;
    std::vector<int>::iterator result;

    // ELITISM
    if (pParent->elitism) {
      // get best individual and add to new generation
      int index = 0;
      ind = prePop->getBestIndividual(index);
      newInd = curPop->addIndividual(ind);
      indKept.push_back(index);

      std::vector<std::string> parents = ind->getParents();
      if (parents.size() == 2) {
        newInd->addParent(parents[0]);
        newInd->addParent(parents[1]);
      }

      std::vector<std::string> partners = ind->getPartners();
      for (unsigned int x = 0; x < partners.size(); x++) {
        newInd->addPartner(partners[x]);
      }

      newInd->setEvaluate(false);
      kept++;
    }
    while (kept < this->nKeep) {
      // select individual from population
      i = pGaSelection->select();
      result = std::find(indKept.begin(), indKept.end(), i);
      // if not selected already then add ---> bug if pop size is less than nKeep --> finite loop
      if (result == indKept.end()) {
        ind = prePop->getIndividual(i, 0, 1);
        if (ind) {
          newInd = curPop->addIndividual(ind);
          indKept.push_back(i);

          std::vector<std::string> parents = ind->getParents();
          if (parents.size() == 2) {
            newInd->addParent(parents[0]);
            newInd->addParent(parents[1]);
          }
          std::vector<std::string> partners = ind->getPartners();
          for (unsigned int x = 0; x < partners.size(); x++) {
            newInd->addPartner(partners[x]);
          }

          newInd->setEvaluate(false);
          kept++;
        }
        else {return;}
      }
    }
}

// ============================================================
// Function : crossover()
// ------------------------------------------------------------
// nCrossover individuals of the new generation are produced 
// by recombination
// Each time two parents are selected only one child is produced
// Two parents are allowed to produce nChild offspring
// ============================================================
void gaOperators::crossover(gaPopulation* curPop, gaPopulation* prePop)
{
    bool successful = false;
    std::string popName = curPop->getName();

    gaIndividual* ind1;
    gaIndividual* ind2;
    int i1 = 0;
    int i2 = 0;
    int crossover = 0;
    int failures = 0;

    while (crossover < this->nCrossover) {
      i1 = pGaSelection->select();
      i2 = pGaSelection->select();
      if (i1 == i2) {
        continue;
      }
      else {
        ind1 = prePop->getIndividual(i1, 0, 1);
        ind2 = prePop->getIndividual(i2, 0, 1);

        std::stringstream ss1;
        ss1 << curPop->getNumIndividuals();
        std::string nInd = ss1.str().c_str();
        std::string newIndName = popName + "_c" + nInd;
        successful = pGaCrossOver->reproduce(ind1, ind2, curPop, newIndName);
        if (successful) {
          crossover++;
        }
        else {
          failures++;
        }
        successful = false;
      }
      if (failures > 1000) { // maybe a percentage of the # of individuals
        return;
      }
    }
}

// ============================================================
// Function : mutate()
// ------------------------------------------------------------
// nMutate individuals of the new generation are produced
// by mutation
// ============================================================
void gaOperators::mutate(gaPopulation* curPop, gaPopulation* prePop)
{
    bool successful = false;
    std::string popName = curPop->getName();
    gaIndividual* ind1 = 0;
    gaIndividual* ind2 = 0;
    gaIndividual* newInd = 0;
    int i1 = 0;
    int i2 = 0;
    int mutate = 0;

    while (mutate < this->nMutate) {
      i1 = pGaSelection->select();
      i2 = pGaSelection->select();
      if (i1 == i2) {
        continue;
      }
      else {
        ind1 = prePop->getIndividual(i1, 0, 1);
        ind2 = prePop->getIndividual(i2, 0, 1);
        if (!ind1 or !ind2) {
          std::cout << " Error in gaOperator::mutate " << std::endl;
          throw MTKException( " Error in gaOperator::mutate ");
        }

        std::stringstream ss1;
        ss1 << curPop->getNumIndividuals();
        std::string nInd = ss1.str().c_str();
        std::string newIndName = popName + mutantName + nInd;

        successful = pGaCrossOver->reproduce(ind1, ind2, curPop, newIndName);
        if (successful) {
          newInd = curPop->getLastAdded();
          if (!newInd) {
            std::cout << " Error in gaOperator::mutate " << std::endl;
            throw MTKException( " Error in gaOperator::mutate ");
          }

          successful = pGaMutate->mutation(newInd);
#ifdef DEBUG
          std::cout << " gaOperators::mutate -> Adding " << newIndName << std::endl;
#endif
          mutate++;
        }
        successful = false;
      }
    }
}

// ============================================================
// Function : average()
// ------------------------------------------------------------
// nAverage individuals of the new generation are produced
// by averaging
// ============================================================
void gaOperators::average(gaPopulation* curPop, gaPopulation* prePop)
{
    if (this->nAverage == 0) return;
    bool successful = false;
    std::string popName = curPop->getName();

    gaIndividual* ind1;
    gaIndividual* ind2;
    int i1 = 0;
    int i2 = 0;
    int average = 0;

    while (average < this->nAverage) {
      i1 = pGaSelection->select();
      i2 = pGaSelection->select();
      if (i1 == i2) {
        continue;
      }
      else {
        ind1 = prePop->getIndividual(i1, 0, 1);
        ind2 = prePop->getIndividual(i2, 0, 1);

        std::stringstream ss1;
        ss1 << curPop->getNumIndividuals();
        std::string nInd = ss1.str().c_str();
        std::string newIndName = popName + "_a" + nInd;

        successful = pGaAverage->average(ind1, ind2, curPop, newIndName);
        if (successful) {
          average++;
        }
        successful = false;
      }
    }
}

// ============================================================
// Function : removeRedundant()
// ------------------------------------------------------------
// Removes clones from the population
// ============================================================
void gaOperators::removeRedundant(gaPopulation* curPop, gaPopulation* prePop)
{
    if (chreDiff == 0.0) {
      std::cout << " gaOperators::removeRedundant chreDiff is zero" << std::endl;
    }

    bool theSame = false;
    int clones = 0;
    gaIndividual* ind1;
    gaIndividual* ind2;
    gaIndividual* ind3;

    curPop->getAbsValues();
    std::vector<gaIndividual*> toBeDeleted;
    std::vector<gaIndividual*>::iterator vecIterator;

    for (int i = 0; i < curPop->getNumIndividuals(); i++) {
      ind1 = curPop->getIndividual(i, 0, 1);
      for (int j = i+1; j < curPop->getNumIndividuals(); j++) {
        ind2 = curPop->getIndividual(j, 0, 1);
        if ( (ind1->getAbsValue() - ind2->getAbsValue()) < this->chreDiff) {
          theSame = ind1->compare(ind2);
          if (theSame) {
            if (!ind1->getEvaluate()) {
              ind3 = ind2;
            }
            else {
              ind3 = ind1;
            }
            vecIterator = std::find(toBeDeleted.begin(), toBeDeleted.end(), ind3);
            if (vecIterator == toBeDeleted.end()) {
              toBeDeleted.push_back(ind3);
              clones++;
            }
            ind3 = 0;
          }
          theSame = false;
        }
      }
    }

    if (!toBeDeleted.empty()){
      std::sort( toBeDeleted.begin(), toBeDeleted.end() );
      for (unsigned int n = 0; n < toBeDeleted.size(); n++) {
        curPop->delIndividual(toBeDeleted[n]);
      }
    }
    //std::cout << " removeRedundant : number of individuals = " << curPop->getNumIndividuals() << std::endl;

    if (clones > 0) {
      int tempNmutate = this->nMutate;
      this->nMutate = clones;
      this->mutantName = "_r";
      this->mutate(curPop, prePop);
      this->mutantName = "_m";
      this->nMutate = tempNmutate;
    }
}

} // MTKpp namespace

