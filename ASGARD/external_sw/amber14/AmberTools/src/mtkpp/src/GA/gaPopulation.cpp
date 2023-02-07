/*! 
   \file gaPopulation.cpp
   \brief Container for gaIndividuals
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

#include "gaPopulation.h"
#include "gaRegion.h"
#include "gaIndividual.h"
#include "gaRanNumGen.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : gaPopulation()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaPopulation::gaPopulation(gaRegion *parent):pParent(parent) {}

// ============================================================
// Function : ~gaPopulation()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaPopulation::~gaPopulation() {}

// ============================================================
// Function :addIndividual()
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaIndividual* gaPopulation::addIndividual()
{
    pGaIndividual = new gaIndividual(this);
    if (!pGaIndividual) {
      std::cout << " Error in gaPopulation::addIndividual ... exiting " << std::endl;
      throw MTKException(" Error in gaPopulation::addIndividual ... exiting ");
    }
    this->itsIndividualList.push_back(pGaIndividual);

    pGaIndividual->setId(this->itsIndividualList.size());

    return pGaIndividual;
}

// ============================================================
// Function :addIndividual(gaIndividual)
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaIndividual* gaPopulation::addIndividual(gaIndividual* rhs)
{
    pGaIndividual = new gaIndividual(rhs, this);
    if (!pGaIndividual) {
      std::cout << " Error in gaPopulation::addIndividual ... exiting " << std::endl;
      throw MTKException(" Error in gaPopulation::addIndividual ... exiting ");
    }
    itsIndividualList.push_back(pGaIndividual);

    pGaIndividual->setId(this->itsIndividualList.size());

    return pGaIndividual;
}

// ============================================================
// Function :delIndividual()
// ------------------------------------------------------------
// Delete Individual from Population
// ============================================================
void gaPopulation::delIndividual(gaIndividual* ind)
{
    individualIterator b = std::find(itsIndividualList.begin(), itsIndividualList.end(), ind);

    if (b != itsIndividualList.end()) {
      itsIndividualList.erase(b);
    }
}

// ============================================================
// Function : getIndividual()
// ------------------------------------------------------------
// 
// ============================================================
gaIndividual* gaPopulation::getIndividual(int number, bool id, bool index)
{
    if (id) {
      for (individualIterator c=itsIndividualList.begin(); c != itsIndividualList.end(); c++) {
        pGaIndividual = *c;
        if (pGaIndividual->getId() == number) {
          return pGaIndividual;
        }
      }
    }
    if (index) {
      return itsIndividualList[number];
    }
    return NULL;
}

// ============================================================
// Function : getIndividual()
// ------------------------------------------------------------
// 
// ============================================================
gaIndividual* gaPopulation::getIndividual(std::string name)
{
    for (individualIterator c=itsIndividualList.begin(); c != itsIndividualList.end(); c++) {
      pGaIndividual = *c;
      if (pGaIndividual->getName() == name) {
        return pGaIndividual;
      }
    }
    return NULL;
}

// ============================================================
// Function : getIndividuals()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<gaIndividual*> gaPopulation::getIndividuals()
{
    return this->itsIndividualList;
}

// ============================================================
// Function : getLastAdded()
// ------------------------------------------------------------
// 
// ============================================================
gaIndividual* gaPopulation::getLastAdded()
{
    if (this->itsIndividualList.size() > 0) {
      return this->itsIndividualList.back();
    }
    return NULL;
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// - Initializes population's chromosomes
//
// - If the population size is less that numInds then a
//   population is build from the first individual.
//
// - If population size is less than numInds but not equal 
//   to one, a population is build from the individuals provided.
//
// - If popluation size is equal to numInds then 
//   no initialization is done.
// ============================================================
void gaPopulation::initialize()
{
#ifdef DEBUG
    std::cout << "gaPopulation::initialize" << std::endl;
#endif

    int maxNumInds = pParent->maxInds; // define this at the world or region level

    int curNumInds = this->getNumIndividuals();

    if (curNumInds > maxNumInds) {
      std::cout << "Current Number of individuals is greater than the maximum allowed\n" << std::endl;
      std::cout << "\tIncrease maxInds in your input file." << std::endl;
      return;
    }

    if (curNumInds < maxNumInds) {
      // Initialize new population from Best individual
      if (curNumInds == 1) {
        bestIndividual = itsIndividualList[0];

        for (int i = 1; i < maxNumInds; i++) {
          pGaIndividual = this->addIndividual(bestIndividual);
          if (!pGaIndividual) {
            std::cout << " Error in gaPopulation::initialize ... exiting " << std::endl;
            throw MTKException(" Error in gaPopulation::initialize ... exiting ");
          }
          pGaIndividual->initialize();
          pGaIndividual->setId(i+1);

          std::stringstream ss1;
          ss1 << i;
          std::string str_ind = ss1.str().c_str(); 
          std::string iName = this->itsName + "_i" + str_ind;
          pGaIndividual->setName(iName);
        }
      }
      // Initialize new population from X individuals randomly
      else {
        for (int i = curNumInds; i < maxNumInds; i++) {
          pGaIndividual = this->addIndividual();
          int ranIndex = int(ranNumBetweenZeroAndX(curNumInds));
          bestIndividual = this->getIndividual(ranIndex,0,1);
          pGaIndividual = this->addIndividual(bestIndividual);
          pGaIndividual->initialize();
          pGaIndividual->setId(i+1);

          std::stringstream ss1;
          ss1 << i;
          std::string str_ind = ss1.str().c_str(); 
          std::string iName = this->itsName + "_i"+str_ind;
          pGaIndividual->setName(iName);
        }
      }
    }
}

// ============================================================
// Function : getAbsValue()
// ------------------------------------------------------------
// 
// ============================================================
void gaPopulation::getAbsValues()
{
    double a;
    for (int i = 0; i < this->getNumIndividuals(); i++) {
      pGaIndividual = this->getIndividual(i,0,1);
      a = pGaIndividual->getAbsValue();
    }
}

// ============================================================
// Function : rank()
// ------------------------------------------------------------
// Rank Individuals
// ============================================================
void gaPopulation::rank()
{
#ifdef DEBUG
    puts("gaPopulation::rank");
#endif

    if (pParent->funcDir == "minimize") {
      std::sort(this->itsIndividualList.begin(), this->itsIndividualList.end(), gaIndividual::less);
    }
    else {
      std::sort(this->itsIndividualList.begin(), this->itsIndividualList.end(), gaIndividual::greater);
    }

#ifdef DEBUG
    for (individualIterator c = this->itsIndividualList.begin(); c != this->itsIndividualList.end(); c++) {
      pGaIndividual = *c;
      std::cout << " " << pGaIndividual->getName() << " " << pGaIndividual->getFitness() << std::endl;
    }
#endif
}

// ============================================================ 
// ===                                                      === 
// ===          s e t / g e t  f u n c t i o n s            === 
// ===                                                      === 
// ============================================================ 

// ============================================================
// Function : setId()
// ------------------------------------------------------------
// Set Internal Index
// ============================================================
void gaPopulation::setId(int n)
{
    itsId = n;
}

// ============================================================
// Function : getId()
// ------------------------------------------------------------
// Get Internal Index
// ============================================================
int gaPopulation::getId()
{
    return itsId;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void gaPopulation::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaPopulation::getName()
{
    return this->itsName;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
// 
// ============================================================
gaRegion* gaPopulation::getParent()
{
    return this->pParent;
}

// ============================================================
// Function : getNumIndividuals()
// ------------------------------------------------------------
// 
// ============================================================
int gaPopulation::getNumIndividuals()
{
    if (!itsIndividualList.empty()) {
      return itsIndividualList.size();
    }
    return 0;
}

// ============================================================
// Function : getBestIndividual()
// ------------------------------------------------------------
// Returns best gaIndividual based on fitness
// ============================================================
gaIndividual* gaPopulation::getBestIndividual(int &i)
{
    bestIndividual = itsIndividualList[0];
    double fit = bestIndividual->getFitness();
    int t = 0;
    for (individualIterator c = itsIndividualList.begin(); c != itsIndividualList.end(); c++) {
      pGaIndividual = *c;
      if (pParent->funcDir == "minimize") {
        if (pGaIndividual->getFitness() < fit) {
          bestIndividual = pGaIndividual;
          i = t;
        }
      }
      else {
        if (pGaIndividual->getFitness() > fit) {
          bestIndividual = pGaIndividual; 
          i = t;
        }
      }
      t++;
    }
    return bestIndividual;
}

} // MTKpp namespace

