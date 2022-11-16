/*! 
   \file gaRegion.cpp
   \brief Container for gaPopulations
   \author Martin Peters

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

#include "gaRegion.h"
#include "gaWorld.h"
#include "gaIndividual.h"
#include "gaPopulation.h"
#include "gaOperators.h"

namespace MTKpp
{

// ============================================================
// Function : gaRegion()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaRegion::gaRegion(gaWorld *parent):pParent(parent)
{
    maxInds = 0;
    seed = 0;
    popKeep = 0;

    pKeep = 0.0;
    pCrossover = 0.0;
    pMutate = 0.0;
    pAverage = 0.0;

    maxGens = 0;
    nChild = 0;

    elitism = 0;
    selectionPressure = 0.0;
    chreDiff = 0.0;

    selection = "";
    crossover = "";
    mutate = "";
    average = "";

    itsName = "";
    itsId = 0;
    curPop = 0;

    pGaOperators = new gaOperators(this);
}

// ============================================================
// Function : ~gaRegion()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaRegion::~gaRegion() {}

// ============================================================
// Function : setup()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setup(const int& maxInds, const int& seed, const int& maxGens,
               const int& elitism, const double& selectionPressure, const int& curPop,
               const int& popKeep, const int& nChild, const double& chreDiff, 
               const double& pKeep, const double& pCrossover, const double& pMutate, const double& pAverage,
               const std::string& selection, const std::string& crossover,
               const std::string& mutate, const std::string& average, const std::string& funcDir,
               std::vector<double> maxParameters, std::vector<double> minParameters,
               std::vector<double> stepSize)
{
    this->maxInds = maxInds;
    this->seed = seed;
    this->maxGens = maxGens;
    this->elitism = elitism;
    this->selectionPressure = selectionPressure;
    this->curPop = curPop;
    this->popKeep = popKeep;
    this->nChild = nChild;
    this->chreDiff = chreDiff;
    this->pKeep = pKeep;
    this->pCrossover = pCrossover;
    this->pMutate = pMutate;
    this->pAverage = pAverage;
    this->selection = selection;
    this->crossover = crossover;
    this->mutate = mutate;
    this->average = average;
    this->funcDir = funcDir;
    this->maxParameters = maxParameters;
    this->minParameters = minParameters;
    this->stepSize = stepSize;
}

// ============================================================
// Function :addPopulation()
// ------------------------------------------------------------
// Add Population to Region
// ============================================================
gaPopulation* gaRegion::addPopulation()
{
    int id = 1;
    std::string name = "pop";

    pGaPopulation = new gaPopulation(this);

    if (this->itsPopulationList.size() > 0) {
      gaPopulation* pCurPop = this->itsPopulationList[0];
      id = this->itsPopulationList.size();
      //id = pCurPop->getId() + 1;
      name = pCurPop->getName();

      std::stringstream ss1;
      ss1 << id;
      std::string popId = ss1.str().c_str();
      name = name + popId;
    }

    pGaPopulation->setId(id);
    pGaPopulation->setName(name);
    this->itsPopulationList.push_back(pGaPopulation);
    return pGaPopulation;
}

// ============================================================
// Function :delPopulation()
// ------------------------------------------------------------
// Delete Population from Region
// ============================================================
void gaRegion::delPopulation(gaPopulation* pop)
{
    populationIterator b = std::find(itsPopulationList.begin(), itsPopulationList.end(), pop);

    if (b != itsPopulationList.end()) {
      itsPopulationList.erase(b);
    }
}

// ============================================================
// Function : getPopulation()
// ------------------------------------------------------------
// Return Populations by number or id
// ============================================================
gaPopulation* gaRegion::getPopulation(const int &number, bool id, bool index)
{
    if (id) {
      for (populationIterator c=itsPopulationList.begin(); c != itsPopulationList.end(); c++) {
        pGaPopulation = *c;
        if (pGaPopulation->getId() == id) {
          return pGaPopulation;
        }
      }
    }
    if (index) {
      return itsPopulationList[number];/// NEED SOME ERROR HANDLING HERE
    }
    return NULL;
}

// ============================================================
// Function : getPopulations()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<gaPopulation*> gaRegion::getPopulations()
{
    return this->itsPopulationList;
}

// ============================================================
// Function : getNumPopulations()
// ------------------------------------------------------------
// Returns the number of populations stored in the region
// ============================================================
int gaRegion::getNumPopulations()
{
    if (!itsPopulationList.empty()) {
      return itsPopulationList.size();
    }
    return 0;
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// This only initializes the current population.
// ============================================================
void gaRegion::initialize()
{
#ifdef DEBUG
    std::cout << "gaRegion::initialize" << std::endl;
#endif

    gaPopulation* pCurGaPopulation = this->itsPopulationList[0];

    pCurGaPopulation->initialize();

    pGaOperators->setup();
}

// ============================================================
// Function : nextGeneration()
// ------------------------------------------------------------
// Consider the following parameters:
//
//      Keep       50%
//      Crossover  20%
//      Mutate     18%
//      Average    12%
//
//   50% of the individuals are eliminated based on fitness,
//   then these empty slot are generated from mutation of
//   individuals (18%), or cross-breeding (20%) and averaging
//   (12%) of pairs.
//
// ============================================================
void gaRegion::nextGeneration()
{
    typedef std::vector<gaIndividual*>::iterator individualIterator;

#ifdef DEBUG
    std::cout << "gaRegion::nextGeneration #pops = "
              << this->getNumPopulations() << std::endl;
#endif

    // Add a new generation
    pGaPopulation = this->addPopulation();
    int numPop = this->itsPopulationList.size();
    pPreviousGaPop = itsPopulationList[numPop-2];

#ifdef DEBUG
    std::cout << "Number of Individuals :" 
              << pGaPopulation->getNumIndividuals() << " "
              << pPreviousGaPop->getNumIndividuals() << std::endl;

    std::cout << " gaRegion::nextGeneration:keep" << std::endl;
#endif

    // Keep
    pGaOperators->keep(pGaPopulation, pPreviousGaPop);

#ifdef DEBUG
    std::cout << "After Keep: Number of Individuals :"
              << pGaPopulation->getNumIndividuals() << " "
              << pPreviousGaPop->getNumIndividuals() << std::endl;

    std::vector<gaIndividual*> indList = pGaPopulation->getIndividuals();
    for (individualIterator c = indList.begin(); c != indList.end(); c++) {
      gaIndividual* pGaIndividual = *c;
      std::cout << pGaIndividual->getName() << " ";
    }
    std::cout << " " << std::endl;
    std::cout << " gaRegion::nextGeneration:crossover" << std::endl;
#endif

    // Cross Over (recombination)
    pGaOperators->crossover(pGaPopulation, pPreviousGaPop);

#ifdef DEBUG
    std::cout << "After Crossover: Number of Individuals :"
              << pGaPopulation->getNumIndividuals() << " " 
              << pPreviousGaPop->getNumIndividuals() << std::endl;

    indList = pGaPopulation->getIndividuals();
    for (individualIterator c = indList.begin(); c != indList.end(); c++) {
      gaIndividual* pGaIndividual = *c;
      std::cout << pGaIndividual->getName() << " ";
    }
    std::cout << " " << std::endl;
    std::cout << " gaRegion::nextGeneration:mutate" << std::endl;
#endif

    // Mutation
    pGaOperators->mutate(pGaPopulation, pPreviousGaPop);

#ifdef DEBUG
    std::cout << "After Mutate: Number of Individuals :"
              << pGaPopulation->getNumIndividuals() << " "
              << pPreviousGaPop->getNumIndividuals() << std::endl;

    indList = pGaPopulation->getIndividuals();
    for (individualIterator c = indList.begin(); c != indList.end(); c++) {
      gaIndividual* pGaIndividual = *c;
      std::cout << pGaIndividual->getName() << " ";
    }
    std::cout << " " << std::endl;
    std::cout << " gaRegion::nextGeneration:average" << std::endl;
#endif

    // Average
    pGaOperators->average(pGaPopulation, pPreviousGaPop);

#ifdef DEBUG
    std::cout << "After Average: Number of Individuals :"
              << pGaPopulation->getNumIndividuals() << " "
              << pPreviousGaPop->getNumIndividuals() << std::endl;

    indList = pGaPopulation->getIndividuals();
    for (individualIterator c = indList.begin(); c != indList.end(); c++) {
      gaIndividual* pGaIndividual = *c;
      std::cout << pGaIndividual->getName() << " ";
    }
    std::cout << " " << std::endl;
    std::cout << " gaRegion::nextGeneration:remove redundant" << std::endl;
#endif

    // Remove Redundant and Fill up remaining slots with mutants
    pGaOperators->removeRedundant(pGaPopulation, pPreviousGaPop);

#ifdef DEBUG
    std::cout << "After removeRedundant: Number of Individuals :"
              << pGaPopulation->getNumIndividuals() << " "
              << pPreviousGaPop->getNumIndividuals() << std::endl;

    indList = pGaPopulation->getIndividuals();
    for (individualIterator c = indList.begin(); c != indList.end(); c++) {
      gaIndividual* pGaIndividual = *c;
      std::cout << pGaIndividual->getName() << " ";
    }
    std::cout << " " << std::endl;
#endif
/*
    if (itsPopulationList.size() > static_cast<unsigned int> (this->popKeep)) {
      gaPopulation* ft = itsPopulationList.front();
      delPopulation(ft);
    }
*/
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
void gaRegion::setId(int n)
{
    itsId = n;
}

// ============================================================
// Function : getId()
// ------------------------------------------------------------
// Get Internal Index
// ============================================================
int gaRegion::getId()
{
    return itsId;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaRegion::getName()
{
    return this->itsName;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
// 
// ============================================================
gaWorld* gaRegion::getParent()
{
    return pParent;
}

// ============================================================
// Function : setFitness()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setFitness(double d)
{
    this->popConv.push_back(d);
}

// ============================================================
// Function : getFitness()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<double> gaRegion::getFitness()
{
    return this->popConv;
}

// ============================================================
// Function : setMinParameters()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setMinParameters(std::vector<double> m)
{
#ifdef DEBUG
    std::cout << "gaRegion::setMinParameters" << std::endl;
#endif
    this->minParameters = m;
}

// ============================================================
// Function : getMinParameter()
// ------------------------------------------------------------
// 
// ============================================================
double gaRegion::getMinParameter(const int& chr, const int& gen, const int& pos)
{
    int index1 = 0;
    int index2 = 0;
    for (int i = 0; i < pParent->chrPerInd; i++) {
      int g = pParent->genePerChr[i];
      for (int j = 0; j < g; j++) {
        int l = pParent->geneSizes[index1];
        index1++;
        for (int k = 0; k < l; k++) {
          if (i == chr and j == gen and k == pos) {
            return this->minParameters[index2];
          }
          index2++;
        }
      }
    }
    std::cout << "gaRegion::getMinParameter: OUT OF RANGE " << std::endl;
    return 0.0;
}

// ============================================================
// Function : setMaxParameters()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setMaxParameters(std::vector<double> m)
{
    this->maxParameters = m;
}

// ============================================================
// Function : getMaxParameter()
// ------------------------------------------------------------
// 
// ============================================================
double gaRegion::getMaxParameter(const int& chr, const int& gen, const int& pos)
{
    int index1 = 0;
    int index2 = 0;
    for (int i = 0; i < pParent->chrPerInd; i++) {
      int g = pParent->genePerChr[i];
      for (int j = 0; j < g; j++) {
        int l = pParent->geneSizes[index1];
        index1++;
        for (int k = 0; k < l; k++) {
          if (i == chr and j == gen and k == pos) {
            return this->maxParameters[index2];
          }
          index2++;
        }
      }
    }
    std::cout << "gaRegion::getMaxParameter: OUT OF RANGE " << std::endl;
    return 0.0;
}

// ============================================================
// Function : setStepSizes()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setStepSizes(std::vector<double> s)
{
    this->stepSize = s;
}

// ============================================================
// Function : getStepSize()
// ------------------------------------------------------------
// 
// ============================================================
double gaRegion::getStepSize(const int& chr, const int& gen, const int& pos)
{
    int index1 = 0;
    int index2 = 0;
    for (int i = 0; i < pParent->chrPerInd; i++) {
      int g = pParent->genePerChr[i];
      for (int j = 0; j < g; j++) {
        int l = pParent->geneSizes[index1];
        index1++;
        for (int k = 0; k < l; k++) {
          if (i == chr and j == gen and k == pos) {
            return this->stepSize[index2];
          }
          index2++;
        }
      }
    }
    std::cout << "gaRegion::getStepSize: OUT OF RANGE " << std::endl;
    return 0.0;
}

// ============================================================
// Function : getMaxGens()
// ------------------------------------------------------------
// 
// ============================================================
int gaRegion::getMaxGens()
{
    return this->maxGens;
}

// ============================================================
// Function : setGeneValues()
// ------------------------------------------------------------
// 
// ============================================================
void gaRegion::setGeneValues(std::vector<std::vector<double> > gVs)
{
    this->geneValues = gVs;
}

// ============================================================
// Function : getGeneValues()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<double> gaRegion::getGeneValues(const int& gen)
{
    return this->geneValues[gen];
}

} // MTKpp namespace
