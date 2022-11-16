/*! 
   \file gaIndividual.cpp
   \brief Container for gaChromosomes
   \author Martin Peters

   Container for gaChromosomes

   $Date: 2010/03/29 20:24:52 $
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

#include "gaIndividual.h"
#include "gaPopulation.h"
#include "gaChromosome.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : gaIndividual()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaIndividual::gaIndividual(gaPopulation* parent)
{
    this->pParent = parent;
    this->itsName = "";
    this->itsFitness = 0.0;
    this->bEvaluate = true;
}

// ============================================================
// Function : gaIndividual()
// ------------------------------------------------------------
// Copy Constructor for the class.
// ============================================================
gaIndividual::gaIndividual(gaIndividual* rhs, gaPopulation* parent)
{
    for (int i = 0; i < rhs->getNumChromosomes(); i++) {
      pGaChromosome = this->addChromosome(rhs->getChromosome(i,0,1));
    }
    this->pParent = parent;
    this->itsName = rhs->getName();
    this->itsFitness = rhs->getFitness();
    this->bEvaluate = true;
}

// ============================================================
// Function : ~gaIndividual()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaIndividual::~gaIndividual() {}

// ============================================================
// Function :addChromosome()
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaChromosome* gaIndividual::addChromosome()
{
    pGaChromosome = new gaChromosome(this);

    itsChromosomeList.push_back(pGaChromosome);

    pGaChromosome->setId(itsChromosomeList.size()-1);

    return pGaChromosome;
}

// ============================================================
// Function :addChromosome(gaChromosome*)
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaChromosome* gaIndividual::addChromosome(gaChromosome* rhs)
{
    pGaChromosome = new gaChromosome(rhs);

    itsChromosomeList.push_back(pGaChromosome);

    return pGaChromosome;
}

// ============================================================
// Function :delChromosome()
// ------------------------------------------------------------
// Delete Region from World
// ============================================================
void gaIndividual::delChromosome(gaChromosome* chr)
{
    chromosomeIterator b = std::find(itsChromosomeList.begin(), itsChromosomeList.end(), chr);

    if (b != itsChromosomeList.end()) {
      itsChromosomeList.erase(b);
    }
}

// ============================================================
// Function : getChromosome()
// ------------------------------------------------------------
// 
// ============================================================
gaChromosome* gaIndividual::getChromosome(int number, bool id, bool index)
{
    if (id) {
      for (chromosomeIterator c=itsChromosomeList.begin(); c != itsChromosomeList.end(); c++) {
        pGaChromosome = *c;
        if (pGaChromosome->getId() == number) {
          return pGaChromosome;
        }
      }
    }
    if (index) {
      return itsChromosomeList[number];
    }
    return NULL;
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// Initialize 
// ============================================================
void gaIndividual::initialize()
{
#ifdef DEBUG
    std::cout << "gaIndividual::initialize" << std::endl;
#endif

    for (chromosomeIterator c = itsChromosomeList.begin(); c != itsChromosomeList.end(); c++) {
      pGaChromosome = *c;
      if (!pGaChromosome) {
        std::cout << " Error in gaIndividual::initialize ... exiting " << std::endl;
        throw MTKException(" Error in gaIndividual::initialize ... exiting ");
      }
      pGaChromosome->initialize();
    }
}

// ============================================================
// Function : getAbsValue()
// ------------------------------------------------------------
// 
// ============================================================
double gaIndividual::getAbsValue()
{
    this->abs_value = 0.0;
    for (chromosomeIterator c=itsChromosomeList.begin(); c != itsChromosomeList.end(); c++) {
      pGaChromosome = *c;
      this->abs_value +=  pGaChromosome->getAbsValue();
    }
    return this->abs_value;
}

// ============================================================
// Function : compare()
// ------------------------------------------------------------
// 
// ============================================================
bool gaIndividual::compare(gaIndividual* rhs)
{
    bool c = false;
    for (unsigned int i = 0; i < this->itsChromosomeList.size(); i++) {
      c = this->itsChromosomeList[i]->compare(rhs->getChromosome(i,0,1));
      if (!c) {
        return c;
      }
    }
    return true;
}

// ============================================================
// Function : printToScreen()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::printToScreen()
{
    std::cout << this->itsName << " " << this->itsFitness << std::endl;
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
void gaIndividual::setId(int n)
{
    itsId = n;
}

// ============================================================
// Function : getId()
// ------------------------------------------------------------
// Get Internal Index
// ============================================================
int gaIndividual::getId()
{
    return itsId;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaIndividual::getName()
{
    return this->itsName;
}

// ============================================================
// Function : setFitness()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::setFitness(double score)
{
    this->itsFitness = score;
}

// ============================================================
// Function : getFitness()
// ------------------------------------------------------------
// 
// ============================================================
double gaIndividual::getFitness()
{
    return this->itsFitness;
}

// ============================================================
// Function : setEvaluate()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::setEvaluate(bool b)
{
    this->bEvaluate = b;
}

// ============================================================
// Function : getEvaluate()
// ------------------------------------------------------------
// 
// ============================================================
bool gaIndividual::getEvaluate()
{
    return this->bEvaluate;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
// 
// ============================================================
gaPopulation* gaIndividual::getParent()
{
    return pParent;
}

// ============================================================
// Function : getNumChromosomes()
// ------------------------------------------------------------
// Get Internal Index
// ============================================================
int gaIndividual::getNumChromosomes()
{
    if (!itsChromosomeList.empty()) {
      return itsChromosomeList.size();
    }
    return 0;
}

// ============================================================
// Function : getGeneticInformation()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<double> gaIndividual::getGeneticInformation()
{
    std::vector<double> geneticInfo;
    if (!itsChromosomeList.empty()) {
      for (chromosomeIterator c=itsChromosomeList.begin(); c != itsChromosomeList.end(); c++) {
        pGaChromosome = *c;
        std::vector<double> chrInfo = pGaChromosome->getGeneticInformation();
        for (unsigned int i = 0; i < chrInfo.size(); i++) {
          geneticInfo.push_back(chrInfo[i]);
        }
      }
    }
    return geneticInfo;
}

// ============================================================
// Function : setGene()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::setGene(int chr, int gen, gaIndividual* ind)
{
    gaChromosome* lhsChr = this->getChromosome(chr, 0, 1);
    gaChromosome* rhsChr = ind->getChromosome(chr, 0, 1);

    lhsChr->setGene(gen, 0, 1, rhsChr);
}

// ============================================================
// Function : mutateGene()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::mutateGene(int chr, int gen)
{
    gaChromosome* pChr = this->getChromosome(chr, 0, 1);
    if (!pChr) {
      std::cout << " gaIndividual::mutateGene ... exiting" << std::endl;
      throw MTKException(" gaIndividual::mutateGene ... exiting");
    }
    pChr->mutateGene(gen, 0, 1);
}

// ============================================================
// Function : averageGene()
// ------------------------------------------------------------
// 
// ============================================================
void gaIndividual::averageGene(int chr, int gen, gaIndividual* ind)
{
    gaChromosome* lhsChr = this->getChromosome(chr, 0, 1);
    gaChromosome* rhsChr = ind->getChromosome(chr, 0, 1);

    lhsChr->averageGene(gen, 0, 1, rhsChr);
}

// ============================================================
// Function : addParent()
// ------------------------------------------------------------
// Add parent to parent list
// ============================================================
void gaIndividual::addParent(std::string s)
{
    this->itsParents.push_back(s);
}

// ============================================================
// Function : getParents()
// ------------------------------------------------------------
// Add parent to parent list
// ============================================================
std::vector<std::string> gaIndividual::getParents()
{
    return this->itsParents;
}

// ============================================================
// Function : addPartner()
// ------------------------------------------------------------
// Add partner to partner list
// ============================================================
void gaIndividual::addPartner(std::string s)
{
    this->itsPartners.push_back(s);
}

// ============================================================
// Function : getPartners()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<std::string> gaIndividual::getPartners()
{
    return this->itsPartners;
}

// ============================================================
// Function : numPartners()
// ------------------------------------------------------------
// Returns number of partners
// ============================================================
int gaIndividual::numPartners()
{
    return this->itsPartners.size();
}

// ============================================================
// Function : mated()
// ------------------------------------------------------------
// Number of times it has mated with s
// ============================================================
int gaIndividual::mated(std::string s)
{
    if (this->itsPartners.size() > 0) {
      return count(this->itsPartners.begin(), this->itsPartners.end(), s);
    }
    return 0;
}

} // MTKpp namespace
