/*! 
   \file gaChromosome.cpp
   \brief Container for gaGenes
   \author Martin Peters

   Container for gaGenes

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

#include "gaChromosome.h"
#include "gaIndividual.h"
#include "gaGene.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : gaChromosome()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaChromosome::gaChromosome(gaIndividual *parent):pParent(parent) {}

// ============================================================
// Function : gaChromosome()
// ------------------------------------------------------------
// Copy Constructor for the class.
// ============================================================
gaChromosome::gaChromosome(gaChromosome* rhs)
{
    for (int i = 0; i < rhs->getNumGenes(); i++) {
      pGaGene = this->addGene(rhs->getGene(i,0,1));
    }
    pParent = rhs->getParent();
    itsId = pParent->getNumChromosomes()-1;
}

// ============================================================
// Function : ~gaChromosome()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaChromosome::~gaChromosome() {}

// ============================================================
// Function :addGene()
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaGene* gaChromosome::addGene()
{
    pGaGene = new gaGene(this);

    itsGeneList.push_back(pGaGene);

    pGaGene->setId(itsGeneList.size()-1);

    return pGaGene;
}

// ============================================================
// Function :addGene(gaGene*)
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaGene* gaChromosome::addGene(gaGene* rhs)
{
    pGaGene = new gaGene(rhs);

    itsGeneList.push_back(pGaGene);

    pGaGene->setId(itsGeneList.size()-1);

    return pGaGene;
}

// ============================================================
// Function :delGene()
// ------------------------------------------------------------
// Delete Gene from Chromosome
// ============================================================
void gaChromosome::delGene(gaGene* gen)
{
    geneIterator b = std::find(itsGeneList.begin(), itsGeneList.end(), gen);

    if (b != itsGeneList.end()) {
      itsGeneList.erase(b);
    }
}

// ============================================================
// Function : getGene()
// ------------------------------------------------------------
//
// ============================================================
gaGene* gaChromosome::getGene(int number, bool id, bool index)
{
    if (id) {
      for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
        pGaGene = *c;
        if (pGaGene->getId() == number) {
          return pGaGene;
        }
      }
    }
    if (index) {
      return itsGeneList[number]; // error handling
    }
    return NULL;
}

// ============================================================
// Function : setGene()
// ------------------------------------------------------------
//
// ============================================================
void gaChromosome::setGene(int number, bool id, bool index, gaChromosome* rhsChr)
{
    gaGene* rhsGene;
    if (id) {
      for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
        pGaGene = *c;
        if (pGaGene->getId() == number) {
          rhsGene = rhsChr->getGene(number, id, index);
          pGaGene->setGene(rhsGene);
        }
      }
    }
    if (index) {
      pGaGene = itsGeneList[number]; // error handling
      if (!pGaGene) {
        std::cout << " Error in gaChromosome::setGene " << std::endl;
      }
      rhsGene = rhsChr->getGene(number, id, index);
      if (!rhsGene) {
        std::cout << " Error in gaChromosome::setGene " << std::endl;
      }
      pGaGene->setGene(rhsGene);
    }
}

// ============================================================
// Function : mutateGene()
// ------------------------------------------------------------
//
// ============================================================
void gaChromosome::mutateGene(int number, bool id, bool index)
{
    if (id) {
      for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
        pGaGene = *c;
        if (pGaGene->getId() == number) {
          pGaGene->mutate();
        }
      }
    }
    if (index) {
      pGaGene = itsGeneList[number]; // error handling
      if (!pGaGene) {
        throw MTKException(" Error in gaChromosome::mutateGene ...  ");
      }
      pGaGene->mutate();
    }
}

// ============================================================
// Function : averageGene()
// ------------------------------------------------------------
//
// ============================================================
void gaChromosome::averageGene(int number, bool id, bool index, gaChromosome* rhsChr)
{
    gaGene* rhsGene;
    if (id) {
      for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
        pGaGene = *c;
        if (pGaGene->getId() == number) {
          rhsGene = rhsChr->getGene(number, id, index);
          pGaGene->average(rhsGene);
        }
      }
    }
    if (index) {
      pGaGene = itsGeneList[number]; // error handling
      rhsGene = rhsChr->getGene(number, id, index);
      pGaGene->average(rhsGene);
    }
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// Initialize
// ============================================================
void gaChromosome::initialize()
{
    for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
      pGaGene = *c;
      pGaGene->initialize();
    }
}

// ============================================================
// Function : getAbsValue()
// ------------------------------------------------------------
//
// ============================================================
double gaChromosome::getAbsValue()
{
    double a = 0.0;
    for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
      pGaGene = *c;
      a += pGaGene->getAbsValue();
    }
    return a;
}

// ============================================================
// Function : compare()
// ------------------------------------------------------------
//
// ============================================================
bool gaChromosome::compare(gaChromosome* rhs)
{
    bool cmp = false;
    gaGene* rhsGene;
    for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
      pGaGene = *c;
      rhsGene = rhs->getGene(pGaGene->getId(),0,1);
      cmp = pGaGene->compare( rhsGene );
      if (!cmp) {
        return cmp;
      }
    }
    return true;
}

// ============================================================
// Function : setId()
// ------------------------------------------------------------
// Set Internal Index
// ============================================================
void gaChromosome::setId(int n)
{
    itsId = n;
}

// ============================================================
// Function : getId()
// ------------------------------------------------------------
// Get Internal Index
// ============================================================
int gaChromosome::getId()
{
    return itsId;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
//
// ============================================================
void gaChromosome::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
//
// ============================================================
std::string gaChromosome::getName()
{
    return this->itsName;
}


// ============================================================
// Function : getParent()
// ------------------------------------------------------------
//
// ============================================================
gaIndividual* gaChromosome::getParent()
{
    return pParent;
}

// ============================================================
// Function : getNumGenes()
// ------------------------------------------------------------
//
// ============================================================
int gaChromosome::getNumGenes()
{
    if (!itsGeneList.empty()) {
      return itsGeneList.size();
    }
    return 0;
}

// ============================================================
// Function : getGeneticInformation()
// ------------------------------------------------------------
//
// ============================================================
std::vector<double> gaChromosome::getGeneticInformation()
{
    std::vector<double> geneticInfo;
    if (!itsGeneList.empty()) {
      for (geneIterator c = itsGeneList.begin(); c != itsGeneList.end(); c++) {
        pGaGene = *c;
        std::vector<double> geneInfo = pGaGene->getGeneticInformation();
        for (unsigned int i = 0; i < geneInfo.size(); i++) {
          geneticInfo.push_back(geneInfo[i]);
        }
      }
    }
    return geneticInfo;
}

} // MTKpp namespace

