/*! 
   \file gaGene.cpp
   \brief Class to handle genes
   \author Martin Peters

   Class to handle genes

   $Date: 2007/09/14 10:28:10 $
   $Revision: 1.6 $

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

#include "gaGene.h"
#include "gaWorld.h"
#include "gaRegion.h"
#include "gaPopulation.h"
#include "gaIndividual.h"
#include "gaChromosome.h"
#include "gaRanNumGen.h"

namespace MTKpp
{

// ============================================================
// Function : gaGene()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaGene::gaGene(gaChromosome *parent):pParent(parent) 
{
    curRegion = pParent->getParent()->getParent()->getParent();
}

// ============================================================
// Function : gaGene()
// ------------------------------------------------------------
// Copy Constructor for the class.
// ============================================================
gaGene::gaGene(gaGene* rhs)
{
    for (int i=0; i < rhs->getNumBits(); i++) {
      this->itsBits.push_back(rhs->getBit(i));
    }
    this->pParent = rhs->getParent();
    curRegion = pParent->getParent()->getParent()->getParent();
    this->itsId = pParent->getNumGenes()-1;
    this->itsName = rhs->getName();
}

// ============================================================
// Function : ~gaGene()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaGene::~gaGene() {}

// ============================================================
// Function : gaGene()
// ------------------------------------------------------------
// 
// ============================================================
void gaGene::setGene(gaGene* rhs)
{
    this->itsBits.clear();
    for (int i = 0; i < rhs->getNumBits(); i++) {
      this->itsBits.push_back(rhs->getBit(i));
    }
}

// ============================================================
// Function : addBit()
// ------------------------------------------------------------
// 
// ============================================================
void gaGene::addBit(double b)
{
    this->itsBits.push_back(b);
}

// ============================================================
// Function : getBit()
// ------------------------------------------------------------
// Get Bit
// ============================================================
double gaGene::getBit(const int &i)
{ // TO-DO  Error handling
    return itsBits[i];
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// Initialize
// ============================================================
void gaGene::initialize()
{
    bool setValues = false;

    // (works for 1 chr per inds with, N genes in the chr)
    if (curRegion->getParent()->getChrPerInd() == 1) {
      std::vector<int> vI = curRegion->getParent()->getGenePerChr();
      if (vI.size() == 1) {
        if (vI[0] == 1) {
          setValues = true;
        }
      }
    }

    if (setValues) {
      std::vector<double> possibleValues = curRegion->getGeneValues(this->itsId);

      bool d = true;
      double newValue = 0.0;
      if (possibleValues.size() == 1) {
        newValue = possibleValues[0];
        d = false;
      }
      while (d) {
        int rPos = int(ranNumBetweenZeroAndX(possibleValues.size()));
        if (!(possibleValues[rPos] == itsBits[0])) {
          newValue = possibleValues[rPos];
          d = false;
        }
      }
      itsBits[0] = newValue;
#ifdef DEBUG
      std::cout << " gaGene::initialize:old bit = " << itsBits[0]
                << " new bit = " <<  newValue << std::endl;
#endif
    }
    else {
      int rPos = int(ranNumBetweenZeroAndX(itsBits.size()));
#ifdef DEBUG
      std::cout << " gaGene::initialize:rPos = " << rPos << std::endl;
#endif
      int curChr = this->getParent()->getId();

      double test = 0.0;
      bool successful = false;

      while (!successful) {
        test = itsBits[rPos] + (2.0 * (ranNumBetweenZeroAndOne() - 0.5)
                                    * curRegion->getStepSize(curChr, this->itsId, rPos));
        if ( test < curRegion->getMinParameter(curChr, this->itsId, rPos)
             or test > curRegion->getMaxParameter(curChr, this->itsId, rPos)) {
          successful = false;
        }
        else {
#ifdef DEBUG
          std::cout << " gaGene::initialize:old bit = " << itsBits[rPos]
                    << " new bit = " << test << std::endl;
#endif
          itsBits[rPos] = test;
          successful = true;
        }
      }
    }
}

// ============================================================
// Function : getAbsValue()
// ------------------------------------------------------------
// absValue
// ============================================================
double gaGene::getAbsValue()
{
    double a = 0.0;
    for (unsigned int i = 0; i < this->itsBits.size(); i++) {
      a += std::abs(this->itsBits[i]);
    }
    return a;
}

// ============================================================
// Function : compare()
// ------------------------------------------------------------
// 
// ============================================================
bool gaGene::compare(gaGene* rhs)
{
    for (unsigned int i = 0; i < this->itsBits.size(); i++) {
      if ( (std::abs(this->itsBits[i]) - std::abs(rhs->getBit(i))) > curRegion->chreDiff) {
        return false;
      }
    }
    return true;
}

// ============================================================
// Function : mutate()
// ------------------------------------------------------------
// mutate
// ============================================================
void gaGene::mutate()
{
    bool setValues = false;

    // (works for 1 chr per inds with, N genes in the chr)
    if (curRegion->getParent()->getChrPerInd() == 1) {
      std::vector<int> vI = curRegion->getParent()->getGenePerChr();
      if (vI.size() == 1) {
        if (vI[0] == 1) {
          setValues = true;
        }
      }
    }

    if (setValues) {
      std::vector<double> possibleValues = curRegion->getGeneValues(this->itsId);

      bool d = true;
      double newValue = 0.0;
      if (possibleValues.size() == 1) {
        newValue = possibleValues[0];
        d = false;
      }
      while (d) {
        int rPos = int(ranNumBetweenZeroAndX(possibleValues.size()));
        if (!(possibleValues[rPos] == itsBits[0])) {
          newValue = possibleValues[rPos];
          d = false;
        }
      }
      itsBits[0] = newValue;
#ifdef DEBUG
      std::cout << " gaGene::initialize:old bit = " << itsBits[0]
                << " new bit = " <<  newValue << std::endl;
#endif
    }
    else {
      int rPos = int(ranNumBetweenZeroAndX(itsBits.size()));

      int curChr = this->getParent()->getId();

      double test = 0.0;
      bool successful = false;
      while (! successful) {
        test = itsBits[rPos] + 2.0 * (ranNumBetweenZeroAndOne()-0.5)
                                   * curRegion->getStepSize(curChr, this->itsId, rPos);
        if (test < curRegion->getMinParameter(curChr, this->itsId, rPos)
             or test > curRegion->getMaxParameter(curChr, this->itsId, rPos)) {
          successful = false;
        }
        else {
          itsBits[rPos] = test;
          successful = true;
        }
      }
    }
}

// ============================================================
// Function : average()
// ------------------------------------------------------------
// 
// ============================================================
void gaGene::average(gaGene* rhs)
{
    for (int i = 0; i < rhs->getNumBits(); i++) {
      this->itsBits[i] = (this->itsBits[i] + rhs->getBit(i)) / 2;
    }
}

// ============================================================
// Function : printToScreen()
// ------------------------------------------------------------
// 
// ============================================================
void gaGene::printToScreen()
{
    for (int i = 0; i < this->getNumBits();i++) {
      std::cout << this->getBit(i) << " ";
    }
    std::cout << " " <<std::endl; 
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
void gaGene::setId(int n)
{
    itsId = n;
}

// ============================================================
// Function : getId()
// ------------------------------------------------------------
// Get Internal Index
// ============================================================
int gaGene::getId()
{
    return itsId;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void gaGene::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaGene::getName()
{
    return this->itsName;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
// 
// ============================================================
gaChromosome* gaGene::getParent()
{
    return pParent;
}

// ============================================================
// Function : getNumBits()
// ------------------------------------------------------------
// Get Number of Bits
// ============================================================
int gaGene::getNumBits()
{
    if (!itsBits.empty()) {
      return itsBits.size();
    }
    return 0;
}

// ============================================================
// Function : getGeneticInformation()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<double> gaGene::getGeneticInformation()
{
    std::vector<double> geneticInfo;
    if (!itsBits.empty()) {
      for (unsigned int i = 0; i < itsBits.size(); i++) {
        geneticInfo.push_back(itsBits[i]);
      }
    }
    return geneticInfo;
}

} // MTKpp namespace
