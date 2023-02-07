/*! 
   \file gaWorld.cpp
   \brief Container for gaRegions
   \author Martin Peters

   This is the main class of the GA.
   Container for gaRegions

   $Date: 2007/09/14 10:28:10 $
   $Revision: 1.9 $

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

#include "gaWorld.h"
#include "gaRegion.h"
#include "gaOutput.h"
#include "gaVersion.h"

namespace MTKpp
{

// ============================================================
// Function : gaWorld()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaWorld::gaWorld(std::string name):itsName(name)
{
    itsName = name;
    pGaRegion = 0;
    chrPerInd = 0;
    outputFile = "";
    restartFile = "";
    convergFile = "";
    levelOfOutput = 1;
    status = 0;
    this->setAuthor(GA_AUTHOR);
    time (&startTime);
    pGaOutput = new gaOutput(this);
}

// ============================================================
// Function : ~gaWorld()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaWorld::~gaWorld() 
{
    delete pGaOutput;
}

// ============================================================
// Function :addRegion()
// ------------------------------------------------------------
// Add Region to World
// ============================================================
gaRegion* gaWorld::addRegion()
{
    pGaRegion = new gaRegion(this);

    itsRegionList.push_back(pGaRegion);

    return pGaRegion;
}

// ============================================================
// Function :delRegion()
// ------------------------------------------------------------
// Delete Region from World
// ============================================================
void gaWorld::delRegion(gaRegion* reg)
{
    regionIterator b = std::find(itsRegionList.begin(), itsRegionList.end(), reg);

    if (b != itsRegionList.end()) {
          itsRegionList.erase(b);
    }
}

// ============================================================
// Function : getRegion()
// ------------------------------------------------------------
// 
// ============================================================
gaRegion* gaWorld::getRegion(int number, bool id, bool index)
{
    if (id) {
      for (regionIterator c=itsRegionList.begin(); c != itsRegionList.end(); c++) {
        pGaRegion = *c;

        if (pGaRegion->getId() == number) {
          return pGaRegion;
        }
      }
    }
    if (index) {
      return itsRegionList[number]; /// NEED SOME ERROR HANDLING HERE
    }
    return NULL;
}

// ============================================================
// Function : getRegion()
// ------------------------------------------------------------
// 
// ============================================================
gaRegion* gaWorld::getRegion(std::string name)
{
    for (regionIterator c=itsRegionList.begin(); c != itsRegionList.end(); c++) {
      pGaRegion = *c;

      if (pGaRegion->getName() == name) {
        return pGaRegion;
      }
    }
    return NULL;
}

// ============================================================
// Function : getRegions()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<gaRegion*> gaWorld::getRegions()
{
    return this->itsRegionList;
}

// ============================================================
// Function : getNumRegions()
// ------------------------------------------------------------
// 
// ============================================================
int gaWorld::getNumRegions()
{
    if (!itsRegionList.empty()) {
      return itsRegionList.size();
    }
    return 0;
}

// ============================================================
// Function : setup()
// ------------------------------------------------------------
// 
// ============================================================
void gaWorld::setup(const int& chrPerInd, std::vector<int> iGenePerChr, std::vector<int> iGeneSizes,
                    const std::string outputFile, const std::string restartFile, const std::string convergFile, const int outputLevel)
{
#ifdef DEBUG
    std::cout << "gaWorld::setup" << std::endl;
#endif

    this->chrPerInd = chrPerInd;
    if (iGenePerChr.size() > 0) {
      this->genePerChr = iGenePerChr;
    }
    if (iGeneSizes.size() > 0) {
      this->geneSizes = iGeneSizes;
    }
    this->outputFile = outputFile;
    this->restartFile = restartFile;
    this->convergFile = convergFile;
    this->levelOfOutput = outputLevel;
    //pGaOutput->setProgramName(progName);
    pGaOutput->openOutputFile();

    if (this->chrPerInd == 1) {
      for (unsigned int i = 0; i < this->genePerChr.size(); i++) {
        if (this->genePerChr[i] == 1) {
          pGaOutput->writeError("Each chromosome must have more than one gene");
        }
      }
    }
}

// ============================================================
// Function : getGaOutput()
// ------------------------------------------------------------
// 
// ============================================================
gaOutput* gaWorld::getGaOutput()
{
    return this->pGaOutput;
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
// Initialize Region
// ============================================================
void gaWorld::initialize()
{
#ifdef DEBUG
    std::cout << "gaWorld::initialize" << std::endl;
#endif

    for (regionIterator c=itsRegionList.begin(); c != itsRegionList.end(); c++) {
      pGaRegion = *c;
      pGaRegion->initialize();
    }
}

// ============================================================
// Function : rank()
// ------------------------------------------------------------
// Rank Region
// ============================================================
void gaWorld::rank()
{
   std::cout << " DO SOME REGION RANKING HERE PLEASE " << std::endl;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
// 
// ============================================================
void gaWorld::setName(std::string n)
{
    this->itsName = n;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaWorld::getName()
{
    return this->itsName;
}

// ============================================================
// Function : getChrPerInd()
// ------------------------------------------------------------
// 
// ============================================================
int gaWorld::getChrPerInd()
{
    return this->chrPerInd;
}

// ============================================================
// Function : setGenePerChr()
// ------------------------------------------------------------
// 
// ============================================================
void gaWorld::setGenePerChr(std::vector<int> x)
{
#ifdef DEBUG
    std::cout << "gaWorld::setGenePerChr" << std::endl;
#endif
    this->genePerChr.clear();
    this->genePerChr.resize(0);
    for (unsigned int i = 0; i < x.size(); i++) {
     this->genePerChr.push_back(x[i]);
    }
}

// ============================================================
// Function : getGenePerChr()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<int> gaWorld::getGenePerChr()
{
    return this->genePerChr;
}

// ============================================================
// Function : setGeneSizes()
// ------------------------------------------------------------
// 
// ============================================================
void gaWorld::setGeneSizes(std::vector<int> x)
{
#ifdef DEBUG
    std::cout << "gaWorld::setGeneSizes" << std::endl;
#endif

    this->geneSizes.clear();
    this->geneSizes.resize(0);

    for (unsigned int i = 0; i < x.size(); i++) {
      this->geneSizes.push_back(x[i]);
    }
}

// ============================================================
// Function : getGeneSizes()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<int> gaWorld::getGeneSizes()
{
    return this->geneSizes;
}

// ============================================================
// Function : getOutputFileName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaWorld::getOutputFileName()
{
    return this->outputFile;
}

// ============================================================
// Function : getRestartFileName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaWorld::getRestartFileName()
{
    return this->restartFile;
}

// ============================================================
// Function : getConvergFileName()
// ------------------------------------------------------------
// 
// ============================================================
std::string gaWorld::getConvergFileName()
{
    return this->convergFile;
}

// ============================================================
// Function : getLevelOfOutput()
// ------------------------------------------------------------
// 
// ============================================================
int gaWorld::getLevelOfOutput()
{
    return this->levelOfOutput;
}

// ============================================================
// Function : setAuthor()
// ------------------------------------------------------------
// 
// ============================================================
void gaWorld::setAuthor(std::string name)
{
    this->authors.push_back(name);
}

// ============================================================
// Function : getAuthors()
// ------------------------------------------------------------
// 
// ============================================================
std::vector<std::string> gaWorld::getAuthors()
{
    return this->authors;
}

// ============================================================
// Function : setStatus()
// ------------------------------------------------------------
// 
// ============================================================
void gaWorld::setStatus(int status)
{
    this->status = status;
}

// ============================================================
// Function : getStatus()
// ------------------------------------------------------------
// 
// ============================================================
int gaWorld::getStatus()
{
    return this->status;
}

} // MTKpp namespace
