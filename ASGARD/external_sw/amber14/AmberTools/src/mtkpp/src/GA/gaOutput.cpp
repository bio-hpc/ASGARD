/*! 
   \file gaOutput.cpp
   \brief Class to handle the GA output
   \author Martin Peters

   $Date: 2010/03/29 20:24:52 $
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
#include <string.h>
#include "gaOutput.h"
#include "gaWorld.h"
#include "gaRegion.h"
#include "gaPopulation.h"
#include "gaIndividual.h"
#include "gaChromosome.h"
#include "gaGene.h"
#include "gaVersion.h"

#include "Utils/constants.h"
#include "Utils/copyright.h"

namespace MTKpp
{

// ============================================================
// Function : gaOutput()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaOutput::gaOutput(gaWorld *w):myWorld(w) 
{
    this->programName = "MTK++::GA";
}

// ============================================================
// Function : ~gaOutput()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaOutput::~gaOutput() {}

// ============================================================
// Function : openOutputFile()
// ------------------------------------------------------------
// Open output file
// ============================================================
void gaOutput::openOutputFile()
{
    std::string filename = this->myWorld->getOutputFileName();

    this->outputFileStream.open(filename.c_str());

    if (!this->outputFileStream) {
      std::cout << "\n UNABLE TO OPEN GA OUTPUT FILE"
                << "\nFILENAME = " << myWorld->getOutputFileName() << std::endl;
      return;
    }
    if (!this->myWorld) return;

    // HEADER DATA
    this->prtHeader(this->outputFileStream);
}

// ============================================================
// Function : writeInput()
// ------------------------------------------------------------
// Write inputted data to output file
// ============================================================
void gaOutput::writeInput()
{
    if (!this->myWorld) return;

    // INPUT DATA
    this->prtWorld(this->outputFileStream);

    // INPUT OPTIONS
    this->prtOptions(this->outputFileStream);

    // FLUSH OUTPUT
    this->outputFileStream.flush();
}

// ============================================================
// Function : writeResults()
// ------------------------------------------------------------
// Write results to the output file
// ============================================================
void gaOutput::writeResults()
{
    if (outputFileStream.is_open()) {
      // WORLD DATA
      this->prtWorld(this->outputFileStream);

      // TIMINGS
      this->prtTail(this->outputFileStream);

      // CLOSE FILE
      this->outputFileStream.close();
    }
    else {
      std::cout << " Output file is not open, unable to write results to output file" << std::endl;
    }
}

// ============================================================
// Function : writeConvergence()
// ------------------------------------------------------------
// Write convergence data
// ============================================================
void gaOutput::writeConvergence()
{
    std::string filename = this->myWorld->getConvergFileName();

    this->convergFileStream.open(filename.c_str());

    if (!convergFileStream) {
      std::cout << "\n UNABLE TO OPEN GA CONVERGENCE FILE"
                << "\nFILENAME = " << myWorld->getConvergFileName() << std::endl;
      return;
    }
    if (!this->myWorld) return;

    std::vector<gaRegion*> regions = myWorld->getRegions();
    for (regionIterator c = regions.begin(); c != regions.end(); c++) {
      gaRegion* reg = *c;
      this->convergFileStream << "#Region:" << reg->getName() << std::endl;
      std::vector<double> converg = reg->getFitness();
      if (converg.size() > 0) {
        this->convergFileStream << "#Gen Fitness" << std::endl;
        for (unsigned int i = 0; i < converg.size(); i++) {
          this->convergFileStream << i+1 << " "  << converg[i] << std::endl;
        }
      }
    }
    this->convergFileStream.close();
}


// ============================================================
// Function : setProgramName()
// ------------------------------------------------------------
//
// ============================================================
void gaOutput::setProgramName(std::string progName)
{
    if (progName != "MTK++::GA") {
      progName = "MTK++::GA:" + progName;
    }
    this->programName = progName;
}

// ============================================================
// Function : prtHeader()
// ------------------------------------------------------------
// Write header to output file
// ============================================================
void gaOutput::prtHeader(std::ostream& os)
{
    char temp[81];
    std::string outline = "..............................................................................";

    copyright(os);
    os << "|                                                                              |"
       << std::endl;

    // program name
    std::string tempString = "Program:";
    sprintf(temp,"|%-78s|", (tempString.c_str()));
    os << temp << std::endl;

    int s = outline.size()-programName.size();
    std::string sub = outline.substr(0,s);
    sub = sub + programName;
    sprintf(temp,"|%78s|", sub.c_str());
    os << temp << std::endl;
    os << "|                                                                              |" << std::endl;

    // AUTHORS
    tempString = "Authors:";
    sprintf(temp,"|%-78s|", (tempString.c_str()));
    os << temp << std::endl;

    std::vector<std::string> authors = this->myWorld->getAuthors();
    for (unsigned int a = 0; a < authors.size(); a++) {
      s = outline.size() - authors[a].size();
      sub = outline.substr(0,s);
      sub = sub + authors[a];
      sprintf(temp,"|%78s|", sub.c_str());
      os << temp << std::endl;
    }
    os << "|                                                                              |" << std::endl;

    // VERSION
    tempString = "Version:";
    sprintf(temp,"|%-78s|", (tempString.c_str()));
    os << temp << std::endl;

    //s = outline.size() - strlen(VERSION);
    std::string strVERSION = VERSION;
    s = outline.size() - strVERSION.length();
    sub = outline.substr(0,s);
    sub = sub + strVERSION;
    sprintf(temp,"|%78s|", sub.c_str());
    os << temp << std::endl;
    os << "|                                                                              |" << std::endl;

    // DATE/TIME
    tempString = "Date:";
    sprintf(temp,"|%-78s|", (tempString.c_str()));
    os << temp << std::endl;

    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    //s = outline.size() - strlen(asctime (timeinfo));
    std::string ti(asctime (timeinfo));
    s = outline.size() - ti.length();
    sub = outline.substr(0,s);
    sub = sub + std::string(asctime (timeinfo));
    sprintf(temp,"|%78s|", sub.substr(0, sub.size()-1).c_str());
    os << temp << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "+------------------------------------------------------------------------------+\n" << std::endl;
}

// ============================================================
// Function : prtWorld()
// ------------------------------------------------------------
// Writes gaWorld to stream os
// ============================================================
void gaOutput::prtWorld(std::ostream& os)
{
    os << "+-----------------------------------GA WORLD-----------------------------------+" << std::endl;

    os << "World: " << myWorld->getName() << std::endl;

    std::vector<gaRegion*> regions = myWorld->getRegions();
    for (regionIterator c = regions.begin(); c != regions.end(); c++) {
      gaRegion* reg = *c;
      os << " Region: " << reg->getName() << std::endl;

      std::vector<gaPopulation*> pops = reg->getPopulations();
      for (populationIterator d = pops.begin(); d != pops.end(); d++) {
        gaPopulation* pop = *d;
        os << "  Population: " << pop->getName() << " id = " << pop->getId() << std::endl;

        std::vector<gaIndividual*> indList = pop->getIndividuals();
        for (individualIterator c = indList.begin(); c != indList.end(); c++) {
          gaIndividual* ind = *c;
          os << "   Individual: " << ind->getName() << " id = " << ind->getId()
             << " fitness = " << ind->getFitness() << std::endl;

          for (int l = 0; l < ind->getNumChromosomes(); l++) {
            gaChromosome* chr = ind->getChromosome(l,0,1);
            os << "    Chromosome: " << chr->getName() << std::endl;

            for (int t = 0; t < chr->getNumGenes(); t++) {
              gaGene* gene = chr->getGene(t,0,1);
              os << "     Gene: " << gene->getName() << std::endl;
              os << "      " << (*gene);
            }
          }
        }
      }
    }
    os << "+------------------------------------------------------------------------------+\n" << std::endl;
}

// ============================================================
// Function : prtOptions()
// ------------------------------------------------------------
// 
// ============================================================
void gaOutput::prtOptions(std::ostream& os)
{
    typedef std::vector<gaRegion*>::iterator regionIterator;
//    char temp[81];

    os << "+--------------------------------INPUT OPTIONS---------------------------------+" << std::endl;

    // WORLD OPTIONS
    os << "World: " << myWorld->getName() << std::endl;
    os << " Output file name = " << myWorld->getOutputFileName().c_str() << std::endl;
    os << " Restart file name = " << myWorld->getRestartFileName().c_str() << std::endl;
    os << " Convergence file name = " << myWorld->getConvergFileName().c_str() << std::endl;
    os << " Output verbosity = " << myWorld->getLevelOfOutput() << std::endl;
    os << " Chromomsomes per individual = " << myWorld->getChrPerInd() << std::endl;

    for (int i = 0; i < myWorld->chrPerInd; i++) {
      if (myWorld->genePerChr[i] == 1) {
        os << "  Chromosome #" << i+1 << " has " << myWorld->genePerChr[i] << " Gene " << std::endl;
      }
      else {
        os << "  Chromosome #" << i+1 << " has " << myWorld->genePerChr[i] << " Genes " << std::endl;
      }
    }

    // REGION OPTIONS
    std::vector<gaRegion*> regions = myWorld->getRegions();
    for (regionIterator c = regions.begin(); c != regions.end(); c++) {
      gaRegion* reg = *c;
      os << " Region: " << reg->getName() << std::endl;
      os << "  Maximum Number of Individuals = " << reg->maxInds << std::endl;
      os << "  Random Number Generator Seed = " << reg->seed << std::endl;
      os << "  Maximum Number of Generations = " << reg->maxGens << std::endl;
      os << "  Current Population = " << reg->curPop << std::endl;
      os << "  Number of Generations to Save = " << reg->popKeep << std::endl;
      os << "  Keep Percentage = " << reg->pKeep << std::endl;
      os << "  Crossover Percentage = " << reg->pCrossover <<  std::endl;
      os << "  Mutate Percentage = " << reg->pMutate << std::endl;
      os << "  Average Percentage = " << reg->pAverage << std::endl;
      os << "  Number of Children per pair = " << reg->nChild << std::endl;
      os << "  Chromosome/Gene Difference value = " << reg->chreDiff << std::endl;
      os << "  Selection Method used = " << reg->selection << std::endl;
      if (reg->elitism) {
        os << "  Elitism is on" << std::endl;
      }
      else {
        os << "  Elitism is off" << std::endl;
      }
      os << "  Selection Pressure = " << reg->selectionPressure << std::endl;
      os << "  Crossover Method used = " << reg->crossover << std::endl;
      os << "  Mutate Method used = " << reg->mutate << std::endl;
      os << "  Average Method used = " << reg->average << std::endl;
      os << "  Maximum Parameters = ";
      for (unsigned int i = 0; i < reg->maxParameters.size(); i++) {
        os << reg->maxParameters[i] << " ";
      }
      os << " " << std::endl;

      os << "  Minimum Parameters = ";
      for (unsigned int i = 0; i < reg->minParameters.size(); i++) {
        os << reg->minParameters[i] << " ";
      }
      os << " " << std::endl;

      os << "  Step Sizes = ";
      for (unsigned int i = 0; i < reg->stepSize.size(); i++) {
        os << reg->stepSize[i] << " ";
      }
      os << " " << std::endl;
    }
    os << "+------------------------------------------------------------------------------+\n" << std::endl;
}

// ============================================================
// Function : prtTail()
// ------------------------------------------------------------
// Prints end of output file
// ============================================================
void gaOutput::prtTail(std::ostream& os)
{
    char temp[81];
    std::cout << " gaOutput::prtTail " << os << std::endl;

    os << "+-----------------------------------TIMINGS------------------------------------+" << std::endl;
    os << "|                                                                              |" << std::endl;
    time (&this->myWorld->endTime);
    int diffTime = (int) difftime(this->myWorld->endTime, this->myWorld->startTime);

    if ( diffTime < 60 ) {
      sprintf(temp,"|%16s%-62d|", " # Seconds    = ", diffTime);
      os << temp << std::endl;
    }
    else if ( diffTime < 3600 ) {
      int min = (int) diffTime / 60;
      int sec = (int) diffTime % 60;
      sprintf(temp,"|%16s%-62d|", " # Minutes    = ", min);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Seconds    = ", sec);
      os << temp << std::endl;
    }
    else if ( diffTime < 86400 ) { // DAY
      int hours = (int) diffTime / 3600;
      int hourRemainder = (int) diffTime % 3600;
      int min = (int) hourRemainder / 60;
      int sec = (int) diffTime % 60;
      sprintf(temp,"|%16s%-62d|", " # Hours      = ", hours);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Minutes    = ", min);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Seconds    = ", sec);
      os << temp << std::endl;
    }
    else if ( diffTime < 31536000 ) { // YEAR
      int days = (int) diffTime / 86400;
      int daysRemainder = (int) diffTime % 86400;
      int hours = (int) daysRemainder / 3600;
      int hourRemainder = (int) (diffTime - 86400) % 3600;
      int min = (int) hourRemainder / 60;
      int sec = (int) diffTime % 60;
      sprintf(temp,"|%16s%-62d|", " # Days       = ", days);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Hours      = ", hours);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Minutes    = ", min);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Seconds    = ", sec);
      os << temp << std::endl;
    }
    else {
      int years = (int) diffTime / 31536000;
      int yearsRemainder = (int) diffTime % 31536000;
      int days = (int) yearsRemainder / 86400;
      int daysRemainder = (int) diffTime % 86400;
      int hours = (int) daysRemainder / 3600;
      int hourRemainder = (int) (diffTime - 86400) % 3600;
      int min = (int) hourRemainder / 60;
      int sec = (int) diffTime % 60;
      sprintf(temp,"|%16s%-62d|", " # Years      = ", years);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Days       = ", days);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Hours      = ", hours);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Minutes    = ", min);
      os << temp << std::endl;
      sprintf(temp,"|%16s%-62d|", " # Seconds    = ", sec);
      os << temp << std::endl;
    }
    os << "|                                                                              |" << std::endl;
    os << "+------------------------------------------------------------------------------+\n" << std::endl;
}

// ============================================================
// Function : prtWarning()
// ------------------------------------------------------------
// Prints warning message
// ============================================================
void gaOutput::prtWarning(std::ostream& os, std::string warning)
{
    myWorld->setStatus(1);
    char temp[81];
    if (os) {
      os << "+-----------------------------------WARNING------------------------------------+" << std::endl;
      os << "|                                                                              |" << std::endl;
      sprintf(temp,"| %-76s |", warning.c_str());
      os << temp << std::endl;
      os << "|                                                                              |" << std::endl;
      os << "+------------------------------------------------------------------------------+\n" << std::endl;
    }
}

// ============================================================
// Function : writeWarning()
// ------------------------------------------------------------
// Write warning message to output file
// ============================================================
void gaOutput::writeWarning(std::string warning)
{
    this->prtWarning(this->outputFileStream, warning);
}

// ============================================================
// Function : prtError()
// ------------------------------------------------------------
// Prints error message
// ============================================================
void gaOutput::prtError(std::ostream& os, std::string error)
{
    char temp[81];
    if (os) {
      os << "+------------------------------------ERROR-------------------------------------+" << std::endl;
      os << "|                                                                              |" << std::endl;
      sprintf(temp,"| %-76s |", error.c_str());
      os << temp << std::endl;
      os << "|                                                                              |" << std::endl;
      os << "+------------------------------------------------------------------------------+\n" << std::endl;
    }
    myWorld->setStatus(2);
}

// ============================================================
// Function : writeError()
// ------------------------------------------------------------
// Write error message to output file
// ============================================================
void gaOutput::writeError(std::string error)
{
    this->prtError(this->outputFileStream, error);

    // FLUSH OUTPUT
    this->outputFileStream.flush();
}

} // MTKpp namespace

