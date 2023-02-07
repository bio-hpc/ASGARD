/*! 
   \file gaWorld.h
   \brief Container for gaRegions
   \author Martin Peters

   This is the main class of the GA.
   Container for gaRegions

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

#ifndef GAWORLD_H
#define GAWORLD_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <time.h>

#include "Utils/constants.h"

namespace MTKpp
{

class gaRegion;
class gaOutput;

// ============================================================
// Class : gaWorld()
// ------------------------------------------------------------
/*! 
   \class gaWorld
   \brief Container for gaRegions
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaWorld
{
    friend class gaRegion;
    friend class gaOutput;
public:

    /*!
       \brief gaWorld Constructor
       \param name gaWorld
    */
    gaWorld(std::string name = "world");

    //! gaWorld Destructor
    virtual ~gaWorld();

    /*!
       \brief Add gaRegion to the gaWorld
       \return gaRegion pointer
    */
    gaRegion* addRegion();

    /*!
       \brief Delete gaRegion in the gaWorld
       \param reg gaRegion to be deleted
    */
    void delRegion(gaRegion* reg);

    /*!
       \brief Get gaRegion in the gaWorld
       \param n number of gaRegion, either id or index
       \param id bool return by id
       \param index bool return by index
       \return gaRegion pointer
    */
    gaRegion* getRegion(int n,bool id, bool index);

    /*!
       \brief Get gaRegion in the gaWorld
       \param name name of gaRegion
       \return gaRegion pointer
    */
    gaRegion* getRegion(std::string name);

    /*!
       \brief Get gaRegions in the gaWorld
       \return vector of gaRegion pointers
    */
    std::vector<gaRegion*> getRegions();

    /*!
       \brief Get number of gaRegions in the gaWorld
       \return number of gaRegions
    */
    int getNumRegions();

    /*!
       \brief Set up gaWorld
       \param chrPerInd number of gaChromosome per gaIndividual
       \param genePerChr number of gaGenes per gaChromosomes
       \param geneSizes size of each gaGene
       \param outputFile Output File Name
       \param restartFile Restart File Name
       \param convergFile Convergence File Name
       \param outputLevel Output verbosity control, levels include:
         -# 1  (default)
         -# 2
         -# 3
    */
    void setup(const int& chrPerInd, std::vector<int> genePerChr, std::vector<int> geneSizes,
               const std::string outputFile, const std::string restartFile, const std::string convergFile, const int outputLevel);

    /*!
       \brief Get gaOutput pointer
       \return gaOutput pointer
    */
    gaOutput* getGaOutput();

    /*!
       \brief Initialize gaWorld
    */
    void initialize();

    /*!
       \brief Rank gaRegions
    */
    void rank();

    //---------------//
    // - GET / SET - //
    //---------------//

    /*!
       \brief Set the gaWorld name
       \param name gaWorld name
    */
    void setName(std::string name);

    /*!
       \brief Get the gaWorld name
       \return gaWorld name
    */
    std::string getName();

    /*!
       \brief Get the number of gaChromosomes per gaIndividual
       \return number of gaChromosomes per gaIndividual
    */
    int getChrPerInd();

    /*!
       \brief Set the number of gaGenes per gaChromosome
       \param gPc number of gaGenes per gaChromosome
    */
    void setGenePerChr(std::vector<int> gPc);

    /*!
       \brief Get the number of gaGenes per gaChromosome
       \return number of gaGenes per gaChromosome
    */
    std::vector<int> getGenePerChr();

    /*!
       \brief Set the gaGenes sizes
       \param gS gaGene sizes
    */
    void setGeneSizes(std::vector<int> gS);

    /*!
       \brief Get the gaGenes sizes
       \return gaGene sizes
    */
    std::vector<int> getGeneSizes();

    /*!
       \brief Get the output file name
       \return output file name
    */
    std::string getOutputFileName();

    /*!
       \brief Get the restart file name
       \return restart file name
    */
    std::string getRestartFileName();

    /*!
       \brief Get the convergence file name
       \return convergence file name
    */
    std::string getConvergFileName();

    /*!
       \brief Get the level of output
       \return level of output
    */
    int getLevelOfOutput();

    /*!
       \brief Set author
       \param name Name and email address of author

       \code
         example: "john doe (johndoe@ga.com)"
       \endcode
    */
    void setAuthor(std::string name);

    /*!
       \brief Get authors
       \return vector of author names
    */
    std::vector<std::string> getAuthors();

    /*!
       \brief Set status
       \param status
         -# 0 ok
         -# 1 warnings
         -# 2 error
    */
    void setStatus(int status);

    /*!
       \brief Get status
       \return status
         -# 0 ok
         -# 1 warnings
         -# 2 error
    */
    int getStatus();

protected:

    //! gaRegion iterator
    typedef std::vector<gaRegion*>::iterator regionIterator;

    //! gaRegion vector
    std::vector<gaRegion*>   itsRegionList;

    //! gaRegion pointer
    gaRegion*                pGaRegion;

    //! gaOutput pointer
    gaOutput*                pGaOutput;

    //! gaWorld name
    std::string              itsName;

    //! number of gaChromosomes per gaIndividual
    int                      chrPerInd;

    //! number of gaGenes per gaChromosome
    std::vector<int>         genePerChr;

    //! gaGene sizes
    std::vector<int>         geneSizes;

    //! Output file name
    std::string              outputFile;

    //! Restart file name
    std::string              restartFile;

    //! Convergence file name
    std::string              convergFile;

    //! Level of output to be produced
    int                      levelOfOutput;

    //! Energy function Authors
    std::vector<std::string> authors;

    /*!
       \brief Status of the program
         -# 0 ok
         -# 1 warnings
         -# 2 error
    */
    int                      status;

    //! Start time
    time_t                   startTime;

    //! End time
    time_t                   endTime;
};

} // MTKpp namespace

#endif // GAWORLD_H
