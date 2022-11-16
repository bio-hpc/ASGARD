/*! 
   \file gaAverage.h
   \brief Performs chromosome averaging
   \author Martin Peters

   $Date: 2010/03/29 20:24:52 $
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

#ifndef GAAVERAGE_H
#define GAAVERAGE_H

#include <iostream>
#include <string>
#include <vector>

namespace MTKpp
{

class gaOperators;
class gaIndividual;
class gaPopulation;

// ============================================================
// Class : gaAverage()
// ------------------------------------------------------------
/*! 
   \class gaAverage
   \brief Performs chromosome averaging
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaAverage
{
public:

    /*!
       \brief gaAverage Constructor
       \param parent gaOperators pointer
       \param nChild number of gaIndividuals to be produced by averaging
       \param aver Type of averaging procedure employed
       \param chrPerInd Number of gaChromosomes per gaIndividual
       \param genePerChr Number of gaGenes per gaChromosome
    */
    gaAverage(gaOperators *parent,std::vector<int> genePerChr, int nChild = 1, 
              std::string aver = "single-gene", int chrPerInd = 1);

    //! gaAverage destructor
    virtual ~gaAverage();

    /*!
       \brief Main function of gaAverage
       \param ind1 First gaIndividual
       \param ind2 Second gaIndividual
       \param pop Current gaPopulation
       \param name New gaIndividual name
       \return successful operation
    */
    bool average(gaIndividual* ind1, gaIndividual* ind2, gaPopulation* pop, std::string name);

protected:

    /*!
       \brief Single gene averaging
       \param newInd First gaIndividual
       \param parent Second gaIndividual
    */
    bool __singleGene(gaIndividual* newInd, gaIndividual* parent);

    /*!
       \brief Multiple gene averaging
       \param newInd First gaIndividual
       \param parent Second gaIndividual
    */
    bool __multipleGene(gaIndividual* newInd, gaIndividual* parent);

private:

    //! gaOperator pointer
    gaOperators*                  pParent;

    //! Number of offspring produced by each pair of parents
    int                           nChild;

    //! Type of Averaging to be carried out
    std::string                   aver;

    //! Number of gaChromsomes per gaIndividual
    int                           chrPerInd;

    //! Number of gaGenes per gaChromsome
    std::vector<int>              genePerChr;
};

} // MTKpp namespace

#endif // GAAVERAGE_H
