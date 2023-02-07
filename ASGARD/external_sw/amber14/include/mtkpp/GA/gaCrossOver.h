/*! 
   \file gaCrossOver.h
   \brief Performs recombination of gaChromosomes
   \author Martin Peters

   CrossOver Functions Available:

   1) random single gene crossover per chromosome

   2) random multiple gene crossover per chromosome

   3) random loci crossover over all chromosomes

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

#ifndef GACROSSOVER_H
#define GACROSSOVER_H

#include <iostream>
#include <string>
#include <vector>

#include "Utils/constants.h"

namespace MTKpp
{

class gaOutput;
class gaPopulation;
class gaIndividual;

// ============================================================
// Class : gaCrossOver()
// ------------------------------------------------------------
/*! 
   \class gaCrossOver
   \brief Performs recombination of gaChromosomes
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaCrossOver
{
public:

    /*!
       \brief gaCrossOver Constructor
       \param out gaOutput pointer
       \param nChild number of gaIndividuals to be produced by averaging
       \param crossover Type of crossover procedure employed
       \param chrPerInd Number of gaChromosomes per gaIndividual
       \param genePerChr Number of gaGenes per gaChromosome
    */
    gaCrossOver(gaOutput *out, std::vector<int> genePerChr, int nChild = 1, 
                std::string crossover = "random-locus", int chrPerInd = 1);

    //! gaCrossOver Destructor
    virtual ~gaCrossOver();

    /*!
       \brief Main function of gaCrossOver
       \param ind1 First gaIndividual
       \param ind2 Second gaIndividual
       \param pop Current gaPopulation
       \param name New gaIndividual name
       \return successful operation
    */
    bool reproduce(gaIndividual* ind1, gaIndividual* ind2, gaPopulation* pop, std::string name);

protected:

    /*!
       \brief Single gene crossover per chromosome
       \param newInd First gaIndividual
       \param parent Second gaIndividual

       \code
         {} == Individual
         || == chromosome
         [] == gene

         ind1: {|[aaaa][cccc]||[ffff][tttt]|}
                                             --> {|[aaaa][dddd]||[hhhh][tttt]|}
         ind2: {|[bbbb][dddd]||[hhhh][jjjj]|}
       \endcode

    */
    bool __singleGene(gaIndividual* newInd, gaIndividual* parent);

    /*!
       \brief Multiple gene crossover per chromosome
       \param newInd First gaIndividual
       \param parent Second gaIndividual

       \code
         {} == Individual
         || == chromosome
         [] == gene

        ind1: {|[aa][bb][cc]||[dd][ee][ff]|}
                                            --> {|[aa][bb][xx]||[gg][ee][ii]|}
        ind2: {|[zz][yy][xx]||[gg][hh][ii]|}
       \endcode
    */
    bool __multipleGene(gaIndividual* newInd, gaIndividual* parent);

    /*!
       \brief Random locus crossover per chromosome
       \param newInd First gaIndividual
       \param parent Second gaIndividual

       \code
         {} == Individual
         || == chromosome
         [] == gene

                    |                  |
        ind1: {|[aa]|[bb][cc]||[dd][ee]|[ff]|}
                    |                  |      --> {|[aa][yy][xx]||[dd][ee][ii]}
        ind2: {|[zz]|[yy][xx]||[gg][hh]|[ii]|}
                    |                  |
        \endcode
    */
    bool __randomLocus(gaIndividual* newInd, gaIndividual* parent);

protected:

    //! gaOperator pointer
    gaOutput*                     pGaOutput;

    //! Number of offspring produced by each pair of parents
    int                           nChild;

    //! Type of crossover to be carried out
    std::string                   crossover;

    //! Number of gaChromsomes per gaIndividual
    int                           chrPerInd;

    //! Number of gaGenes per gaChromsome
    std::vector<int>              genePerChr;
};

} // MTKpp namespace

#endif // GACROSSOVER_H
