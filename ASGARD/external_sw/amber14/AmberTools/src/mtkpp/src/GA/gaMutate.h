/*! 
   \file gaMutate.h
   \brief Performs mutation of gaChromosomes
   \author Martin Peters

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

#ifndef GAMUTATE_H
#define GAMUTATE_H

#include <iostream>
#include <string>
#include <vector>

#include "Utils/constants.h"

namespace MTKpp
{

class gaIndividual;
class gaOperators;

// ============================================================
// Class : gaMutate()
// ------------------------------------------------------------
/*! 
   \class gaMutate
   \brief Performs mutation of gaChromosomes
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaMutate
{
public:

    /*!
       \brief gaMutate Constructor
       \param parent gaOperators pointer
       \param mutate Type of mutation procedure employed
       \param chrPerInd Number of gaChromosomes per gaIndividual
       \param genePerChr Number of gaGenes per gaChromosome
    */
    gaMutate(gaOperators *parent, std::vector<int> genePerChr, 
             std::string mutate = "single-gene", int chrPerInd = 1);

    //! gaMutate Destructor
    virtual ~gaMutate();

    /*!
       \brief Main function of gaMutate
       \param ind gaIndividual pointer
       \return successful operation
    */
    bool mutation(gaIndividual* ind);

protected: // functions

    /*!
       \brief Single gene mutation
       \param ind gaIndividual pointer
    */
    bool __singleGene(gaIndividual* ind);

    /*!
       \brief Multiple gene mutation
       \param ind gaIndividual pointer
    */
    bool __multipleGene(gaIndividual* ind);

protected: // data

    //! gaOperator pointer
    gaOperators*                  pParent;

    //! Type of mutation to be carried out
    std::string                   mutate;

    //! Number of gaChromsomes per gaIndividual
    int                           chrPerInd;

    //! Number of gaGenes per gaChromsome
    std::vector<int>              genePerChr;
};

} // MTKpp namespace

#endif // GAMUTATE_H
