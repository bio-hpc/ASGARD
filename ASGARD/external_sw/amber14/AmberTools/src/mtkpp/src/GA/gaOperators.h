/*! 
   \file gaOperators.h
   \brief This class performs the GA operations including crossover
          mutation, averaging, etc.
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

#ifndef GAOPERATORS_H
#define GAOPERATORS_H

#include <iostream>
#include <string>
#include <vector>

#include "Utils/constants.h"

namespace MTKpp
{

class gaOutput;
class gaRegion;
class gaPopulation;
class gaSelection;
class gaCrossOver;
class gaMutate;
class gaAverage;

// ============================================================
// Class : gaOperators()
// ------------------------------------------------------------
/*! 
   \class gaOperators
   \brief This class performs the GA operations including crossover, 
          mutation, averaging, etc.
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaOperators
{
public:

    /*!
       \brief gaOperators Constructor
       \param parent gaRegion pointer
    */
    gaOperators(gaRegion *parent = 0);

    //! gaOperators Destructor
    virtual ~gaOperators();

    /*!
       \brief Set up all GA operators including selection, crossover, mutate, averaging

       Calculates nKeep, nCrossover, nMutate, and nAverage
    */
    void setup();

    /*!
       \brief nKeep individuals of the new generation are individuals from the previous generation
       \param curPop Current gaPopulation
       \param prePop Previous gaPopulation
    */
    void keep(gaPopulation* curPop, gaPopulation* prePop);

    /*!
       \brief nCrossover individuals of the new generation are produced by recombination

              Each time two parents are selected only one child is produced

              Two parents are allowed to produce nChild offspring
       \param curPop Current gaPopulation
       \param prePop Previous gaPopulation
    */
    void crossover(gaPopulation* curPop, gaPopulation* prePop);

    /*!
       \brief nMutate individuals of the new generation are produced by mutation
       \param curPop Current gaPopulation
       \param prePop Previous gaPopulation
    */
    void mutate(gaPopulation* curPop, gaPopulation* prePop);

    /*!
       \brief nAverage individuals of the new generation are produced by averaging
       \param curPop Current gaPopulation
       \param prePop Previous gaPopulation
    */
    void average(gaPopulation* curPop, gaPopulation* prePop);

    /*!
       \brief Removes clones from the population
       \param curPop Current gaPopulation
       \param prePop Previous gaPopulation
    */
    void removeRedundant(gaPopulation* curPop, gaPopulation* prePop);

protected:

    //! gaRegion pointer
    gaRegion*    pParent;

    //! gaSelection pointer
    gaSelection* pGaSelection;

    //! gaCrossOver pointer
    gaCrossOver* pGaCrossOver;

    //! gaMutate pointer
    gaMutate*    pGaMutate;

    //! gaAverage pointer
    gaAverage*   pGaAverage;

    //! gaOutput pointer
    gaOutput*    pGaOutput;

    //! nKeep individuals of the new generation are individuals from the previous generation
    int nKeep;

    //! nCrossover individuals of the new generation are produced by recombination
    int nCrossover;

    //! nMutate individuals of the new generation are produced by mutation
    int nMutate;

    //! nAverage individuals of the new generation are produced by averaging
    int nAverage;

    //! chromosome must difference greater than this quantity to be saved when removing clones
    double chreDiff;

    //! mutant name
    std::string mutantName;
};

} // MTKpp namespace

#endif // GAOPERATORS_H
