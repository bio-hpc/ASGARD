/*! 
   \file gaSelection.h
   \brief Selects gaIndividuals from gaPopulation based on fitness
   \author Martin Peters

   $Date: 2010/03/29 20:24:52 $
   $Revision: 1.5 $

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

#ifndef GASELECTION_H
#define GASELECTION_H

#include <iostream>
#include <string>
#include <vector>

namespace MTKpp
{

class gaOperators;
class gaGaussian;

// ============================================================
// Class : gaSelection()
// ------------------------------------------------------------
/*! 
   \class gaSelection
   \brief Selects gaIndividuals from gaPopulation based on fitness
 
    Selection Functions Available:

    1) Random

    2) Semi-random

    3) Proportional fitness / Roulette wheel

    4) Tournament

    5) Truncation

   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaSelection
{
public:

    /*!
       \brief gaSelection Constructor
       \param parent gaOperators pointer
       \param selection Type of selection procedure employed
       \param maxInds Maximum gaIndividuals in the gaPopulation
    */
    gaSelection(gaOperators *parent = 0, std::string selection = "semi-random", int maxInds = 0);

    //! gaSelection Destructor
    virtual ~gaSelection();

    //! 
    int select();

    /*!
       \brief Create a gaGaussian
       \param start Start value for the gaussian
       \param end End value for the gaussian
       \param selPress The standard deviation of the gaussian
    */
    void setGaussian(int start, int end, double selPress);

protected: // Functions

    //! Returns uniform random individuals
    int __random();

    //! Returns semi-random individuals weighted by a normal gaussian distribution
    int __semi_random();

    /*! 
       \brief returns individuals using proportional fitness/roulette wheel
       \todo implement
    */
    int __rouletteWheel();

    /*! 
       \brief returns individuals by perform a tournament (size N) of fitness
       \todo implement
    */
    int __tournament();

    /*! 
       \brief Used by the keep operator a certain percentage of the 
              generation is allowed to live without competition or selection
       \todo implement
    */
    int __truncation();

protected: // Data

    //! gaOperators pointer
    gaOperators*                  pParent;

    //! gaGaussian pointer
    gaGaussian*                   pGaGaussian;

    //! Selection method used
    std::string                   selection;

    //! Maximum number of gaIndividuals in the gaPopulation
    int                           maxInds;
};

} // MTKpp namespace

#endif // GASELECTION_H
