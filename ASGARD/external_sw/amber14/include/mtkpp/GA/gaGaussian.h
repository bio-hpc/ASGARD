/*! 
   \file gaGaussian.h
   \brief Contains gaussian information used by gaSelection
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

#ifndef GAGAUSSIAN_H
#define GAGAUSSIAN_H

#include <iostream>
#include <string>
#include <vector>
#include <string>

#include "Utils/constants.h"

namespace MTKpp
{

// ============================================================
// Class : gaGaussian()
// ------------------------------------------------------------
/*! 
   \class gaGaussian
   \brief Contains gaussian information used by gaSelection
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaGaussian
{
public:

    /*!
       \brief gaGaussian Constructor
    */
    gaGaussian();

    //! gaGaussian Destructor
    virtual ~gaGaussian();

    /*!
       \brief Set up gaGaussian
       \param st Start value for the gaussian
       \param en End value for the gaussian
       \param stDe Standard deviation of the gaussian (called selection pressure in input file)
       \param numInd Number of gaChromosomes per gaIndividual
    */
    void setup(const int &st, const int &en, const double &stDe, const int &numInd);

    /*!
       \brief Set up gaGaussian
       \param randomNum Picks gaIndividual using the gaussian and the random number
    */
    int pick(double randomNum);

protected:

    //! Generates the gaussian using the parameter provided to setup()
    void makeGaussian();

    //! Values that make up the gaussian
    std::vector<double>           gaussian;

    //! Start value
    int                           start;

    //! End value
    int                           end;

    //! Number of gaIndividuals in the population
    int                           numIndividuals;

    //! Standard deviation of the gaussian
    double                        stdDev;

    //! 1/(stdDev*2.506628)
    double                        invStdDevPi2;
};

} // MTKpp namespace

#endif // GAGAUSSIAN_H
