/*! 
   \file amberNonBonded.h
   \brief AMBER non-bonded energy and gradient
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
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

#ifndef AMBERNONBONDED_H
#define AMBERNONBONDED_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

#include "Utils/constants.h"

namespace MTKpp
{

class amber;

// ============================================================
// Class : amberNonBonded()
// ------------------------------------------------------------
/*! 
   \class amberNonBonded

   \brief AMBER non-bonded energy and gradient

   \author Martin Peters

   \version 0.1

   \date 2006

   \section AMBERNonBondedEnG AMBER Non-Bonded Energy And Gradient

*/
// ============================================================
class amberNonBonded
{
public:
    /*!
       \brief amberNonBonded Constructor
    */
    amberNonBonded();

    /*!
       \brief amberNonBonded Constructor
       \param pAmber parent potential
    */
    amberNonBonded(amber* pAmber);

    /*!
       \brief amberNonBonded Destructor
    */
    virtual ~amberNonBonded();

    /*!
       \brief Calculate amberNonBonded Energy
       \return amberNonBonded energy
    */
    double calculateE();

    /*!
       \brief Calculate amberNonBonded Energy and Gradient
       \return amberNonBonded energy
    */
    double calculateG();

    /*!
       \brief Decompose the amberNonBonded Energy
    */
    double decompose();

protected:
    amber*         pAmber;

    void setToOne(int v[], int s);

};

} // MTKpp namespace

#endif // AMBERNONBOND_H

