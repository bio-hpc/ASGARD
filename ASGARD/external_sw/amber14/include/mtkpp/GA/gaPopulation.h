/*! 
   \file gaPopulation.h
   \brief Container for gaIndividuals
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

#ifndef GAPOPULATION_H
#define GAPOPULATION_H

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <sstream>

#include "Utils/constants.h"

namespace MTKpp
{

class gaRegion;
class gaIndividual;

// ============================================================
// Class : gaPopulation()
// ------------------------------------------------------------
/*! 
   \class gaPopulation
   \brief Container for gaIndividuals
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaPopulation
{
public:

    /*!
       \brief gaPopulation Constructor
       \param parent gaRegion pointer
    */
    gaPopulation(gaRegion *parent = 0);

    //! gaRegion Destructor
    virtual ~gaPopulation();

    /*!
       \brief Add gaIndividual to gaPopulation
       \return gaIndividual pointer
    */
    gaIndividual* addIndividual();

    /*!
       \brief Add gaIndividual to gaPopulation
       \param ind gaIndividual pointer
       \return gaIndividual pointer
    */
    gaIndividual* addIndividual(gaIndividual* ind);

    /*!
       \brief Delete gaIndividual from gaPopulation
       \param ind gaIndividual to be deleted
       \return gaIndividual pointer
    */
    void delIndividual(gaIndividual* ind);

    /*!
       \brief Get gaIndividual from gaPopulation
       \param number gaIndividual number, either id or index
       \param id bool return by id
       \param index bool return by index
       \return gaIndividual pointer
    */
    gaIndividual* getIndividual(int number, bool id, bool index);

    /*!
       \brief Get gaIndividual from gaPopulation
       \param name gaIndividual name
       \return gaIndividual pointer
    */
    gaIndividual* getIndividual(std::string name);

    /*!
       \brief Get gaIndividual list from gaPopulation
       \return list of gaIndividual
    */
    std::vector<gaIndividual*> getIndividuals();

    /*!
       \brief Get last gaIndividual to be added to gaPopulation
       \return gaIndividual pointer
    */
    gaIndividual* getLastAdded();

    /*!
       \brief Setup gaPopulation
       \param maxInds Maximum number of gaIndividuals
    */
    void setup(const int& maxInds);

    /*!
       \brief Initialize gaPopulation
    */
    void initialize();

    /*!
       \brief Get absolute values for all gaIndividuals in gaPopulation
    */
    void getAbsValues();

    /*!
       \brief Rank all gaIndividuals in gaPopulation
    */
    void rank();

    //-------------//
    // - get/set - //
    //-------------//

    /*!
       \brief Set id for gaPopulation
       \param id gaPopulation id
    */
    void setId(int id);

    /*!
       \brief Get id of gaPopulation
       \return id of gaPopulation
    */
    int getId();

    /*!
       \brief Set name for gaPopulation
       \param name gaPopulation name
    */
    void setName(std::string name);

    /*!
       \brief Get name of gaPopulation
       \return name of gaPopulation
    */
    std::string getName();

    /*!
       \brief Get gaRegion which gaPopulation is a member of
       \return gaRegion pointer
    */
    gaRegion* getParent();

    /*!
       \brief Get number of gaIndividuals in gaPopulation
       \return number of gaIndividuals in gaPopulation
    */
    int getNumIndividuals();

    /*!
       \brief Get best gaIndividual in gaPopulation in terms of fitness
       \param i index of best gaIndividual
       \return gaIndividual pointer
    */
    gaIndividual* getBestIndividual(int &i);

protected:

    //! gaIndividual pointer
    typedef std::vector<gaIndividual*>::iterator individualIterator;

    //! gaIndividual list
    std::vector<gaIndividual*>    itsIndividualList;

    //! gaRegion pointer
    gaRegion*                     pParent;

    //! gaIndividual pointer
    gaIndividual*                 pGaIndividual;

    //! gaIndividual pointer
    gaIndividual*                 bestIndividual;

    //! gaPopulation id
    int                           itsId;

    //! gaPopulation name
    std::string                   itsName;

};

} // MTKpp namespace

#endif // GAPOPULATION_H
