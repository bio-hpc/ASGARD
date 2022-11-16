/*! 
   \file gaChromosome.h
   \brief Container for gaGenes
   \author Martin Peters

   Container for gaGenes

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

#ifndef GACHROMOSOME_H
#define GACHROMOSOME_H

#include <iostream>
#include <string>
#include <vector>

#include "Utils/constants.h"

namespace MTKpp
{

class gaIndividual;
class gaGene;

// ============================================================
// Class : gaChromosome()
// ------------------------------------------------------------
/*! 
   \class gaChromosome
   \brief Container for gaGenes
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaChromosome
{
public:

    /*!
       \brief gaChromosome Constructor
       \param parent gaIndividual pointer
    */
    gaChromosome(gaIndividual *parent = 0);

    /*!
       \brief gaChromosome Copy Constructor
       \param rhs gaIndividual pointer
    */
    gaChromosome(gaChromosome* rhs);

    //! gaIndividual Destructor
    virtual ~gaChromosome();

    /*!
       \brief Add gaGene to the gaChromosome
       \return gaGene pointer
    */
    gaGene* addGene();

    /*!
       \brief Add gaGene to the gaChromosome
       \param rhs gaGene pointer
       \return gaGene pointer
    */
    gaGene* addGene(gaGene* rhs);

    /*!
       \brief Delete gaGene from gaChromosome
       \param gen gaGene pointer
    */
    void delGene(gaGene* gen);

    /*!
       \brief Set gaGene in gaChromosome
       \param n number of gaGene either id or index
       \param id bool return by id
       \param index bool return by index
       \param rhs gaChromosome pointer
    */
    void setGene(int n, bool id, bool index, gaChromosome* rhs);

    /*!
       \brief Get gaGene from gaChromosome
       \param n number of gaGene either id or index
       \param id bool return by id
       \param index bool return by index
       \return gaGene pointer
    */
    gaGene* getGene(int n, bool id, bool index);

    /*!
       \brief Mutate gaGene
       \param n number of gaGene either id or index
       \param id bool return by id
       \param index bool return by index
    */
    void mutateGene(int n, bool id, bool index);

    /*!
       \brief Average gaGene
       \param n number of gaGene either id or index
       \param id bool return by id
       \param index bool return by index
       \param rhs gaChromosome pointer
    */
    void averageGene(int n, bool id, bool index, gaChromosome* rhs);

    /*!
       \brief Initialize gaChromosome
    */
    void initialize();

    /*!
       \brief Get the absolute value of all the gaGenes summed up
       \return Absolute value
    */
    double getAbsValue();

    /*!
       \brief Compate two gaChromosomes
       \param rhs other gaChromosome
       \return true/false
    */
    bool compare(gaChromosome* rhs);

    //-------------//
    // - get/set - //
    //-------------//

    /*!
       \brief Set id of gaChromosome
       \param id gaChromosome id
    */
    void setId(int id);

    /*!
       \brief Get id of gaChromosome
       \return gaChromosome id
    */
    int getId();

    /*!
       \brief Set name of gaChromosome
       \param name gaChromosome name
    */
    void setName(std::string name);

    /*!
       \brief Get name of gaChromosome
       \return gaChromosome name
    */
    std::string getName();

    /*!
       \brief Get gaIndividual which gaChromosome is a member of
       \return gaIndividual pointer
    */
    gaIndividual* getParent();

    /*!
       \brief Get number of genes in gaChromosome
       \return number of genes in gaChromosome
    */
    int getNumGenes();

    /*!
       \brief Get genetic information of the gaChromosome
       \return genetic info
    */
    std::vector<double> getGeneticInformation();

protected:

    //! gaGene iterator
    typedef std::vector<gaGene*>::iterator geneIterator;

    //! gaGene vector
    std::vector<gaGene*>          itsGeneList;

    //! gaIndividual pointer
    gaIndividual*                 pParent;

    //! gaGene pointer
    gaGene*                       pGaGene;

    //! gaChromosome id
    int                           itsId;

    //! gaChromosome name
    std::string                   itsName;

};

} // MTKpp namespace

#endif // GACHROMOSOME_H
