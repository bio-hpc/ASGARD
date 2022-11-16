/*! 
   \file gaIndividual.h
   \brief Container for gaChromosomes
   \author Martin Peters

   Container for gaChromosomes

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

#ifndef GAINDIVIDUAL_H
#define GAINDIVIDUAL_H

#include <iostream>
#include <string>
#include <vector>

#include "Utils/constants.h"

namespace MTKpp
{

class gaPopulation;
class gaChromosome;

// ============================================================
// Class : gaIndividual()
// ------------------------------------------------------------
/*! 
   \class gaIndividual
   \brief Container for gaChromosomes
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaIndividual
{
public:

    /*!
       \brief gaIndividual Constructor
       \param parent gaPopulation pointer
    */
    gaIndividual(gaPopulation *parent = 0);

    /*!
       \brief gaIndividual Copy Constructor
       \param rhs gaIndividual pointer
       \param parent gaPopulation pointer
    */
    gaIndividual(gaIndividual* rhs, gaPopulation *parent);

    //! gaIndividual Destructor
    virtual ~gaIndividual();

    /*!
       \brief Add gaChromosome to the gaIndividual
       \return gaChromosome pointer
    */
    gaChromosome* addChromosome();

    /*!
       \brief Add chromosome of the individual
       \param chr gaChromosome pointer
       \return gaChromosome pointer
    */
    gaChromosome* addChromosome(gaChromosome* chr);

    /*!
       \brief Delete chromosome from individual
       \param chr gaChromosome pointer
    */
    void delChromosome(gaChromosome* chr);

    /*!
       \brief Get chromosome from individual
       \param n number of chromosome either id or index
       \param id bool return by id
       \param index bool return by index
       \return gaChromosome pointer
    */
    gaChromosome* getChromosome(int n, bool id, bool index);

    /*!
       \brief Initialize gaIndividual
    */
    void initialize();

    /*!
       \brief Get the absolute value of all the gaGenes summed up
       \return Absolute value
    */
    double getAbsValue();

    /*!
       \brief Compare two gaIndividual with respect to its chromosomal make up
       \param rhs gaIndividual pointer
       \return true/false
    */
    bool compare(gaIndividual* rhs);

    //bool operator<(const gaIndividual* rhs);
    //bool operator==(const gaIndividual* rhs);

    /*!
       \brief Compares two gaIndividuals based on fitness
       \param lhs first gaIndividual
       \param rhs Second gaIndividual
       \return boolean

       After this function is defined the STL sort() function in gaPopulation can be used
    */
    static bool less(const gaIndividual *lhs, const gaIndividual *rhs)
    {
        return lhs->itsFitness < rhs->itsFitness;
    }

    /*!
       \brief Compares two gaIndividuals based on fitness
       \param lhs first gaIndividual
       \param rhs Second gaIndividual
       \return boolean

       After this function is defined the STL sort() function in gaPopulation can be used
    */
    static bool greater(const gaIndividual *lhs, const gaIndividual *rhs)
    {
        return lhs->itsFitness > rhs->itsFitness;
    }

    /*!
       \brief Formatted print of the gaIndividual
    */
    void printToScreen();

    //---------------//
    // - GET / SET - //
    //---------------//

    /*!
       \brief Set the individual id
       \param id individual int
    */
    void setId(int id);

    /*!
       \brief Get gaIndividual id
       \return gaIndividual id
    */
    int getId();

    /*!
       \brief Set gaIndividual name
       \param name gaIndividual name
    */
    void setName(std::string name);

    /*!
       \brief Get gaIndividual name
       \return gaIndividual name
    */
    std::string getName();

    /*!
       \brief Set the gaIndividual fitness value
       \param f fitness value
    */
    void setFitness(double f);

    /*!
       \brief Get fitness value
       \return gaIndividual fitness value
    */
    double getFitness();

    /*!
       \brief Set the gaIndividual evaluate value
       \param b evaluate value
    */
    void setEvaluate(bool b);

    /*!
       \brief Get evaluate value
       \return gaIndividual evaluate value
    */
    bool getEvaluate();

    /*!
       \brief Get gaPopulation which gaIndividual is a member of
       \return gaPopulation pointer
    */
    gaPopulation* getParent();

    /*!
       \brief Get number of gaChromosomes in gaIndividual
       \return number of gaChromosomes
    */
    int getNumChromosomes();

    /*!
       \brief Get single vector of all genetic information of the gaIndividual
       \return genetic info
    */
    std::vector<double> getGeneticInformation();

    /*!
       \brief Set gaGene
       \param chr gaChromosome index in gaIndividual
       \param gen gaGene index in gaChromosome
       \param ind gaIndividual that copied
    */
    void setGene(int chr, int gen, gaIndividual* ind);

    /*!
       \brief Mutate gaGene
       \param chr gaChromosome index in gaIndividual
       \param gen gaGene index in gaChromosome
    */
    void mutateGene(int chr, int gen);

    /*!
       \brief Average gaGene
       \param chr gaChromosome index in gaIndividual
       \param gen gaGene index in gaChromosome
       \param ind gaIndividual that copied
    */
    void averageGene(int chr, int gen, gaIndividual* ind);

    /*!
       \brief Add parent gaIndividual
       \param p parent name
    */
    void addParent(std::string p);

    /*!
       \brief Get parents
       \return parents
    */
    std::vector<std::string> getParents();

    /*!
       \brief List of partner which gaIndividual has mated with
       \param partnerName partner name
    */
    void addPartner(std::string partnerName);

    /*!
       \brief Get partners
       \return partners
    */
    std::vector<std::string> getPartners();

    /*!
       \brief Returns number of partners
       \return number of partners
    */
    int numPartners();

    /*!
       \brief Has it mated with a certain gaIndividual
       \param partnerName partner name
       \return number of times it has mated with partnerName
    */
    int mated(std::string partnerName);

protected:

    //! gaChromosome iterator
    typedef std::vector<gaChromosome*>::iterator chromosomeIterator;

    //! gaChromosome vector
    std::vector<gaChromosome*>    itsChromosomeList;

    //! gaPopulation pointer
    gaPopulation*                 pParent;

    //! gaChromosome pointer
    gaChromosome*                 pGaChromosome;

    //! gaIndividual id in gaPopulation
    int                           itsId;

    //! gaIndividual name
    std::string                   itsName;

    //! gaIndividual genetic absolute value
    double                        abs_value;

    //! gaIndividual fitness vales
    double                        itsFitness;

    //! evaluate gaIndividual
    bool                          bEvaluate;

    //! partners vector
    std::vector<std::string>      itsPartners;

    //! partners vector
    std::vector<std::string>      itsParents;

};

} // MTKpp namespace

#endif // GAINDIVIDUAL_H
