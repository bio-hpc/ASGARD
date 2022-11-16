/*! 
   \file gaGene.h
   \brief Class to handle genes
   \author Martin Peters

   Class to handle genes

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

#ifndef GAGENE_H
#define GAGENE_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>

#include "Utils/constants.h"

namespace MTKpp
{

class gaRegion;
class gaChromosome;

// ============================================================
// Class : gaGene()
// ------------------------------------------------------------
/*! 
   \class gaGene
   \brief Class to handle genes
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaGene
{
public:

    /*!
       \brief gaGene Constructor
       \param parent gaChromosome pointer
    */
    gaGene(gaChromosome *parent = 0);

    /*!
       \brief gaGene Copy Constructor
    */
    gaGene(gaGene*);

    //! gaIndividual Destructor
    virtual ~gaGene();

    /*!
       \brief Set gene equal to rhs
       \param rhs gaGene pointer
    */
    void setGene(gaGene* rhs);

    /*!
       \brief Add bit to gene
       \param b bit to be added
    */
    void addBit(double b);

    /*!
       \brief Get bit from gene
       \param i index of bit in gene
       \return bit
    */
    double getBit(const int& i);

    /*!
       \brief Initialize gaGene
    */
    void initialize();

    /*!
       \brief Get the absolute value of all bits summed together
       \return absolute vales
    */
    double getAbsValue();


    /*!
       \brief Compare two gaGenes
       \param rhs Second gaGene
       \return true/false
    */
    bool compare(gaGene* rhs);


    /*!
       \brief Mutate gaGene
    */
    void mutate();


    /*!
       \brief Average gaGene, updates this gaGene
       \param rhs Second gaGene
    */
    void average(gaGene* rhs);

    /*!
       \brief Formatted print of the gaGene
    */
    void printToScreen();

    friend std::ostream& operator<< (std::ostream& os, const gaGene& g) {
      for (unsigned int i = 0; i < g.itsBits.size(); i++) {
        os << g.itsBits[i] << " ";
      }
      os << "" << std::endl;
      return os;
    }


    //-------------//
    // - GET/SET - //
    //-------------//

    /*!
       \brief Set id of gaGene
       \param id gaGene id
    */
    void setId(int id);

    /*!
       \brief Get id of gaGene
       \return id of gaGene
    */
    int getId();

    /*!
       \brief Set name of gaGene
       \param name gaGene name
    */
    void setName(std::string name);

    /*!
       \brief Get name of gaGene
       \return name of gaGene
    */
    std::string getName();

    /*!
       \brief Get gaChromosome which gaGene is a member of
       \return gaChromosome pointer
    */
    gaChromosome* getParent();

    /*!
       \brief Get number of bits in gaGene
       \return number of bits in gaGene
    */
    int getNumBits();

    /*!
       \brief Get genetic information of the gaGene
       \return genetic info
    */
    std::vector<double> getGeneticInformation();

protected:

    //! gaRegion pointer
    gaRegion*                     curRegion;

    //! gaChromosome pointer
    gaChromosome*                 pParent;

    //! gaGene bits
    std::vector<double>           itsBits;

    //! gaGene id
    int                           itsId;

    //! gaGene name
    std::string                   itsName;

};

} // MTKpp namespace

#endif // GAGENE_H
