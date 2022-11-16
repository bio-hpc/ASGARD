/*!
   \file gaParser.h
   \brief Parses GA xml files using Qt & xercesc
   \author Martin Peters

   $Date: 2010/03/29 20:39:35 $
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

#ifndef GAPARSER_H
#define GAPARSER_H

#include <iostream>
#include <string>
#include <vector>
#include "baseParser.h"

#ifdef USE_XERCES
#include "domErrorHandler.h"
#endif // USE_XERCES

namespace MTKpp
{

class gaWorld;
class gaRegion;
class gaPopulation;
class gaIndividual;
class gaChromosome;
class gaGene;

// ============================================================
// Class : gaParser()
// ------------------------------------------------------------
/*! 
   \class gaParser
   \brief Reads GA xml files
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class gaParser : public baseParser
{
public:

    /*!
       \brief gaParser Constructor
    */
    gaParser(gaWorld*);

    //! gaParser Destructor
    virtual ~gaParser();

    /*!
       \brief Reads GA world and input xml files
       \param w world xml file
       \param i input xml file
       \return boolean
    */
    int Read(std::string w, std::string i);

    /*!
       \brief Reads GA world xml file
       \param i Input file
       \return boolean
    */
    int ReadWorld(std::string i);

    /*!
       \brief Reads GA input xml file
       \param i Input file
       \return boolean
    */
    int ReadInput(std::string i);

#ifdef USE_XERCES
protected:

    /*!
       \brief gaWorld reader
       \param d dom node
    */
    void gaWorldFiller(DOMNode* d);

    /*!
       \brief gaRegion reader
       \param d dom node
    */
    void gaRegionFiller(DOMNode* d);

    /*!
       \brief gaPopulation reader
       \param d dom node
    */
    void gaPopulationFiller(DOMNode* d);

    /*!
       \brief gaIndividual reader
       \param d dom node
    */
    void gaIndividualFiller(DOMNode* d);

    /*!
       \brief gaChromosome reader
       \param d dom node
    */
    void gaChromosomeFiller(DOMNode* d);

    /*!
       \brief gaGene reader
       \param d dom node
    */
    void gaGeneFiller(DOMNode* d);

    // gaInput functions

    /*!
       \brief gaWorld input reader
       \param d dom node
    */
    void gaWorldInputFiller(DOMNode* d);

    /*!
       \brief gaRegion input reader
       \param d dom node
    */
    void gaRegionInputFiller(DOMNode* d);

    /*!
       \brief gaPopulation input reader
       \param d dom node
    */
    void gaPopulationInputFiller(DOMNode* d);
#endif // USE_XERCES

protected:

    //! gaWorld pointer
    gaWorld*       myWorld;

    //! gaRegion pointer
    gaRegion*      myRegion;

    //! gaPopulation pointer
    gaPopulation*  myPopulation;

    //! gaIndividual pointer
    gaIndividual*  myIndividual;

    //! gaChromosome pointer
    gaChromosome*  myChromosome;

    //! gaGene pointer
    gaGene*        myGene;

    //! Was the world file read in?
    bool bWorldFileRead;
};

} // MTKpp namespace

#endif // GAPARSER_H
