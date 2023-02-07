/*! 
   \file gaOutput.h
   \brief Class to handle the GA output
   \author Martin Peters

   $Date: 2010/03/29 20:24:52 $
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

#ifndef GAOUTPUT_H
#define GAOUTPUT_H

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

#include <stdio.h>
#include <string.h>
#include <time.h>

namespace MTKpp
{

class gaWorld;
class gaRegion;
class gaPopulation;
class gaIndividual;
class gaChromosome;
class gaGene;

// ============================================================
// Class : gaOutput()
// ------------------------------------------------------------
/*! 
   \class gaOutput
   \brief Class to handle the GA output
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaOutput
{
public:

    /*!
       \brief gaOutput Constructor
       \param w gaWorld pointer
    */
    gaOutput(gaWorld* w);

    //! gaOutput destructor
    virtual ~gaOutput();

    /*!
       \brief Open Output File
    */
    void openOutputFile();

    /*!
       \brief Write GA output file
    */
    void writeInput();

    /*!
       \brief Write GA output file
    */
    void writeResults();

    /*!
       \brief Write GA convergence file
    */
    void writeConvergence();


    /*!
       \brief Set program name
       \param progName Program name
    */
    void setProgramName(std::string progName = "MTK++::GA");

    /*!
       \brief Write GA header
    */
    void prtHeader(std::ostream& os);

    /*!
       \brief Print gaWorld information
       \param os Ouput stream
    */
    void prtWorld(std::ostream& os);

    /*!
       \brief Print inputed options
       \param os Ouput stream
    */
    void prtOptions(std::ostream& os);

    /*!
       \brief Print inputed options
       \param os Ouput stream
    */
    void prtTail(std::ostream& os);

    /*!
       \brief Print warning message
       \param os Ouput stream
       \param warning Warning message
    */
    void prtWarning(std::ostream& os, std::string warning);

    /*!
       \brief Write warning message to output file
       \param warning Warning message
    */
    void writeWarning(std::string warning);

    /*!
       \brief Print error message
       \param os Ouput stream
       \param error Error message
    */
    void prtError(std::ostream& os, std::string error);

    /*!
       \brief Write error message to output file
       \param error Error message
    */
    void writeError(std::string error);

protected:

    //! gaWorld pointer
    gaWorld*                      myWorld;

    //! Program name
    std::string                   programName;

    //! Output File Stream
    std::ofstream                 outputFileStream;

    //! Convergence File Stream
    std::ofstream                 convergFileStream;

    //! region iterator
    typedef std::vector<gaRegion*>::iterator regionIterator;

    //! population iterator
    typedef std::vector<gaPopulation*>::iterator populationIterator;

    //! individual iterator
    typedef std::vector<gaIndividual*>::iterator individualIterator;
};

} // MTKpp namespace

#endif // GAOUTPUT_H
