/*!
   \file atomTypeParser.h
   \brief atom type parser
   \author Martin Peters

   Parses atom type xml files using xercesc

   $Date: 2010/03/29 20:39:34 $
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

#ifndef ATOMTYPEPARSER_H
#define ATOMTYPEPARSER_H

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#ifdef USE_XERCES
#include "domErrorHandler.h"
#endif // USE_XERCES

namespace MTKpp
{

class atomTypes;
struct atomTypeTMP;

// ============================================================
// Class : atomTypeParser()
// ------------------------------------------------------------
/*! 
   \class atomTypeParser
   \brief Parses atom type xml files using xercesc
   \author Martin Peters
*/
// ============================================================
class atomTypeParser
{
public:

    /*!
       \brief atomTypeParser Constructor
    */
    atomTypeParser(atomTypes*);

    //! atomTypeParser Destructor
    virtual ~atomTypeParser();
 
    /*!
       \brief Read atom type xml file
       \param i Input file
       \return boolean
    */
    int Read(std::string i);

#ifdef USE_XERCES
protected:

    //! Main reading function
    void atomTypesFiller(DOMNode*);

    //! Individual atom type reading function
    void atomTypeFiller(DOMNode*);

    //! parameter reading function
    void parameterFiller(DOMNode*);
#endif // USE_XERCES

protected:

    //! atomTypes pointer
    atomTypes*          pAtomTypes;

    //! atomType pointer
    atomTypeTMP*           pAtomType;
};

} // MTKpp namespace

#endif // ATOMTYPEPARSER_H

