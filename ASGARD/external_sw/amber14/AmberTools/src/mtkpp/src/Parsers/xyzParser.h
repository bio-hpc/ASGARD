/*!
   \file xyzParser.h
   \brief Parses XYZ files
   \author Martin Peters

   $Date: 2010/02/20 03:19:36 $
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

#ifndef XYZPARSER_H
#define XYZPARSER_H

#include <map>
#include <iostream>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <vector>

#include "baseParser.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
class element;
struct Bond;

// ============================================================
// Class : xyzParser()
// ------------------------------------------------------------
/*! 
   \class xyzParser
   \brief Reads and writes XYZ format files
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class xyzParser : public baseParser
{
public:

    /*!
       \brief xyzParser Constructor
    */
    xyzParser();

    //! xyzParser Destructor
    ~xyzParser();

    /*!
       \brief Read xyz file
       \param i input file
       \param c collection pointer
    */
    void           Read(const std::string &i,collection* c);

    /*!
       \brief Write xyz file
       \param o output file
       \param m Molecule pointer
    */
    void           Write(const std::string &o, molecule* m);

    //! collection pointer
    collection*    pCol;

    //! molecule pointer
    molecule*      pMol;

    //! submolecule pointer
    submolecule*   pSmol;

    //! atom pointer
    atom*          pAtom;
};

} // MTKpp namespace

#endif // XYZPARSER_H

