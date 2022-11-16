/*!
   \file amberParser.h
   \brief Parses amber prmtop/crd files
   \author Martin Peters

   Reads and writes amber prmtop/crd files

   $Date: 2007/09/14 11:17:26 $
   $Revision: 1.2 $

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

#ifndef AMBERPARSER_H
#define AMBERPARSER_H

#include <map>
#include <list>

#include "StringManip.h"
#include "baseParser.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
class vector3d;
struct Bond;

// ============================================================
// Class : amberParser()
// ------------------------------------------------------------
/*!
   \class amberParser
   \brief Reads and writes amber prmtop/crd files
   \author Martin Peters
   \version 0.1
   \date 2007
*/
// ============================================================
class amberParser : public baseParser
{

public:
    /*!
       \brief amberParser Constructor
    */
    amberParser();

    //! amberParser Destructor
    ~amberParser();
    
    /*!
       \brief Write amber prmtop file
       \param inpcrd output inpcrd file
       \param prmtop output prmtop file
       \param c collection pointer
    */
    void           Write(const std::string &inpcrd, const std::string &prmtop, collection* c);

};

} // MTKpp namespace

#endif // AMBERPRMTOPPARSER_H
