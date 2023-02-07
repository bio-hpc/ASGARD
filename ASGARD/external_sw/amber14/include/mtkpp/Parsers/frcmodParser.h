/*!
   \file frcmodParser.h
   \brief Parses AMBER frcmod files
   \author Martin Peters

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.4 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2007  (see AUTHORS file for a list of contributors)

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

#ifndef FRCMODPARSER_H
#define FRCMODPARSER_H

#include <map>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <vector>

#include "baseParser.h"

namespace MTKpp
{

class parameters;
struct atomType;
struct bondParam;
struct angleParam;
struct torsionParam;
struct improperParam;
struct hBondParam;
struct equivalentAtomsParam;

// ============================================================
// Class : frcmodParser()
// ------------------------------------------------------------
/*! 
   \class frcmodParser
   \brief Reads and writes AMBER frcmod format files
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class frcmodParser : public baseParser
{
public:

    /*!
       \brief frcmodParser Constructor
    */
    frcmodParser(parameters* c, std::string groupName);

    //! acParser Destructor
    ~frcmodParser();

    /*!
       \brief Read an AMBER frcmod file
       \param i input file name
    */
    void Read(const std::string &i);

    /*!
       \brief Write an AMBER frcmod file
       \param o output file
       \param p parameter set name
    */
    void Write(const std::string &o, const std::string &p);

protected:

    //! group name
    std::string              groupName;

    //! parameters pointer
    parameters*              pParameters;

    //! atomType pointer
    atomType*                pAtomType;

    //! bondParam pointer
    bondParam*               pBondParam;

    //! angleParam pointer
    angleParam*              pAngleParam;

    //! torsionParam pointer
    torsionParam*            pTorsionParam;

    //! improperParam pointer
    improperParam*           pImproperParam;

    //! hBondParam pointer
    hBondParam*              pHBondParam;

    //! equivalentAtomsParam pointer
    equivalentAtomsParam*    pEquivalentAtomsParam;

};

} // MTKpp namespace

#endif // FRCMODPARSER_H

