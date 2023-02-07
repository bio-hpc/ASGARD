/*!
   \file prepParser.h
   \brief Parses AMBER prep files
   \author Martin Peters

   $Date: 2008/09/24 13:57:16 $
   $Revision: 1.5 $

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

#ifndef PREPPARSER_H
#define PREPPARSER_H

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

class collection;
class molecule;
class submolecule;
class atom;
class element;
struct Bond;

class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;
struct stdBond;
struct stdImproper;
struct stdLoop;
struct stdAlias;
struct stdRing;
struct stdFeature;

// ============================================================
// Class : prepParser()
// ------------------------------------------------------------
/*! 
   \class prepParser
   \brief Reads and writes AMBER prep format files
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class prepParser : public baseParser
{
public:

    /*!
       \brief prepParser Constructor
    */
    prepParser();

    //! prepParser Destructor
    ~prepParser();

    /*!
       \brief open prep file
       \param i input file
    */
    int            openFile(const std::string &i);

    /*!
       \brief read fragment header section
       \param name fragment name
       \param symbol fragment symbol
       \param charge fragment charge
    */
    int           readHeader(std::string &name, std::string &symbol, double &charge);

    /*!
       \brief read fragment main section
    */
    void           readFragment(stdFrag* pStdFrag);

    /*!
       \brief Read prep file
       \param i input file
       \param g stdGroup pointer
    */
    void           Read(const std::string &i, stdGroup* g);

    /*!
       \brief Read prep file
       \param i input file
       \param f stdFrag pointer
    */
    void           Read(const std::string &i, stdFrag* f);

    /*!
       \brief Write prep file
       \param o output file
       \param s stdFrag pointer
    */
    void           Write(const std::string &o, stdFrag* s);

    /*!
       \brief Write prep file
       \param o output file
       \param g stdGroup pointer
    */
    void           Write(const std::string &o, stdGroup* g);

    //! collection pointer
    collection*    pCol;

    //! molecule pointer
    molecule*      pMol;

    //! submolecule pointer
    submolecule*   pSmol;

    //! atom pointer
    atom*          pAtom;

    //! file handler
    std::ifstream iprep;
};

} // MTKpp namespace

#endif // PREPPARSER_H

