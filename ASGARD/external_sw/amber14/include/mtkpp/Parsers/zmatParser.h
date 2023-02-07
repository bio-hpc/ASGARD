/*!
   \file zmatParser.h
   \brief Parses Z-Matrix files
   \author Martin Peters

   Reads and writes guassian files

   $Date: 2007/09/14 11:17:26 $
   $Revision: 1.3 $

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

#ifndef ZMATPARSER_H
#define ZMATPARSER_H

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#include "baseParser.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
class element;
class connections;
struct Bond;
struct Angle;
struct Torsion;
struct Improper;
class vector3d;

// ============================================================
// Class : zmatParser()
// ------------------------------------------------------------
/*! 
   \class zmatParser
   \brief Reads and writes z-matrix format files
   \author Martin Peters
   \version 0.1
   \date 2006
*/
// ============================================================
class zmatParser : public baseParser

{
public:

    /*!
       \brief zmatParser Constructor
    */
    zmatParser();

    //! zmatParser Destructor
    ~zmatParser();

    /*!
       \brief Generates Z-matrix
       \param pMolecule molecule pointer
       \param zmatrix see below
       \param zmatData see below
       \param atomMap A map between the old atom ordering and the new

       zmatrix:
        -# [[symbol, bondIndex, bondName, angleIndex, angleName, torsionIndex, torsionName]]

       zmatrix data:
        -# [bondName:size, angleName:size, torsionName:size]
    */
    int            genZmatrix(molecule* pMolecule,
                              std::vector<std::vector<std::string> > &zmatrix,
                              std::map<std::string, double> &zmatData,
                              std::map<int, int> &atomMap);

    /*!
       \brief Read zmatrix formatted file
       \param i Input file
       \param c collection pointer
       \param zmatrix see below
       \param zmatData see below

       zmatrix:
        -# [[symbol, bondIndex, bondName, angleIndex, angleName, torsionIndex, torsionName]]

       zmatrix data:
        -# [bondName:size, angleName:size, torsionName:size]
    */
    void           Read(const std::string &i, collection* c,
                        std::vector<std::vector<std::string> > &zmatrix,
                        std::map<std::string, double> &zmatData);

    /*!
       \brief Write zmatrix formatted input file
       \param o Output file
       \param m molecule pointer
    */
    void           Write(const std::string &o, molecule* m);

    /*!
       \brief Write zmatrix formatted input file
       \param o Output file
       \param m molecule pointer
       \param coordinates Vector of coordinates which overrides the molecules coordinates
    */
    void           Write(const std::string &o, molecule* m, std::vector< vector3d > &coordinates);

    /*!
       \brief Write zmatrix formatted input file
       \param o Output file
       \param c collection pointer
       \param molID The ID of the molecule in the collection
    */
    void           Write(const std::string &o, collection* c, const int &molID);

};

} // MTKpp namespace

#endif // ZMATPARSER_H
