/*!
   \file molParser.h
   \brief Parses mol files
   \author Martin Peters

   Reads and writes mol files

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

#ifndef MOLPARSER_H
#define MOLPARSER_H

#include <map>
#include <vector>
#include <iostream>
#include <string.h>
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
// Class : molParser()
// ------------------------------------------------------------
/*! 
   \class molParser
   \brief Reads and writes MOL format files 
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class molParser : public baseParser

{
public:

    /*!
       \brief molParser Constructor
    */
    molParser();

    //! molParser Destructor
    ~molParser();

    /*!
       \brief Read MOL formatted file
       \param i Input file
       \param m molecule pointer
       \param bohr Convert coordinates to Bohr or not
    */
    void            Read(const std::string &i, molecule* m, const bool &bohr = false);

    /*!
       \brief Read MOL formatted file
       \param i Input file
       \param c collection pointer
       \param bohr Convert coordinates to Bohr or not
    */
    void            Read(const std::string &i, collection* c, const bool &bohr = false);

    /*!
       \brief Write MOL formatted file
       \param o Output file
       \param m molecule pointer
    */
    void            Write(const std::string &o, molecule* m);

    /*!
       \brief Write MOL formatted file
       \param o Output file
       \param m molecule pointer
       \param coordinates Vector of coordinates which overrides the molecules coordinates
    */
    void            Write(const std::string &o, molecule* m, std::vector< vector3d > &coordinates);

    /*!
       \brief Write MOL formatted file
       \param o Output file
       \param c collection pointer
       \param molID The ID of the molecule in the collection
    */
    void            Write(const std::string &o, collection* c, const int &molID);

protected: // functions

    //! Read a single line of the input file
    /*!
       \param fileline Actual line
       \param pthisAtom AtomLine struct to temporary store atom info
    */
    bool            ReadAtomLine(std::string &fileline, AtomLine* pthisAtom);

    //! Get atom
    /*!
       \param i atom index in itsAtomMap
    */
    int             getAtom(int i);

protected: //data

    //////////////////////////////////////////
    // - pointers to the backbone objects - //
    //////////////////////////////////////////

    //! molecule pointer
    molecule*       pMolecule;

    //! submolecule pointer
    submolecule* pSubMolecule;

    //! atom pointer
    atom*           pAtom;

    //! Bond pointer
    Bond*           pBond;

    //! atom pointer used in bonding section
    atom*           pBondAtom1;

    //! atom pointer used in bonding section
    atom*           pBondAtom2;

    //! internal atom mapping local to molParser
    std::map<int, int>   itsAtomMap;

    // connections pointer
    connections* pConnections;
};

} // MTKpp namespace

#endif // MOLPARSER_H
