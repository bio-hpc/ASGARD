/*!
   \file mol2Parser.h
   \brief Parsers mol2 files
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

#ifndef MOL2PARSER_H
#define MOL2PARSER_H

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

struct Bond;

// ============================================================
// Class : mol2Parser()
// ------------------------------------------------------------
/*! 
   \class mol2Parser
   \brief Reads and writes MOL2 format files
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class mol2Parser : public baseParser
{
public:

    /*!
       \brief baseParser Constructor
    */
    mol2Parser();

    //! baseParser Destructor
    ~mol2Parser();

    /*!
       \brief Read MOL2 file
       \param i input file
       \param c collection pointer
    */
    void            Read(const std::string &i, collection* c);

    /*!
       \brief Read MOL2 file
       \param i input file
       \param m molecule pointer
    */
    void            Read(const std::string &i, molecule* m);


    /*!
       \brief Write MOL2 file
       \param o input file
       \param m molecule pointer
    */
    void            Write(const std::string &o, molecule* m);

    //! get charge type
    std::string getChargeType() {
      return chargeType;
    };

    //! get total charge
    double getTotalCharge() {
      return totalCharge;
    };

    //! has total charge
    bool hasTotalCharge() {
      return hasTotalChargeRemark;
    };

    //! collection pointer
    collection*    pCol;

    //! molecule pointer
    molecule*      pMol;

    //! submolecule
    submolecule*   pSmol;

    //! atom pointer
    atom*          pAtom;

    //! Bond pointer
    Bond*          pBond;

    //! atom pointer, Used in bonding section
    atom*          pBondAtom1;

    //! atom pointer, Used in bonding section
    atom*          pBondAtom2;

    //! Charge type
    std::string chargeType;

    /*!
       \brief Total charge on the chemical specie described by the mol2
    */
    double totalCharge;

    //! has total charge remark
    bool hasTotalChargeRemark;

protected: // functions

    /*!
       \brief Read atom line
       \param fileline file line
       \param pthisAtom atomLine pointer
       \return boolean
    */
    bool            ReadAtomLine(std::string &fileline, AtomLine* pthisAtom);

    /*!
       \brief Get atom index from map
       \param i atom index
       \return atom index
    */
    int             getAtom(int i);

    /*!
       \brief Get the 1-Letter residue code
       \param s 3-Letter code
       \return 1-Letter code
    */
    std::string get1LCode(std::string s);

    /*!
       \brief Asists building the backbone after accumulating the atoms structures as a
        residue or submolecule set.  Adds a submolecule to the current molecule per
        vector of  AtomLine's.
       \param pCollection ?
       \param residueID ?
       \param n ?
       \param atomLineVector ?
    */
    void buildupSubmolecule(collection* pCollection, std::string& residueID, unsigned int n, std::vector<AtomLine*>& atomLineVector);

protected: // data

    //! map used in bonding section of write
    std::map<int, int>  itsAtomMap;

    //! Residue name to 1 letter code
    std::map<std::string, std::string> res21l;

    //! residue name map iterator
    typedef std::map<std::string, std::string>::iterator nameMapIterator;
};

} // MTKpp namespace

#endif // MOL2PARSER_H
