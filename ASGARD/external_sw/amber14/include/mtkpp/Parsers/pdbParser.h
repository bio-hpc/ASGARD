/*!
   \file pdbParser.h
   \brief Parses pdb files
   \author Martin Peters
   \author Duane Williams

   Reads and writes pdb files

   $Date: 2010/05/04 20:23:42 $
   $Revision: 1.18 $

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

#ifndef PDBPARSER_H
#define PDBPARSER_H

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
class metalCenter;
class vector3d;
struct Bond;

/*!
   \struct pdbInfo
   \brief Container for PDB info
   \author Martin Peters
   \date 2007
*/
struct pdbInfo {
    //! resolution
    double resolution;

    //! experimental technique
    std::string expTechnique;
    
    //! Remarks
    std::vector<std::string> remarks;
};


// ============================================================
// Class : pdbParser()
// ------------------------------------------------------------
/*!
   \class pdbParser
   \brief Reads and writes PDB format files
   \author Martin Peters 
   \author Duane Williams
   \version 0.1
   \date 2005
*/
// ============================================================
class pdbParser : public baseParser
{

public:
    /*!
       \brief pdbParser Constructor
    */
    pdbParser();

    //! pdbParser Destructor
    ~pdbParser();

    /*!
       \brief Read PDB file
       \param i input file
       \param c collection pointer
    */
    void           Read(const std::string &i, collection* c);

    /*!
       \brief Read PDB file
       \param i input file
       \param c collection pointer
       \param m molecule pointer
    */
    void           Read(const std::string &i, collection* c, molecule* m);

    /*!
       \brief Update coordinates of molecule from PDB file
       \param i input file
       \param m molecule pointer
    */
    void           updateMolCoords(const std::string &i, molecule* m);

    /*!
       \brief Write PDB file
       \param o output file
       \param m molecule pointer
    */
    void           Write(const std::string &o, molecule* m);

    /*!
       \brief Write PDB file
       \param o output file
       \param m vector of molecules
    */
    void           Write(const std::string &o, std::vector<molecule*> m);

    /*!
       \brief Write PDB file
       \param o output file
       \param c collection pointer
    */
    void           Write(const std::string &o, collection* c);

    /*!
       \brief Write PDB file
       \param o output file
       \param c collection pointer
       \param d d?
    */
    void           Write(const std::string &o, collection* c, const double d);

    /*!
       \brief Write PDB formatted file
       \param o Output file
       \param m molecule pointer
       \param coordinates Vector of coordinates which overrides the molecules coordinates
    */
    void            Write(const std::string &o, molecule* m, std::vector< vector3d > &coordinates);

    /*!
       \brief Write PDB formatted file
       \param o Output file
       \param m metalCenter pointer
    */
    void            Write(const std::string &o, metalCenter* m);

    /*!
       \brief Write PDB file
       \param os output stream
       \param c collection pointer
    */
    void            Write(std::ostream& os, collection* c);

    /*!
       \brief Get Total charge
    */
    double getTotalCharge() {
      return totalCharge;
    };

    /*!
       \brief Has total charge
    */
    bool hasTotalCharge() {
      return hasTotalChargeRemark;
    };

    /*!
       \brief 
    */
    pdbInfo*       itsPdbInfo;

    /*!
       \brief Total charge on the chemical specie described by the pdb
    */
    double totalCharge;
    bool hasTotalChargeRemark;

protected: // functions

    /*!
       \brief Create pdbInfo object
    */
    void           setInfo();

    /*!
       \brief Get the 1-Letter residue code
       \param s 3-Letter code
        \return 1-Letter code
    */
    std::string get1LCode(std::string s);

    /*!
       \brief Get the type of the current atom
       
       - Pro -- All protein atoms
       - Met -- Zn, Ca atoms
       - Wat -- Water atom
       - Lig -- Ligand atoms
       \param s atom name
       \param r residue name
       \param a atom bool
       \param b hetatm bool
       \return atom type
    */
    std::string determineType(std::string s, std::string r, bool a, bool b);

protected: // data

    //////////////////////////////////////////
    // - pointers to the backbone objects - //
    //////////////////////////////////////////

    //! collection pointer
    collection*    pCollection;

    //! molecule pointer
    molecule*      pMolecule;

    //! submolecule pointer
    submolecule*   pSubMolecule;

    //! atom pointer
    atom*          pAtom;

    //! atom pointer
    atom*          pAtom1;

    //! atom pointer
    atom*          pAtom2;

    //! Bond pointer
    Bond*          pBond;

    //! internal atom mapping local to pdbParser
    std::map<int, int>  itsAtomMap;

    //! coordinate pointer
    vector3d*      pCoord1;

    //! coordinate pointer
    vector3d*      pCoord2;

    //! molecule iterator
    typedef std::vector<molecule*>::iterator moleculeIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator atomIterator;

    //! atom map used to generate CONECT info
    typedef std::multimap<atom*,atom*> itsConectMMap_ptrs;

    //! atom list used to generate CONECT info
    std::list<atom*>    listConect_ptrs;

    //! Residue name to 1 letter code
    std::map<std::string, std::string> res21l;

    //! residue name map iterator
    typedef std::map<std::string, std::string>::iterator nameMapIterator;
};

} // MTKpp namespace

#endif // PDBPARSER_H
