/*!
   \file sdfParser.h
   \brief Parses sd files
   \author Martin Peters

   Reads and writes sd files

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.14 $

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

#ifndef SDFPARSER_H
#define SDFPARSER_H

#include <map>
#include <iostream>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <vector>

#include <math.h>
#include <stdexcept>

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

// ============================================================
// Class : sdfParser()
// ------------------------------------------------------------
/*! 
   \class sdfParser
   \brief Reads and writes SDF format files 
   \author Martin Peters
   \version 0.1
   \date 2005

   \section sdfExample Sdf File Example
   - File Contents
    \verbatim
1,2,2,5,5-Pentamethyl-3-imidazoline 3-oxide
  CDK    1/18/06,2:25
NMRShiftDB 2189 http://www.nmrshiftdb.org:8080/portal/_kbohn_Fri Nov 08 10:04:07
 11 11  0  0  0  0  0  0  0  0999 V2000
   -0.0897    1.5477   -0.1004 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0827    0.2946    0.0756 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2603   -0.3941    0.1731 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2226    0.5932   -0.3763 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5189    1.8781   -0.1370 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9639    2.4834    1.1958 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8006    2.8552   -1.2801 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9570    2.4917   -0.2359 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2676   -1.6724   -0.6676 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5939   -0.7107    1.6323 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3909    0.5845    0.5141 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
  2  3  1  0  0  0  0
  1  8  1  0  0  0  0
  3  4  1  0  0  0  0
  3  9  1  0  0  0  0
  4  5  1  0  0  0  0
  3 10  1  0  0  0  0
  5  1  1  0  0  0  0
  4 11  1  0  0  0  0
M  CHG  1   1   1
M  CHG  1   8  -1
M  END
> <Spectrum 13C 0>
132.1;0.0D;1|23.3;0.0Q;8|23.3;0.0Q;9|23.5;0.0Q;5|23.5;0.0Q;6|26.1;0.0Q;10|60.5;0.0S;2|90.0;0.0S;4|

> <Field Strength [MHz]>
0:50.328

> <Solvent>
0:Acetone-D6 ((CD3)2CO)

> <Temperature [K]>
0:298

$$$$
    \endverbatim

   - First line is the compound name

   - Second Line contains initials, program, date, dimension, scale factors, energy, reg name

   - Third line contains comments

   - Fourth line provides Number of atoms(N), Number of bonds(B), atlist, obs, chiral, nentries, obs, obs, obs, obs, nprop

   - Fifth line onwards provides the x, y, z, atom symbol, massdiff, charge, stereo, hcount, stcare, valence, for each atom

   - (N+5)th line starts the bonding info (at1 at2 type stereo topology)

   - Charges are set after the bonding info. e.g M CHG  1  1  1

   - Info after "M   END" is user defined and in this case come from the nmrshiftdb

   - Now take a look at it graphically
    \image html mdlFormat.png
    \image latex mdlFormat.eps
*/
// ============================================================

class sdfParser : public baseParser
{
public:

    /*!
       \brief sdfParser Constructor
    */
    sdfParser();

    //! sdfParser Destructor
    ~sdfParser();

    /*!
       \brief Read SDF formatted files
       \param i Input file
       \param c collection pointer
       \param bohr Convert coordinates to Bohr or not
    */
    void           Read(const std::string &i, collection* c, const bool &bohr = false);

    /*!
       \brief Read a molecule from a sd files
       \param os input stream
       \param m molecule pointer
       \param bohr Convert coordinates to Bohr or not
    */
    int            ReadMolecule(std::ifstream& os, molecule* m, const bool &bohr = false);

    /*!
       \brief Determine the number of molecules in an sd file
       \param i Input file
    */
    int            numMolecules(const std::string &i);

    /*!
       \brief Write SDF formatted file
       \param o Output file
       \param c collection pointer
    */
    void           Write(const std::string &o, collection* c);

    /*!
       \brief Write SDF formatted file
       \param os Output stream
       \param m molecule pointer
    */
    void           Write(std::ostream& os, molecule* m);

protected: // functions

    /*!
       \ Read a single line of the input file
       \param fileline Actual line
       \param pthisAtom AtomLine struct to temporary store atom info
    */
    bool           ReadAtomLine(std::string &fileline, AtomLine* pthisAtom);

    /*!
       \brief Get atom
       \param i atom index in itsAtomMap
    */
    int            getAtom(int i);

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

    //! Bond pointer
    Bond*          pBond;

    //! atom pointer used in bonding section
    atom*           pBondAtom1;

    //! atom pointer used in bonding section
    atom*           pBondAtom2;

    //! internal atom mapping local to molParser
    std::map<int, int>  itsAtomMap;

    // connections pointer
    connections* pConnections;

    //! molecule iterator
    typedef std::vector<molecule*>::iterator moleculeIterator;

    //! submolecule iterator
    typedef std::vector<submolecule*>::iterator subMoleculeIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator atomIterator;

    //! Bond Map Iterator
    typedef std::map<int, Bond*>::iterator BondMapIterator;

    //! Bond Map
    std::map<int, Bond*> moleculeBondMap;
};

} // MTKpp namespace

#endif // SDFPARSER_H
