/*!
   \file hydrophobize.h
   \brief Determines the hydrophobic groups in a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:43:22 $
   $Revision: 1.4 $

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

#ifndef HYDROPHOBIZE_H
#define HYDROPHOBIZE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

namespace MTKpp
{

class molecule;
class atom;
class vector3d;
struct Bond;

// ============================================================
// Struct : hydrophobe
// ------------------------------------------------------------
/*!
   \struct hydrophobe
   \brief Container for hydrophobic group info
   \author Martin Peters
   \date 2007
*/
// ============================================================
struct hydrophobe
{
    //! Atoms in hydrophobic group
    std::vector<atom*> atoms;
};

// ============================================================
// Class : hydrophobize()
// ------------------------------------------------------------
/*! 
   \class hydrophobize

   \brief Determines the hydrophobic groups in a molecule

   \author Martin Peters

   \version 0.1

   \date 2007
*/
// ============================================================
class hydrophobize
{
public:
    /*!
       \brief hydrophobize Constructor
    */
    hydrophobize(molecule *parent);

    //! hydrophobize Destructor
    virtual ~hydrophobize();

    /*!
       \brief Determines the hydrophobic groups in a molecule
       \return success
    */
    int run();

protected: // DATA

    //---------------//
    // -  POINTERS - //
    //---------------//

    //! molecule pointer
    molecule*                pParent;

    //! atom pointer
    atom*                    pAtom;

    //! Bond pointer
    Bond*                    pBond;

    //! molecule atoms
    std::vector<atom*>       atomList;

    //! number of atoms
    unsigned int             nAtoms;

    //!
    std::vector<int>         atHydrophobicities;

    //! atomic number of each atom
    std::vector<int>         atNumbers;

    //! group which each atom is a part of
    std::vector<int>         atGroups;

    //! period which each atom is a part of
    std::vector<int>         atPeriods;

    //! atom symbols
    std::vector<std::string> atSymbols;

    //! atom hybridizations
    std::vector<int>         atHybridizations;

    //! formal charge of each atom
    std::vector<int>         formalCharges;

    //! bonding indices
    std::vector<std::vector<int> > bdAtoms;

    //! Highest bond order each atom contains
    std::vector<int>         bondOrders;

    //! bonds
    std::map<int, Bond*> bonds;

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator BondMapIterator;

    //! O-H groups
    std::vector<int>         OHs;

    //! N-H group
    std::vector<int>         NHs;

    //! S-H group
    std::vector<int>         SHs;

};

} // MTKpp namespace

#endif // HYDROPHOBIZE_H
