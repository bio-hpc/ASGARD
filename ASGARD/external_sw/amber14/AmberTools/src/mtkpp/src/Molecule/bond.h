/*!
   \file bond.h
   \brief Container for bond information
   \author Martin Peters

   $Date: 2010/03/29 20:42:27 $
   $Revision: 1.10 $

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

#ifndef BOND_h
#define BOND_h

#include <iostream>
#include <string>

namespace MTKpp
{

class atom;
struct bondParam;

// ============================================================
// Struct : Bond
// ------------------------------------------------------------
/*! 
   \struct Bond

   \brief Container for Bond info

   \author Martin Peters

   \date 2005

    \section BondTypeDefinitions Bond Type Definitions
      - 0 = Undefined
      - 1 = Single
      - 2 = Double
      - 3 = Triple
      - 4 = Aromatic
      - 5 = Single or Double
      - 6 = Single or Aromatic
      - 7 = Double or Aromatic
      - 8 = Any type

    \section BondStereoDefinitions Bond Stereo Definitions
      - Definitions for Single Bonds
       -# 0 = Not stereo
       -# 1 = Up
       -# 4 = Either
       -# 6 = Down
      - Definitions for Double Bonds
        -# 0 = Use x,y,z coords from atom block to determine cis or trans
        -# 3 = Either cis or trans

    \section BondTopologyDefinitions Bond Topology Definitions
      - 0 = Either
      - 1 = Ring
      - 2 = Chain

    \section BondKindDefinitions Bond Kind Definitions
      - 0 = Undefined
      - 1 = Polar
*/
// ============================================================
struct Bond
{
    //! First Atom in Bond
    atom* atom1;

    //! Second Atom in Bond
    atom* atom2;

    /*!
        Bond Type Definitions
        - 0 = Undefined
        - 1 = Single
        - 2 = Double
        - 3 = Triple
        - 4 = Aromatic
        - 5 = Single or Double
        - 6 = Single or Aromatic
        - 7 = Double or Aromatic
        - 8 = Any type
    */
    int type;

    /*!
        Bond Type Definitions
        - 0 = Undefined
        - 1 = Single
        - 2 = Double
        - 3 = Triple
        - 6 = NO2 bond
        - 7 = Aromatic and single
        - 8 = Aromatic and double
        - 9 = CO2 or CS2 double
        - 10 = other Aromatic
    */
    int bccType;

    /*!
        Bond Stereo
        - Definitions for Single Bonds
         -# 0 = Not stereo
         -# 1 = Up
         -# 4 = Either
         -# 6 = Down
        - Definitions for Double Bonds
         -# 0 = Use x,y,z coords from atom block to determine cis or trans
         -# 3 = Either cis or trans
    */
    int stereo;

    /*!
        Bond Topology Definitions
        - 0 = Either
        - 1 = Ring
        - 2 = Chain
    */
    int topology;

    /*!
        Bond Kind Definitions
        - 0 = Undefined
        - 1 = Polar
        - 2 = Amide
        - 3 = disulfide
    */
    int kind;

    //! Length of Bond
    double size;

    //! Pointer to Bond Parameters
    bondParam* pBondParam;

    /*!
       Rotatable bonds Definitions
       - 0 Not Rotatable
       - 1 Rotatable
       - 2 Rotatable but C2 Symmetric
    */
    int rotatable;
};

} // MTKpp namespace

#endif // BOND_H


