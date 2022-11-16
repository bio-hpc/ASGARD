/*!
   \file fingerPrint.h
   \brief Generates very simple molecular fingerprints
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

#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
struct Bond;
class rings;
class vector3d;

struct ring;

// ============================================================
// Class : fingerPrint()
// ------------------------------------------------------------
/*! 
   \class fingerPrint

   \brief Generates very simple molecular fingerprints

   \author Martin Peters

   \version 0.1

   \date 2006

   \section moleculeFingerPrint_sec Molecular Fingerprints
    There is 3 types of molecular fingerprints used in MTK++.
    - simple
    - fragment
    - 1-D string


   \section simpleFingerPrint_sec Simple Fingerprints
    The simple fingerprint is a very crude method of fingerprinting.

    The fingerprint encodes information about atoms, bonds, and rings.

    The vector looks like the following:
       [atom info, bond type, # of rings[ring info]]
       [H through to I, B-H, C-H, N-H, O-H, S-H, B-C, B=C, B-O, B-N, B-O, B-F, 
        B-S, B-Cl, B-Br, B-I, C-C, C=C, C%C, N-N, N=N, C-N, C=N, C%N, N-O, N=0, 
        N-P, N-Se, N=Se, O-O, C-O, C=O, O-Si, O-S, O=S, O-Se, O=Se, C-F, S-S, 
        C-S, C=S, S-N, C-Cl, P-P, P-C, P-O, P=O, P-S, P-Se, Se-Se, C-Se, C=Se, 
        N-Se, # rings*[size, number of hetero atoms, planarity, aromaticity]
       ]

    The length of the vector depends on the number of rings present in 
    the molecule or fragment.


   \section fragmentFingerPrints Fragment Finger Prints


   \section OneDFingerPrints 1-D Finger Prints


*/
// ============================================================
class fingerPrint
{
public:
    /*!
       \brief fingerPrint Constructor
    */
    fingerPrint();

    //! fingerPrint Destructor
    virtual ~fingerPrint();

    /*!
       \brief Generate simple fingerprint
       \param pMolecule molecule pointer
       \param fp fingerprint
    */
    void generateSimpleFP(molecule* pMolecule, std::vector<unsigned int> &fp);

    /*!
       \brief Compare simple fingerprint
       \param fp1 molecule fingerprint
       \param fp2 fragment fingerprint
       \return
          -# -1 error
          -#  0 no match
          -#  1 match
    */
    int compareSimpleFP(std::vector<unsigned int> &fp1, std::vector<unsigned int> &fp2);

    /*!
       \brief Generate fragment fingerprint
       \param pMolecule molecule pointer
       \param fp fingerprint
    */
    void generateFragmentFP(molecule* pMolecule, std::vector<unsigned int> fp);

    /*!
       \brief Generate 1D fingerprint
       \param pMolecule molecule pointer
       \param fp fingerprint
    */
    void generate1DFP(molecule* pMolecule, std::vector<unsigned int> fp);

protected:

    //! atom pointer
    atom*          pAtom1;

    //! atom pointer
    atom*          pAtom2;

    //! atom pointer
    atom*          pAtom3;

    //! atom pointer
    atom*          pAtom4;

    //! Bond pointer
    Bond*          pBond;

    //! Bond pointer
    Bond*          pBond1;

    //! Bond pointer
    Bond*          pBond2;
};

} // MTKpp namespace

#endif // FINGERPRINT_H
