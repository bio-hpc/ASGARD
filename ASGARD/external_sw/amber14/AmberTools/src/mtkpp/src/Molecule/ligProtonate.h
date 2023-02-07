/*!
   \file ligProtonate.h
   \brief Protonates small molecule
   \author Martin Peters

   $Date: 2010/03/29 20:44:26 $
   $Revision: 1.6 $

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

#ifndef LIGPROTONATE_H
#define LIGPROTONATE_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Utils/constants.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
struct Bond;
class parameters;
struct bondParam;
struct angleParam;
struct torsionParam;

class vector3d;

// ============================================================
// Class : ligProtonate()
// ------------------------------------------------------------
/*!
   \class ligProtonate

   \brief Class to add hydrogens to small molecules

   \author Martin Peters

   \version 0.1

   \date 2006

   \todo Definition of Hydrogen Bond Donors and Acceptor
     The perception of hydrogen bond donors and acceptors is performed in two steps:
     1. A simple default rule is applied to determine the donors and acceptors.
     2. Then, a list of exceptions to the default rules is applied.
     Default Rule
     - Any nitrogen or oxygen atom with at least one pair of free electrons is considered as an acceptor.
     - Any hydrogen atom that is bonded to a nitrogen or oxygen atom is considered as a donor.
     Exceptions
     - The nitrogen atoms of amides, sulfonamides, and equivalent structures are not considered as acceptors.
     - Only one nitrogen atom of guanidines, amidines, and equivalent groups is considered as an acceptor.
     - The nitrogen atoms of aromatic amines are not considered as acceptors.
     - The nitrogen atoms of five-membered aromatic rings like pyrrole are not considered as acceptors.
     - The sulfur atoms of thiocarbonyl compounds are considered as acceptors.

   \section protonateIntro_sec Protonation
   The functions contained in protonate adds hydrogen to unfilled atoms.


   \section protonateAlgorithm_sec Algorithm
    - Add Hydrogens to terminal atoms

    - Add Hydrogens to Ring systems

    - Add Hydrogens to the remaining atoms

   \section protonateDistance_sec Get Best Distance
     \code

        a--H

     \endcode
    - C-H = 1.09 Ang

    - N-H = 1.008 Ang

    - O-H = 0.95 Ang

    - S-H = 1.008 Ang

   \section protonateAngle_sec Get Best Angle
     \code

             b
            /
        H--a

     \endcode
     - a-b bond type
      -# single 109.47 deg
      -# double 120.0 deg
      -# triple 180.0 deg

     - a is a member of a ring
      -# angle = ((360 - ((ringSize-2)*180 )/ringSize)/2)
      -# 5 member ring --> angle equals 126 deg
      -# 6 member ring --> angle equals 120 deg

   \section protonateTorsion_sec Get Best Torsion
     \code

             b--c
            /
        H--a
            \
             d

     \endcode
     - a-b is not a ring bond
      -# a sp3 hybridized [60, 180, 300].
        180 if a is terminal. Test which dihedrals are available and add accordingly
      -# a sp2 [0, 180]
      -# a sp [180]

     - a-b is a ring bond
      -# aromatic [180]
      -# nonaromatic
        - sp3, sp2, sp are same are non ring atoms
*/
// ============================================================
class ligProtonate
{
public:

    /*!
       \brief ligProtonate Constructor
    */
    ligProtonate();

    //! ligProtonate Destructor
    virtual ~ligProtonate();

    /*!
       \brief add hydrogens
       \param s submolecule pointer
    */
    void addHydrogens(submolecule* s);

protected:

    /*!
       \brief Determine the best distance
       \param a atom pointer
       \return The best distance
    */
    double getBestDistance(atom* a);

    /*!
       \brief Determine the best angle
       \param a atom pointer
       \param b atom pointer
       \param r ring size
       \return The best available angle
    */
    double getBestAngle(atom* a, atom* b, const int& r = 0);

    /*!
       \brief Determine the best torsion
       \param a atom pointer
       \param b atom pointer
       \param c atom pointer
       \return The best available torsions
    */
    std::vector<double> getBestTorsions(atom* a, atom* b, atom* c);

    /*!
       \brief Determine the best torsion to add hydrogens to Hbond donors

       \verbatim
        D--H
            .
             .
              A--AA

         E_hb = - cos^2(A_DHA) * e^(-(R_HA-2.0)^2)  --> from andrew

         AA-A ... D > 90 deg
         D ... A < 3.5 ang
         D-H ... A > 90 deg
       \endverbatim

       \param a atom pointer
       \param b atom pointer
       \param c atom pointer
       \param dist distance b/w H atom and the donor to which it's being added
       \param angle
       \param acceptors vector of all hbond acceptors in the molecule
       \return The best available torsions
    */
    std::vector<double> getHBDonorTorsions(atom* a, atom* b, atom* c,
                double dist, double angle, std::vector<atom*> acceptors);

protected: // DATA
    //! collection pointer
    collection*    pCol;

    //! collection boolean
    bool           bCol;

    //! molecule pointer
    molecule*      pMol;

    //! molecule boolean
    bool           bMol;

    //! submolecule list
    std::vector<submolecule*> subMoleculeList;

    //! atom list
    std::vector<atom*>        atomList;

    //! submolecule pointer
    submolecule*   pSubMolecule;

    //! submolecule pointer
    submolecule*   pSubMoleculeMinus1;

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

    //! vector3d pointer
    vector3d*      coord1;

    //! vector3d pointer
    vector3d*      coord2;

};

} // MTKpp namespace

#endif // LIGPROTONATE_H

