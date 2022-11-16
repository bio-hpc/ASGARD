/*!
   \file proProtonate.h
   \brief Protonates proteins
   \author Martin Peters
   \author Andrew Wollacott

   $Date: 2010/03/29 20:44:27 $
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

#ifndef PROPROTONATE_H
#define PROPROTONATE_H

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

class vector3d;

class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;
struct stdBond;
struct stdImproper;
struct stdLoop;
struct stdAlias;

class parameters;
struct bondParam;
struct angleParam;
struct torsionParam;

// ============================================================
// Class : proProtonate()
// ------------------------------------------------------------
/*!
   \class proProtonate

   \brief Class to add hydrogens to proteins

   \author Martin Peters

   \version 0.1

   \date 2006
*/
// ============================================================
class proProtonate
{
    friend class protonate;
public:

    /*!
       \brief protonate Constructor
    */
    proProtonate();

    //! proProtonate Destructor
    virtual ~proProtonate();

    /*!
       \brief add hydrogens by library
       \param s1 submolecule 1 pointer
       \param f1 standard residue 1
    */
       void addHydrogens(submolecule* s1, stdFrag* f1);

protected:

    /*!
       \brief Builds the missing atoms in the first residue
       \param pSubMolecule submolecule pointer
       \param missingChainAtoms The missing standard atoms which need to be built
       \param first3Atoms The first 3 main chain atoms in the molecule
    */
    void buildMissingAtoms(submolecule* pSubMolecule,
                           std::vector<stdAtom*> missingChainAtoms,
                           std::vector<atom*> first3Atoms);

    /*!
       \brief Builds the 3 dummy atom for protonating the first standard residue
       \param first3Atoms The first 3 main chain atoms in the molecule
    */
    void buildDummyAtoms(std::vector<atom*> first3Atoms);

    /*!
       \brief Determine the best torsion
       \param pStdFrag stdFrag pointer
       \param stdAt1 stdAtom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \return The best available torsions
    */
    std::vector<double> getBestTorsions(stdFrag* pStdFrag, stdAtom* stdAt1, atom* at2, atom* at3, atom* at4);

    /*!
       \brief Optimize polar hydrogens
    */
    void optimizePolarHs();

protected: // DATA
    //! collection pointer
    collection*    pCol;

    //! collection boolean
    bool bCol;

    //! molecule pointer
    molecule*      pMol;

    //! molecule boolean
    bool bMol;

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

    //! stdLibrary pointer
    stdLibrary*    pStdLibrary;

    //! stdGroup pointer
    stdGroup*      pStdGroup;

    //! stdFrag pointer
    stdFrag*       pStdFrag;

    //! stdFrag pointer
    stdFrag*       pStdFragMinus1;

    //! stdAtom pointer
    stdAtom*       pStdAtom;

    //! stdAtom pointer
    stdAtom*       pStdAtom1;

    //! stdAtom pointer
    stdAtom*       pStdAtom2;

    //! stdAtom pointer
    stdAtom*       pStdAtom3;

    //! stdAtom pointer
    stdAtom*       pStdAtom4;

    //! dummy atom 1
    atom*          dummyAtom1;

    //! dummy atom 2
    atom*          dummyAtom2;

    //! dummy atom 3
    atom*          dummyAtom3;

    //! parameters
    parameters*    pParam;

    //!
    std::vector<atom*> prev3Atoms;
};

} // MTKpp namespace

#endif // PROPROTONATE_H
