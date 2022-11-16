/*!
   \file protonate.h
   \brief Protonates a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:44:27 $
   $Revision: 1.16 $

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

#ifndef PROTONATE_H
#define PROTONATE_H

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

class proProtonate;
class ligProtonate;
class watProtonate;

// ============================================================
// Class : protonate()
// ------------------------------------------------------------
/*!
   \class protonate

   \brief Class to add hydrogens
*/
// ============================================================
class protonate
{
public:

    /*!
       \brief protonate Constructor
       \param col collection pointer
    */
    protonate(collection *col = 0);

    /*!
       \brief protonate Constructor
       \param mol molecule pointer
    */
    protonate(molecule *mol = 0);

    //! protonate Destructor
    virtual ~protonate();

    /*!
       \brief Adds hydrogens
    */
    void run();

    /*!
       \brief Optimize polar hydrogen in a collection
    */
    void optimizePolarHs();

protected: // FUNCTIONS

    /*!
       \brief Initialize
    */
    void initialize();

    /*!
       \brief Adds hydrogens
    */
    void runCol();

    /*!
       \brief Adds hydrogens
    */
    void runMol(molecule* pMolecule);

    /*
      \brief Attempt to build the missing heavy atoms of a fragment
      \param pSubMolecule submolcule pointer
    */
    int buildMissingHeavyAtoms(submolecule* pSubMolecule);

protected: // DATA
    //! proProtonate
    proProtonate*  pPro;

    //! ligProtonate
    ligProtonate*  pLig;

    //! watProtonate
    watProtonate*  pWat;

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
    std::vector<atom*> atomList;

    //! submolecule pointer
    submolecule*   pSubMolecule;

    //! submolecule pointer
    submolecule*   pSubMoleculeMinus1;

    //! stdLibrary pointer
    stdLibrary*    pStdLibrary;

    //! stdGroup pointer
    stdGroup*      pStdGroup;

    //! stdFrag pointer
    stdFrag*       pStdFrag;

    //! stdFrag pointer
    stdFrag*       pStdFragMinus1;

    //! parameters
    parameters*    pParam;

    //!
    std::vector<atom*> prev3Atoms;
};

} // MTKpp namespace

#endif // PROTONATE_H

