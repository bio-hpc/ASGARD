/*!
   \file atomTyper.h
   \brief Atom Types molecules
   \author Martin Peters

   Assigns standard fragments to each residue/fragment

   Assigns standard atoms to each atom

   $Date: 2010/03/29 20:42:27 $
   $Revision: 1.11 $

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

#ifndef ATOMTYPER_H
#define ATOMTYPER_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;

class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;

class parameters;
struct atomType;

// ============================================================
// Class : atomTyper()
// ------------------------------------------------------------
// Class atomTyper
// ============================================================
class atomTyper
{
public:

    /*!
       \brief atomTyper constructor
    */
    atomTyper();

    /*!
       \brief atomTyper constructor
    */
    atomTyper(int b);

    /*!
       \brief atomTyper destructor
    */
    virtual ~atomTyper();

    /*!
       \brief Atom type by library
       \param c collection pointer
    */
    void atomTypeByLib(collection* c);

    /*!
       \brief Atom type by library
       \param m molecule pointer
    */
    void atomTypeByLib(molecule* m);

    /*!
       \brief Atom type by connectivity
       \param c collection pointer
    */
    void atomTypeByCon(collection* c);

    /*!
       \brief Assign Standard Properties
       \param m molecule pointer
    */
    void assignStdProperties(molecule* m);

protected:

    /*!
       \brief Perceive Histidine residues
    */
    void perceiveHistidines();

    /*!
       \brief Perceive Cysteine residues
    */
    void perceiveCysteines();

    /*!
       \brief Assigns all rings in a Molecule
       \param pMolecule Molecule pointer
    */
    void assignStdRings(molecule* pMolecule);

    /*!
       \brief Assigns all functional groups in a Molecule
       \param pMolecule Molecule pointer
    */
    void assignStdFunctionalGroups(molecule* pMolecule);

    /*!
       \brief Assigns all hydrophobic groups in a Molecule
       \param pMolecule Molecule pointer
    */
    void assignStdHydrophobicGroups(molecule* pMolecule);

    //! molecule iterator
    typedef std::vector<molecule*>::iterator moleculeIterator;

    //! submolecule iterator
    typedef std::vector<submolecule*>::iterator subMoleculeIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator atomIterator;

    //! molecule vector
    std::vector<molecule*>    moleculeList;

    //! submolecule vector
    std::vector<submolecule*> subMoleculeList;

    //! atom vector
    std::vector<atom*>        atomList;

    //! collection pointer
    collection*  pCollection;

    //! molecule pointer
    molecule*    pMolecule;

    //! submolecule pointer
    submolecule* pSubMolecule;

    //! atom pointer
    atom*        pAtom;

    //! standard library pointer
    stdLibrary* pStdLibrary;

    //! standard group pointer
    stdGroup*   pStdGroup;

    //! standard fragment pointer
    stdFrag*    pStdFrag;

    //! standard atom pointer
    stdAtom*    pStdAtom;

    //! parameters pointer
    parameters*  pParameters;

    //! Perceive Histidines boolean
    int bPerceiveHistidines;

    //! Perceive Cysteines boolean
    int bPerceiveCysteines;
};

} // MTKpp namespace

#endif // ATOMTYPER_H
