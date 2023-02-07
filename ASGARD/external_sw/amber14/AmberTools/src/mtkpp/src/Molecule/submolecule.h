/*!
   \file submolecule.h
   \brief Container for atoms and bonds
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
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

#ifndef SUBMOLECULE_H
#define SUBMOLECULE_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

#include "Utils/constants.h"

namespace MTKpp
{

class collection;
class molecule;
class atom;
class vector3d;

struct Bond;
struct Angle;
struct Torsion;
struct Improper;

class stdFrag;
struct stdAtom;

// ============================================================
// Class : submolecule()
// ------------------------------------------------------------
/*! 
   \class submolecule
   \brief Container for atoms and bonds
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class submolecule
{
    friend class molecule;
public:

    /*!
       \brief submolecule Constructor
       \param parent molecule pointer
    */
    submolecule(molecule *parent = 0);

    //! submolecule Destructor
    virtual ~submolecule();

    /*!
       \brief Compares two submolecule based on collection index
       \param lhs first idObject
       \param rhs Second idObject
       \return boolean

       After this function is defined the STL sort() function can be used
    */
    static bool less(submolecule* lhs, submolecule* rhs)
    {
        return lhs->getColIndex() < rhs->getColIndex();
    }

    /*!
       \brief Compares two submolecule based on collection index
       \param lhs first idObject
       \param rhs Second idObject
       \return boolean

       After this function is defined the STL sort() function can be used
    */
    static bool greater(submolecule* lhs, submolecule* rhs)
    {
        return lhs->getColIndex() > rhs->getColIndex();
    }

    /*!
       \brief Copy contents of pSubMol in this submolecule
       \param pSubMol submolecule pointer
    */
    void                     copy(submolecule* pSubMol);

    /*!
       \brief Get molecule which submolecule is apart of
       \return molecule pointer
    */
    molecule*                getParent();

    /*!
       \brief Add atom to the submolecule
       \return atom pointer
    */
    virtual atom*            addAtom();

    /*!
       \brief Copy atom into the submolecule
       \param pAt atom pointer to be copied
       \return atom pointer
    */
    virtual atom*            addAtom(atom* pAt);

    /*!
       \brief Add a Bond between a and b
       \param a atom pointer 1
       \param b atom pointer 2
       \return Bond pointer
    */
    virtual Bond*            addBond(atom* a, atom* b);

    /*!
       \brief Set name of submolecule
       \param name submolecule name
    */
    void                     setName(const std::string& name);

    /*!
       \brief Set 1-Letter name of submolecule
       \param name submolecule 1-Letter name
    */
    void                     set1LName(const std::string& name);

    /*!
       \brief Get name of submolecule
       \return submolecule name
    */
    std::string              getName();

    /*!
       \brief Get 1-Letter name of submolecule
       \return submolecule 1-Letter name
    */
    std::string              get1LName();

    /*!
       \brief Set index of submolecule
       \param i submolecule index
    */
    void                     setIndex(const int& i);

    /*!
       \brief Get index of submolecule
       \return submolecule index
    */
    int                      getIndex();

    /*!
       \brief Set id of submolecule
       \param id submolecule id
    */
    void                     setSubMolId(const int& id);

    /*!
       \brief Get id of submolecule
       \return submolecule id
    */
    int                      getSubMolId();

    /*!
       \brief Set collection index of submolecule
       \param id submolecule collection index
    */
    void                     setColIndex(const int& id);

    /*!
       \brief Get collection index of submolecule
       \return submolecule collection index
    */
    int                      getColIndex();

    /*!
       \brief Set icode of submolecule
       \param i submolecule icode
    */
    void                     setiCode(const std::string& i);

    /*!
       \brief Get icode of submolecule
       \return submolecule icode
    */
    std::string             getiCode();
    
    /*!
       \brief Set the number of atoms in submolecule
       \param n number of atoms
    */
    void                     setNumAtoms(const int& n);

    /*!
       \brief Get number of atoms in the submolecule
       \return number of atoms in submolecule
    */
    int                      getNumAtoms();

    /*!
       \brief Set the number of bonds in submolecule
       \param n number of bonds
    */
    void                     setNumBonds(const int& n);

    /*!
       \brief Get number of bonds in the submolecule
       \return number of bonds in submolecule
    */
    int                      getNumBonds();

    /*!
       \brief Set the standard fragment for this submolecule
       \param f standard fragment (stdFrag) pointer
    */
    void                     setStdFrag(stdFrag* f);

    /*!
       \brief Get standard fragment associated with this submolecule
       \return stdFrag pointer
    */
    stdFrag*                 getStdFrag();

    /*!
       \brief Has a standard fragment
       \return boolean
    */
    bool                     hasStdFrag();

    /*!
       \brief Get all atoms in submolecule
       \return list to atoms
    */
    std::vector<atom*>       getAtomList();

    /*!
       \brief Get all bonds in submolecule
       \return list to bonds
    */
    std::vector<Bond*>       getBondList();

    /*!
       \brief Get an atom in submolecule
       \param i integer index
       \return atom pointer
    */
    atom*                    getAtom(int i);

    /*!
       \brief Get an atom in submolecule
       \param name atom name
       \return atom pointer
    */
    atom*                    getAtom(const std::string& name);

    /*!
       \brief Get atom
       \param pStdAtom standard atom pointer
       \return atom pointer
    */
    atom*                    getAtom(stdAtom* pStdAtom);

    //! Calculate the center of mass
    void                     centerOfMass();

    //! Get the center of mass
    vector3d*                getCenterOfMass();

    /*!
       \brief Get number of stdAtoms in fragment
       \return number of stdAtoms
    */
    int                      numHeavyAtoms();

    /*!
       \brief Get Molecular Weight
       \return molecular weight
    */
    double                   getMolecularWeight();

    /*!
       \brief Print submolecule details to screen
    */
    void                     print();

protected:

    //! submolecules' molecule
    molecule*                pParent;

    //! atom iterator
    typedef std::vector<atom*>::iterator AtomIterator;

    //! submolecule name
    std::string              itsName;

    //! submolecule name
    std::string              its1LName;

    //! iCode from pdb file
    std::string              itsiCode;

    //! Collection index
    int                      itsColIndex;

    //! Molecule index
    int                      itsIndex;

    //! submolecule internal id
    int                      itsSubMolId; // This is a file id when pdb files are read

    //! number of atoms
    int                      itsNumAtoms;

    //! number of bonds
    int                      itsNumBonds;

    // coordinates of its center of mass
    vector3d*                pCenterMass;

    //! list of atoms in submolecule
    std::vector<atom*>       itsAtomList;

    //! list of bonds in submolecule
    std::vector<Bond*>       itsBondList;

    //! atom pointer
    atom*                    pAtom;

    //! Bond pointer
    Bond*                    pBond;

    //! standard fragment pointer
    stdFrag*                 pStdFrag;
};

} // MTKpp namespace

#endif // SUBMOLECULE_H

