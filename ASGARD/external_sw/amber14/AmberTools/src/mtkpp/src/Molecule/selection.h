/*! 
   \file selection.h
   \brief Allows for the selection of certain parts of the collection
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
   $Revision: 1.7 $

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

#ifndef SELECTION_H
#define SELECTION_H

namespace MTKpp
{

// ============================================================
// Class : selection()
// ------------------------------------------------------------
/*! 
   \class selection
   \brief Allows for the selection of certain parts of the collection
   \author Martin Peters
   \version 0.1
   \date 2006
*/
// ============================================================

class selection
{
public:

    /*!
       \brief selection Constructor
       \param parent collection pointer
    */
    selection(collection *parent = 0);

    //! selection Destructor
    virtual ~selection();

    /*!
       \brief Get collection
       \return collection pointer
    */
    collection*              getParent();

    /*!
       \brief parse the selection

       \param seln string to be parsed
       \return success

       \code
       The seln string should have the form:
       /col/mol/submol/atom
         ^    ^    ^     ^
         |    |    |     |++ atom name, number, or name-number
         |    |    |++ submolecule name, number, or name-number
         |    |++ molecule name, number, or name-number
         |++ collection name

       If seln begins with '/' it is assumed that the selection is starting from the 
       top of the structural hierarchy

       If seln doesn't begin with a '/' then it is assumed that the selection is from the
       bottom of the hierarchy up.
       \endcode
    */
    int                      parse(std::string seln);

    /*!
       \brief Get the molecule that's apart of the selection
       \return molecule pointer
    */
    molecule*                getMol();

    /*!
       \brief Get the submolecule that's apart of the selection
       \return submolecule pointer
    */
    submolecule*             getSMol();

    /*!
       \brief Get the atom that's apart of the selection
       \return atom pointer
    */
    atom*                    getAtom();

    /*!
       \brief Return atoms in selection
       \return atoms in selection
    */
    std::vector<atom*>       getAtoms();

    /*!
       \brief Return Selection Type
       \return selection type: 0 == collection, 1 == molecule, 2 == submolecule
    */
    int                      getSelectionType();

    /*!
       \brief Splits string
       \param s string to work on
       \param s2 separator
       \param v vector of string that gets returned
       \param i starting point
    */
    void splitString(std::string& s, const std::string s2, std::vector<std::string>& v, int i);

protected: // Functions
    /*!
       \brief Set molecule
       \param m molecule pointer
    */
    void                     setMol(molecule* m);

    /*!
       \brief Set submolecule
       \param s submolecule pointer
    */
    void                     setSMol(submolecule* s);

    /*!
       \brief Set atom
       \param a atom pointer
    */
    void                     setAtom(atom* a);

protected: // Data

    /*! 
       selection type definitions
       - 1 = molecule
       - 2 = submolecule
       - 3 = atom
    */
    int                      selectionType;

    //! parent collection
    collection*              pParent;

    //! seln molecule
    molecule*                selnMol;

    //! seln submolecule
    submolecule*             selnSMol;

    //! seln atom
    atom*                    selnAtom;

    //! selection name
    std::string              itsName;

    //! selection atoms
    std::vector<atom*>       itsAtoms;

};

} // MTKpp namespace

#endif // SELECTION_H
