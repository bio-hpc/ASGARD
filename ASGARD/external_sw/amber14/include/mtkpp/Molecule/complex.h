/*!
   \file complex.h
   \brief Container for a receptor and ligand
   \author Martin Peters

   $Date: 2010/03/29 20:42:28 $
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

#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "Utils/constants.h"

namespace MTKpp
{
class collection;
class molecule;
class submolecule;
class atom;
class vector3d;

class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;

class elements;

struct funcGroup;
// ============================================================
// Class : complex()
// ------------------------------------------------------------
/*!
   \class complex
   \brief complex
   \author Martin Peters
*/
// ============================================================

class complex
{
    friend class collection;
public:

    /*!
       \brief complex Constructor
       \param pCollection collection pointer
       \param bLigOnly Only consider ligand-ligand interactions
    */
    complex(collection* pCollection, bool bLigOnly);

    //! complex Destructor.
    virtual ~complex();

    /*!
       \brief Get number of atoms
    */
    int getNumAtoms();

    /*!
       \brief Get number of residues
    */
    int getNumResidues();

    /*!
       \brief Get number of atoms in the receptor
    */
    int getNumRecAtoms();

    /*!
       \brief Get number of residues in the receptor
    */
    int getNumRecResidues();

    /*!
       \brief Get number of atoms in the ligand
    */
    int getNumLigAtoms();

    /*!
       \brief Get number of residues in the ligand
    */
    int getNumLigResidues();

    /*!
       \brief Get size of the atom interactions matrix
    */
    int getAtomsMatrixSize();

    /*!
       \brief Get size of the residue interactions matrix
    */
    int getResMatrixSize();

    /*!
       \brief Only store the interactions involing ligand atoms
    */
    bool getLigOnly();

    /*!
       \brief Get the number of receptor-receptor atom interactions
    */
    int getRecRecAtoms();

    /*!
       \brief Get the number or receptor-receptor residue interactions
    */
    int getRecRecResidues();

    /*!
       \brief Get receptor flags arrary
       \param recFlags ligand flags array
    */
    int getRecFlags(int recFlags[]);

    /*!
       \brief Get ligand flags arrary
       \param ligFlags ligand flags array
    */
    int getLigFlags(int ligFlags[]);

    /*!
       \brief Is the atom in a ligand?
       \param a atom index
    */
    bool isAtomLig(int a);

    /*!
       \brief Get new index for atom
       \param a atom index
    */
    int getAtomIndex(int a);

    /*!
       \brief Is the residue in a ligand?
       \param r residue index
    */
    bool isResLig(int r);

    /*!
       \brief Get new index for residue
       \param r residue index
    */
    int getResIndex(int r);

    /*!
       \brief Set the core fragment in the ligand
       \param coreFrag core fragment name
       \return success
    */
    int setCoreFragment(std::string coreFrag);

    /*!
       \brief Set the mapping
       \param m1 internal code
       \param m2 map code
       \return success
    */
    void addMapping(std::string m1, std::string m2);

    /*!
       \brief Get residue mapping
       \param m1 internal code
       \return map code
    */
    std::string getMapping(std::string m1);

    /*!
       \brief Get residue name and id
       \param i index
       \return internal code
    */
    std::string getResNameID(int i);

    /*!
       \brief Set residue vectors
       \return vector of vectors of atom indices
    */
    std::vector<std::vector<int> > getResidues();

    /*!
       \brief Get ligand fragments
       \return vector of ligand fragment indices
    */
    std::vector<int> getLigandFragments();

    /*!
       \brief Get receptor residues
       \return vector of receptor submolecule
    */
    std::vector<submolecule*> getReceptorSubMols();

    /*!
       \brief Get ligand residues
       \return vector of ligand submolecule
    */
    std::vector<submolecule*> getLigandSubMols();

    /*!
       \brief Is the atom in a chain?
       \param c chain index
       \param a atom index
    */
    bool isChain(int c, int a);

protected: // Functions

    /*!
       \brief Define receptor, ligand, and solvent molecules
    */
    void setup();

protected: // Data

    //! Name of complex
    collection*                                  pCollection;

    //! Total number of atoms
    int                                          nAtoms;

    //! Number of atoms in receptor
    int                                          nRecAtoms;

    //! Number of atoms in ligand
    int                                          nLigAtoms;

    //! Size of atom matrices
    int                                          atomMatrixSize;

    //! Total Number of residues
    int                                          nResidues;

    //! Number of residues in receptor
    int                                          nRecResidues;

    //! Number of residues in ligand
    int                                          nLigResidues;

    //! Size of residue matrices
    int                                          resMatrixSize;

    //! Do only lig-lig and rec-lig interactions
    bool                                         bLigOnly;

    //! Number of rec-rec atom interactions
    int                                          nRecRecAtoms;

    //! Number of rec-rec residue interactions
    int                                          nRecRecResidues;

    /*
          rec = [NATOMS] = [1,1,1,0,0,0,0]
          lig = [NATOMS] = [0,0,0,1,1,0,0]
          cof = [NATOMS] = [0,0,0,0,0,0,1]
          sol = [NATOMS] = [0,0,0,0,0,1,0]
       recSol = [NATOMS] = [1,1,1,0,0,1,1]
      indices = [NATOMS] = [1,2,3,5,6,4,7]
        chain = [NATOMS] = [1,1,1,2,2,2,2]
    */

    //! receptor atoms
    int                                          *rec;

    //! ligand atoms
    int                                          *lig;

    //! cofactor atoms
    int                                          *cof;

    //! solvent atoms
    int                                          *sol;

    //! receptor and solvent atoms
    int                                          *recSol;

    //! receptor and solvent atom indices
    int                                          *atomIndices;

    //! receptor residues
    int                                          *resRec;

    //! ligand residues
    int                                          *resLig;

    //! solvent residues
    int                                          *resSol;

    //! receptor and solvent residues
    int                                          *resRecSol;

    //! receptor and solvent residue indices
    int                                          *resIndices;

    //! chain counter
    int                                          chainCounter;

    //! chain indices
    int                                          *chain;

    //! Vector of receptor residues
    std::vector<submolecule*>                    receptor;

    //! ligand molecule
    molecule*                                    pLigand;

    //! Vector of ligand fragments
    std::vector<submolecule*>                    ligand;

    //! Vector of cofactor fragments
    std::vector<submolecule*>                    cofactor;

    //! Vector of solvent residues
    std::vector<submolecule*>                    solvent;

    //! Vector of atom indices
    std::vector<std::vector<int> >               residues;

    //! Indices pointing to ligand fragments in the residues vector
    std::vector<int>                             ligandIndices;

    //! string map iterator
    typedef std::map<std::string, std::string>::iterator stringMapIterator;

    //! vector of residue names and ids
    std::vector<std::string>                     residueNAMEIDs;

    //! Map of residue names and ids
    std::map<std::string, std::string>           residueMap;

    //! cofactor names
    std::vector<std::string>                     cofactors;
};

} // MTKpp namespace

#endif // COMPLEX_H
