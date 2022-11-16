/*!
   \file collection.h
   \brief Container for molecules
   \author Martin Peters

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.20 $

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

#ifndef COLLECTION_H
#define COLLECTION_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "Utils/constants.h"

namespace MTKpp
{

class molecule;
class submolecule;
class atom;
class vector3d;

class elements;

struct Bond;
struct Angle;
struct Torsion;
struct Improper;

// Molecular Mechanics
class stdLibrary;
class parameters;
class atomTyper;
struct bondParam;
struct angleParam;
struct torsionParam;
struct improperParam;
struct hBondParam;
struct LJ612Param;
struct LJ612SE;

// metalloproteins
class metalCenter;
class metalGroup;

// ============================================================
// Class : collection()
// ------------------------------------------------------------
/*!
   \class collection
   \brief Container for molecules
*/
// ============================================================

class collection
{
    friend class molecule;
public:

    //! Collection Constructor
    collection();

    //! Collection Destructor.
    virtual ~collection();

    /*!
       \brief Set collection name
       \param name collection name
    */
    void                     setName(std::string name);

    /*!
       \brief Get collection name
       \return collection name
    */
    std::string              getName();

    /*!
       \brief Add Molecule to Collection
       \return molecule pointer
    */
    virtual molecule*        addMolecule();

    /*!
       \brief Delete Molecule from Collection
       \param pMol molecule pointer
    */
    void                     delMolecule(molecule* pMol);

    /*!
       \brief Delete all Molecules from Collection
    */
    void                     delAllMolecules();

    /*!
       \brief Clears Collection
    */
    void                     clear();

    /*!
       \brief Get last added Molecule from list of molecules
       \return molecule pointer
    */
    molecule*                getMolecule();

    /*!
       \brief Get Molecule from list of molecules based on id
       \param id molecule id
       \return molecule pointer
       \sa getMoleculeList()
    */
    molecule*                getMolecule(int id);

    /*!
       \brief Get Molecule from list of molecules based on name
       \param name molecule name
       \return molecule pointer
       \sa getMoleculeList()
    */
    molecule*                getMolecule(std::string name);

    /*!
       \brief Get last added Molecule from list of molecules
       \return molecule pointer
    */
    molecule*                getLastAddedMolecule();

    /*!
       \brief Get list of molecules
       \return vector of molecule pointers
    */
    virtual std::vector<molecule*>   getMoleculeList();

    /*!
       \brief Get atom
       \param number integer value
       \param atomIndex return atom based on index (default)
       \param fileId return atom based on file id
       \param atomColIndex atom collection index
       \return atom pointer
    */
    atom*                    getAtom(int number, bool atomIndex = true,
                                                 bool fileId = false,
                                                 bool atomColIndex = false);

    /*!
       \brief Get list of atoms
       \return vector of atom pointers
    */
    virtual std::vector<atom*> getAtomList();

    /*!
       \brief Get number of molecules
       \return number of molecules in the collection
    */
    int                      getNumberMolecules();

    /*!
       \brief Get number of submolecules
       \return number of submolecules in the collection
    */
    int                      getNumberSubMolecules();

    /*!
       \brief Get number of atoms in collection
       \return number of atoms in the collection
    */
    int                      getNumAtoms();

    ////  BONDS ////

    /*!
       \brief Get number of bonds in collection
       \return number of bonds in the collection
    */
    int                      getNumBonds();

    /*!
       \brief Get bond list
       \param bonds array of atom indices
    */
    int                      getBonds(int bonds[]);

    /*!
       \brief Get bond parameters
       \param bondParams array of bond parameters
    */
    int                      getBondParams(double bondParams[]);

    /*!
       \brief Get number of bonds that contain a Hydrogen
    */
    int                      getNumBondsWithH();

    /*!
       \brief Get number of bonds that do not contain a Hydrogen
    */
    int                      getNumBondsWithOutH();

    /*!
       \brief Get number of unique bond types used
    */
    int                      getNumUniqueBondTypes();

    ////  ANGLES ////

    /*!
       \brief Get number of angles in collection
       \return number of angles in the collection
    */
    int                      getNumAngles();

    /*!
       \brief Get angle list
       \param angles array of atom indices
    */
    int                      getAngles(int angles[]);

    /*!
       \brief Get angle parameters
       \param angleParams array of angle parameters
    */
    int                      getAngleParams(double angleParams[]);

    /*!
       \brief Get number of angles that contain a Hydrogen
    */
    int                      getNumAnglesWithH();

    /*!
       \brief Get number of angles that do not contain a Hydrogen
    */
    int                      getNumAnglesWithOutH();

    /*!
       \brief Get number of unique angle types used
    */
    int                      getNumUniqueAngleTypes();

    ////  TORSIONS ////

    /*!
       \brief Get number of torsions in collection
       \return number of torsions in the collection
    */
    int                      getNumTorsions();

    /*!
       \brief Get number of MM torsions in collection
       \return number of torsions in the collection
    */
    int                      getNumMMTorsions();

    /*!
       \brief Get torsion atom indices and parameters
       \param torsions array of atom indices
       \param torsionParams array of torsion parameters
    */
    int                      getMMTorsions(int torsions[], double torsionParams[]);

    ////  IMPROPERS ////

    /*!
       \brief Get number of impropers in collection
       \return number of impropers in the collection
    */
    int                      getNumImpropers();

    /*!
       \brief Get number of MM impropers in collection
       \return number of impropers in the collection
    */
    int                      getNumMMImpropers();

    /*!
       \brief Get improper atom indices and parameters
       \param impropers array of atom indices
       \param improperParams array of torsion parameters
    */
    int                      getMMImpropers(int impropers[], double improperParams[]);

    /*!
       \brief Get number of dihedrals that contain a Hydrogen
    */
    int                      getNumDihedralsWithH();

    /*!
       \brief Get number of dihedrals that do not contain a Hydrogen
    */
    int                      getNumDihedralsWithOutH();

    /*!
       \brief Get number of unique dihedral types used
    */
    int                      getNumUniqueDihedralTypes();

    ////  COORDINATES ////

    /*!
       \brief Get all atomic coordinates
       \param coords array of coordinates
    */
    int                      getCoordinates(double coords[]);

    /*!
       \brief Set all atomic coordinates
       \param coords array of coordinates
    */
    int                      setCoordinates(double coords[]);

    //! elements object pointer
    elements*                pElements;

    /////////////////////////////
    // - Molecular Mechanics - //
    /////////////////////////////

    //! Add standard library to the collection
    virtual void             addStdLibrary();

    //! Get Standard Library
    virtual stdLibrary*      getStdLibrary();

    //! Add atom typer object to the collection
    void                     addAtomTyper();

    //! Get atom typer object
    atomTyper*               getAtomTyper();

    //! Add mm parameters object to the collection
    void                     addParameters();

    //! Get mm parameters object
    parameters*              getParameters();

    //! Get number of MM non bonded interactions
    int                      getNumMMnonBondedPairs(double cutoff);

    /*!
       \brief Get MM non bonded pairs
    */
    int                      getMMnonBondedPairs(int nonBonded[],
                              int nonBondedPtrs[], double nonBondedParams[],
                              int nonBonded14Ptrs[], double cutoff);

    /*!
       \brief Get all atomic charges
       \param charges array of charges
       \return success
    */
    int                      getMMCharges(double charges[]);

    /*!
       \brief Get all atomic symbols
       \param symbols char array of symbols
       \return success
    */
    int                      getAtomSymbols(char symbols[]);

    /*!
       \brief Get all atomic names
       \param names char array of names
       \return success
    */
    int                      getAtomNames(char names[]);

    /*!
       \brief Get all atomic masses
       \param masses char array of masses
       \return success
    */
    int                      getAtomMasses(double masses[]);

    /*!
       \brief Get number of unique atom types used
       \return number of unique atom types used
    */
    int                      getNumUniqueAtomTypes();

    /*!
       \brief Get number of unique atom types used
       \return number of unique atom types used
    */
    std::vector<std::string> getUniqueAtomTypes();

    /*!
       \brief Get number of unique bond types used
       \return number of unique bond types used
    */
    std::vector<bondParam*> getUniqueBondTypes();

    /*!
       \brief Get number of unique angle types used
       \return number of unique angle types used
    */
    std::vector<angleParam*> getUniqueAngleTypes();

    /*!
       \brief Get number of unique torsion types used
       \return number of unique torsion types used
    */
    std::vector<torsionParam*> getUniqueTorsionTypes();

    /*!
       \brief Get number of unique improper types used
       \return number of unique improper types used
    */
    std::vector<improperParam*> getUniqueImproperTypes();

    /*!
       \brief Get numeric atom types
       \param a atom types array
       \return success
    */
    int                      getAtomTypes(int a[]);

    /*!
       \brief Get numeric atom types
       \param a atom types array
       \return success
    */
    int                      getAtomTypes(char a[]);

    /*!
       \brief Get L-J Params
       \param r6 L-J R6 array
       \param r12 L-J R12 array
       \return success
    */
    int                      getLJParams(double r6[], double r12[]);

    /*!
       \brief Get number of excluded atoms per atom
       \return success
    */
    int                      getNumExcludedAtoms();

    /*!
       \brief Get number of excluded atoms per atom
       \param e excluded atom list
       \return success
    */
    int                      getNumExcludedAtoms(int e[]);

    /*!
       \brief Get excluded atom list
       \param e excluded atom list
       \return success
    */
    int                      getExcludedAtoms(int e[]);

    /*!
       \brief Get number of excluded atoms per atom
       \return success
    */
    int                      getNumExcluded14Atoms();

    /*!
       \brief Get number of excluded atoms per atom
       \param e excluded atom list
       \return success
    */
    int                      getNumExcluded14Atoms(int e[]);

    /*!
       \brief Get excluded atom list
       \param e excluded atom list
       \return success
    */
    int                      getExcluded14Atoms(int e[]);

    /*!
       \brief Get all residue names
       \param resNames char array of residue names
       \return success
    */
    int                      getResidueNames(char resNames[]);

    /*!
       \brief Get all residue pointers
       \param resPointers int array of residue pointers
       \return success
    */
    int                      getResiduePointers(int resPointers[]);

    /*!
       \brief Set atom index
       \param n atom index
    */
    void                     setAtomIndex(const int& n);

    /*!
       \brief Get atom index
       \return atom index
    */
    int                      getAtomIndex();

    /*!
       \brief Set submolecule index
       \param n index
    */
    void                     setSubMoleculeIndex(const int& n);

    /*!
       \brief Get submolecule index
       \return n index
    */
    int                      getSubMoleculeIndex();

    /*!
       \brief Get total formal charge
       \return formal charge
    */
    virtual int              getFormalCharge();

    /*!
       \brief Get number of nonbonded neighbors within a certain distance
       \param pAt atom pointer
       \param d distance
       \return number of nonbonded neighbors
    */
    int                      getNumNeighbors(atom* pAt, double d);

    /*!
       \brief Get number of heavy atom nonbonded neighbors within a certain distance
       \param pAt atom pointer
       \param d distance
       \return number of nonbonded neighbors
    */
    int                      getNumHeavyNeighbors(atom* pAt, double d);

    /*!
       \brief Get number of atoms in the largest residue
       \return number of atoms in the largest residue
    */
    int                      largestResidueSize();

    /*!
       \brief Renumber entire collection
    */
    void                     renumber();

    //! molecule pointer
    molecule*                pMolecule;

    //! submolecule pointer
    submolecule*             pSubMolecule;

    /////////////////////////////
    // -   Metalloproteins   - //
    /////////////////////////////

    /*!
       \brief Determine if the collection contains metal atoms
    */
    bool                     hasMetal();

    /*!
       \brief Finds and stores a vector of all metal atoms
    */
    void                     findMetals();

    /*!
       \brief Determines the coordination sphere of the metal atoms
    */
    void                     determineMetalEnvironments();

    /*!
       \brief Determine if metal centers share ligating residues
    */
    void                     assignMetalPartners();

    /*!
       \brief Assign parameters to metal centers
    */
    void                     assignMetalParameters();

    /*!
       \brief Determines the coordination sphere of the metal atoms
    */
    std::vector<metalCenter*> getMetalCenters();

protected: // functions

    /*!
       \brief Determine matchings for metal centers
       \param pos Current atom being considered
       \param nBondedAtoms Number of atoms
       \param matchMatrix Match matrix
       \param match storage vector
       \param matchings storage vector for all matches
    */
    void findMatchings(int pos, int nBondedAtoms,
                       int matchMatrix[], std::vector<int> &match,
                       std::vector<std::vector<int> > &matchings);

    /*!
       \brief Determine if the fragment is a subgraph of the molecule
       \param iPos Which atom is being considered
       \param jPos Which atom is being considered
       \param nBondedAtoms Number of atoms
       \param matchMatrix Match matrix
    */
    void updateMatchMatrix(int iPos, int jPos, int nBondedAtoms, int matchMatrix[]);

    /*!
       \param nBondedAtoms Number of atoms
       \param matchMatrix Match matrix
       \param mismatch Whether or not a mismatch occurs or not in the updated match match,
       i.e a row that contains all zeros
    */
    void refineMatchMatrix(int nBondedAtoms, int matchMatrix[], bool &mismatch);

protected: // data

    //! molecule list
    std::vector<molecule*>   itsMoleculeList;

    //! collection name
    std::string              itsName;

    //! atom index
    int                      itsAtomIndex;

    //! submolecule index
    int                      itsSubMoleculeIndex;

    //! standard library pointer
    stdLibrary*              pStdLibrary;

    //! parameters pointer
    parameters*              pParameters;

    //! atomTyper pointer
    atomTyper*               pAtomTyper;

    //! Bond parameter pointer
    bondParam*               pBondParam;

    //! Angle parameter pointer
    angleParam*              pAngleParam;

    //! Torsion parameter pointer
    torsionParam*            pTorsionParam;

    //! Improper parameter pointer
    improperParam*           pImproperParam;

    //! H Bond Param
    hBondParam*              pHBondParam;

    //! LJ 6-12 Param
    LJ612Param*              pLJ612Param;

    //! LJ 6-12 Sigma and Epsilon
    LJ612SE*                 pLJ612SE;

    //! Bond pointer
    Bond*                    pBond;

    //! Angle pointer
    Angle*                   pAngle;

    //! Torsion pointer
    Torsion*                 pTorsion;

    //! Improper pointer
    Improper*                pImproper;

    //! Atom Types used
    std::vector<std::string> atomTypesUsed;

    //! Bond Types used
    std::vector<bondParam*> uniqueBondParams;

    //! Angle Types used
    std::vector<angleParam*> uniqueAngleParams;

    //! Torsion Types used
    std::vector<torsionParam*> uniqueTorsionParams;

    //! Improper Types used
    std::vector<improperParam*> uniqueImproperParams;

    //
    // METALS CENTERS
    //

    //! Available metal ions
    std::vector<std::string> availableMetals;

    //! Metal donor bond distances
    std::map<std::string, double> metalDonorDists;

    //! Metal atoms
    std::vector<atom*> itsMetalAtoms;

    //! Metal coordinations
    std::vector<metalCenter*> itsMetalCenters;

    //! Metal groups
    std::vector<metalGroup*> itsMetalGroups;

    //
    // ITERATORS
    //

    //! molecule iterator
    typedef std::vector<molecule*>::iterator moleculeIterator;

    //! submolecule iterator
    typedef std::vector<submolecule*>::iterator  sMolIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator         AtomIterator;

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator       BondMapIterator;

    //! Angle map iterator
    typedef std::map<ULONG_KIND, Angle*>::iterator AngleMapIterator;

    //! Torsion map iterator
    typedef std::map<ULONG_KIND, Torsion*>::iterator    TorsionMapIterator;

    //! Improper map iterator
    typedef std::map<int, Improper*>::iterator   ImproperMapIterator;

    //! Metal Center map iterator
    typedef std::vector<metalCenter*>::iterator MetalCenterIterator;

    //! Metal Group map iterator
    typedef std::vector<metalGroup*>::iterator MetalGroupIterator;

    //! int vector iterator
    typedef std::vector<unsigned int>::iterator intVectorIterator;

    //! string vector iterator
    typedef std::vector<std::string>::iterator stringVectorIterator;

    //! string:double map iterator
    typedef std::map<std::string, double>::iterator strDbMapIterator;

    //! string:int map iterator
    typedef std::map<std::string, int>::iterator strIntMapIterator;

    //! string:string map iterator
    typedef std::map<std::string, std::string>::iterator strStrMapIterator;
};

} // MTKpp namespace

#endif // COLLECTION_H
