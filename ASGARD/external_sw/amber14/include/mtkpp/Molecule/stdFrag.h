/*!
   \file stdFrag.h
   \brief Container for standard atoms, bonds, angles, torsions and impropers
   \author Martin Peters

   $Date: 2010/04/29 18:59:17 $
   $Revision: 1.17 $

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

#ifndef STDFRAG_h
#define STDFRAG_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdio.h>

#include "Utils/constants.h"

namespace MTKpp
{

class stdGroup;
class vector3d;

/*!
   \struct stdAtom
   \brief Container for standard atom info
*/
struct stdAtom {
    //! atom name
    std::string    identity;

    //! internal fragment index
    int            index;

    //! atom type
    std::string    type;

    /*!
       \brief chain type
       - M => Main chain
       - S => Side, 2 connections
       - E => End, 1 connection
       - B => Branch, 3 connections
       - 3 => 4 connections
       - 4 => 5 connections
       - 5 => 6 connections
       - 6 => 7 connections
    */
    std::string    chain;

    //! atomic charge for MM calculation
    double         atmCharge;

    //! index of atom with which its bonded to
    int            bond12;

    //! bond length
    double         bondLength;

    //! index of atom with which its forming an angle with bond12
    int            bond13;

    //! angle size
    double         bondAngle;

    //! index of atom with which its forming a torsion with bond12 and angle13
    int            bond14;

    //! torsion size
    double         bondTorsion;

    //! atomic number
    int            atNum;

    //! atomic symbol
    std::string    atSymbol;

    //std::string    hybridization;

    /*!
       \brief atom kind definitions
        - 0 = Undefined
        - 1 = Hydrogen
        - 2 = Terminal Heavy Atom
        - 3 = Open Chain Heavy Atom
        - 4 = Closed Chain Heavy Atom
        - 5 = Ring Heavy Atom
        - 6 = Aromatic Ring Heavy Atom
        - 7 = Chain Atom (not terminal, Ring or Hydrogen)
    */
    int            kind;
};

/*!
   \struct stdBond
   \brief Container for standard bond info
*/
struct stdBond {
    //! atom 1 index
    int            atom1;

    //! atom 2 index
    int            atom2;

    /*!
       \brief bond type definitions
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
    int            type;

    /*!
        \brief bond stereo
        - Definitions for Single Bonds
         -# 0 = Not stereo
         -# 1 = Up
         -# 4 = Either
         -# 6 = Down
        - Definitions for Double Bonds
         -# 0 = Use x,y,z coords from atom block to determine cis or trans
         -# 3 = Either cis or trans
    */
    int            stereo;

    /*!
        \brief bond topology definitions
        - 0 = Either
        - 1 = Ring
        - 2 = Chain
    */
    int            topology;

    /*!
        Bond Kind Definitions
        - 0 = Undefined
        - 1 = Polar
    */
    int            kind;

    //! bond length
    double         length;
};

/*!
   \struct stdImproper
   \brief Container for standard improper info
*/
struct stdImproper {
    //! atom 1 index
    int            atom1;

    //! atom 2 index
    int            atom2;

    //! atom 3 index
    int            atom3;

    //! atom 4 index
    int            atom4;
};

/*!
   \struct stdLoop
   \brief Container for standard loop info
*/
struct stdLoop {
    //! atom 1 index
    int            atom1;

    //! atom 2 index
    int            atom2;

    //! see stdBond for details
    int            type;

    //! see stdBond for details
    int            stereo;
};

/*!
   \struct stdAlias
   \brief Container for standard alias info
*/
struct stdAlias {
    //! atom 1 name
    std::string    atom1;

    //! atom 2 name
    std::string    atom2;
};

/*!
   \struct stdRing
   \brief Container for standard ring info
*/
struct stdRing {
    //! atom indices in ring
    std::vector<int>    atoms;

    //! Size of ring
    int                 size;

    /*!
       Is the ring planar?
       - 0 = No
       - 1 = Yes
    */
    int                 planar;

    /*!
       Is the ring aromatic?
       - 0 = No
       - 1 = Yes
    */
    int                 aromatic;

    /*!
       Is the ring a heterocycle?
       - 0 = No
       - 1 = yes
    */
    int                 hetero;

    //! Number of heteroatoms
    int                 nHetero;

    //! Number of Nitrogen atoms
    int                 nNitrogen;

    //! Number of Oxygen atoms
    int                 nOxygen;

    //! Number of Sulfur atoms
    int                 nSulfur;
};

/*!
   \struct stdFeature
   \brief Container for standard feature info
*/
struct stdFeature {
    /*!
       \brief name of feature
        - HBA Hydrogen Bond Acceptor
        - HBD Hydrogen Bond Donor
        - HPB Hydrophobic
        - PCC Positive Charge Center
        - NCC Negative Charge Center
        - PIC Aromatic Center
    */
    std::string         name;

    /*!
       \brief atom indices in feature
    */
    std::vector<int>    atoms;
};

/*!
   \struct stdFuncGroup
   \brief Container for standard functional group

    This functional group container is only used for molecules
*/
struct stdFuncGroup {
    /*!
       \brief name of fragment group, e.g. terminal
    */
    std::string         groupName;

    /*!
       \brief name of fragment, e.g. CH3
    */
    std::string         fragName;

    /*!
       \brief atom indices in functional group
    */
    std::vector<int>    atoms;
};

/*!
   \struct stdConnTorsion
   \brief Container for standard connection torsions
*/
struct stdConnTorsion {
    /*!
       \brief connection atom
    */
    int                 bondAtom;

    /*!
       \brief angle atom
    */
    int                 angleAtom;

    /*!
       \brief torsion atom
    */
    int                 torsionAtom;

    /*!
       \brief torsion value
    */
    double              torsion;
};

/*!
   \struct stdRotBond
   \brief Container for standard rotatable bonds
*/
struct stdRotBond {
    /*!
       \brief first atom
    */
    int                 atom1;

    /*!
       \brief second atom
    */
    int                 atom2;

    /*!
       \brief third atom
    */
    int                 atom3;

    /*!
       \brief fourth atom
    */
    int                 atom4;

    /*!
       \brief torsion values available to the rotatable bond
    */
    std::vector<double>    values;
};

// ============================================================
// Class : stdFrag()
// ------------------------------------------------------------
/*!
   \class stdFrag

   \brief Container for stdAtoms, stdBond, stdImpropers, stdLoops,
          stdAlias'
*/
// ============================================================
class stdFrag
{
public:
    /*!
       \brief stdFrag Constructor
       \param parent stdGroup pointer
    */
    stdFrag(stdGroup *parent = 0);

    /*!
       \brief stdFrag Copy Constructor
       \param parent stdGroup pointer
       \param sf stdFrag pointer
    */
    stdFrag(stdFrag *sf, stdGroup *parent = 0);

    /*!
       \brief Standard Fragment Destructor
    */
    virtual ~stdFrag();

    /*!
       \brief Get group which fragment is apart of
       \return stdGroup pointer
    */
    stdGroup*                getParent();

    /*!
       \brief Set name of fragment
       \param name Name of fragment
    */
    void                     setName(const std::string& name);

    /*!
       \brief Get name of the fragment
       \return fragment name
    */
    std::string              getName();

    /*!
       \brief Set symbol of the fragment
       \param symbol Fragment symbol
    */
    void                     setSymbol(const std::string& symbol);

    /*!
       \brief Get symbol of the fragment
       \return fragment symbol
    */
    std::string              getSymbol();

    /*!
       \brief Set 8 character code of the fragment
       \param code Fragment symbol
    */
    void                     setCode(const std::string& code);

    /*!
       \brief Get code of the fragment
       \return fragment 8 character code
    */
    std::string              getCode();

    /*!
       \brief Set 1 letter character of the fragment
       \param c character
    */
    void                     setCharacter(const std::string& c);

    /*!
       \brief Get 1 letter character of the fragment
       \return character
    */
    std::string              getCharacter();

    /*!
       \brief Set type of the fragment
       \param type Fragment type
    */
    void                     setType(const std::string& type);

    /*!
       \brief Get type of the fragment
       \return fragment type
    */
    std::string              getType();

    /*!
       \brief Set symmetry of the fragment
       \param sym Fragment symmetry
    */
    void                     setSymmetry(const std::string& sym);

    /*!
       \brief Get symmetry of the fragment
       \return fragment symmetry
    */
    std::string              getSymmetry();

    /*!
       \brief Set which fragments this fragment is a subgraph of
       \param g
    */
    void                     setSubGraphs(const std::vector<std::string>& g);

    /*!
       \brief Get list of fragments related to this one
       \return list of fragments related to this one
    */
    std::vector<std::string> getSubGraphs();

    /*!
       \brief Get list of fragments related to this one
       \return list of fragments related to this one
    */
    std::string              getSubGraphStr();

    /*!
       \brief Add a stdAtom to the fragment
       \return stdAtom
    */
    stdAtom*                 addStdAtom();

    /*!
       \brief Add a stdAtom to the fragment
       \param s stdAtom pointer
       \return stdAtom
    */
    stdAtom*                 addStdAtom(stdAtom* s);

    /*!
       \brief Get a stdAtom by atom index
       \param index atom index
       \return stdAtom
    */
    stdAtom*                 getStdAtom(const int& index);

    /*!
       \brief Get a stdAtom by name
       \param name atom name
       \return stdAtom
    */
    stdAtom*                 getStdAtom(const std::string& name);

    /*!
       \brief Has a stdAtom by name
       \param name atom name
       \return stdAtom
    */
    bool                     hasStdAtom(const std::string& name);

    /*!
       \brief Get a list of stdAtoms which are bonded to sAt
       \param sAt stdAtom pointer
       \return list of stdAtoms bonded to sAt
    */
    std::vector<stdAtom*>    getBondedStdAtoms(stdAtom* sAt);

    /*!
       \brief Get number of stdAtoms in fragment
       \return number of stdAtoms
    */
    int                      numStdAtoms();

    /*!
       \brief Get number of stdAtoms in fragment
       \return number of stdAtoms
    */
    int                      numStdHeavyAtoms();

    /*!
       \brief Add a stdBond to the fragment
       \return stdBond
    */
    stdBond*                 addStdBond();

    /*!
       \brief Add a stdBond to the fragment
       \param b stdBond pointer
       \return stdBond
    */
    stdBond*                 addStdBond(stdBond* b);

    /*!
       \brief get a stdBond from the fragment based on atom indices
       \param at1 atom 1 index
       \param at2 atom 2 index
       \return stdBond
    */
    stdBond*                 getStdBond(const int& at1, const int& at2);

    /*!
       \brief get a stdBond from the fragment based on atom indices
       \param at1 atom 1 name
       \param at2 atom 2 name
       \return stdBond
    */
    stdBond*                 getStdBond(const std::string& at1, const std::string& at2);

    /*!
       \brief get a stdBond from the fragment based on atom indices
       \param pAt1 stdAtom 1 pointer
       \param pAt2 stdAtom 2 pointer
       \return stdBond
    */
    stdBond*                 getStdBond(stdAtom* pAt1, stdAtom* pAt2);

    /*!
       \brief Does the fragment contain a bond
       \param at1 atom 1 index
       \param at2 atom 2 index
       \return true/false
    */
    bool                     hasStdBond(const int& at1, const int& at2);

    /*!
       \brief Get the number of stdBonds in the fragment
       \return number of bonds
    */
    int                      numStdBonds();

    /*!
       \brief Get the number of stdBonds which pStdAtom is apart of
       \param pStdAtom stdAtom pointer
       \return number of bonds
    */
    int                      numStdBonds(stdAtom* pStdAtom);

    /*!
       \brief Add a stdImproper to the fragment
       \return stdImproper
    */
    stdImproper*             addStdImproper();

    /*!
       \brief Add a stdImproper to the fragment
       \param i stdImproper pointer
       \return stdImproper
    */
    stdImproper*             addStdImproper(stdImproper* i);

    /*!
       \brief Add a stdLoop to the fragment
       \return stdLoop
    */
    stdLoop*                 addStdLoop();

    /*!
       \brief Add a stdLoop to the fragment
       \param l stdLoop pointer
       \return stdLoop
    */
    stdLoop*                 addStdLoop(stdLoop* l);

    /*!
       \brief get a stdLoop from the fragment based on atom indices
       \param at1 atom 1 index
       \param at2 atom 2 index
       \return stdLoop
    */
    stdLoop*                 getStdLoop(const int& at1, const int& at2);

    /*!
       \brief get a stdLoop from the fragment based on atom indices
       \param pAt1 stdAtom 1 pointer
       \param pAt2 stdAtom 2 pointer
       \return stdLoop
    */
    stdLoop*                 getStdLoop(stdAtom* pAt1, stdAtom* pAt2);

    /*!
       \brief Add a stdAlias to the fragment
       \return stdAlias
    */
    stdAlias*                addStdAlias();

    /*!
       \brief Add a stdAlias to the fragment
       \param a stdAlias pointer
       \return stdAlias
    */
    stdAlias*                addStdAlias(stdAlias* a);

    /*!
       \brief Add a stdRing to the fragment
       \return stdRing
    */
    stdRing*                 addStdRing();

    /*!
       \brief Add a stdRing to the fragment
       \param r stdRing pointer
       \return stdRing
    */
    stdRing*                 addStdRing(stdRing* r);

    /*!
       \brief Add a stdFeature to the fragment
       \return stdFeature
    */
    stdFeature*              addStdFeature();

    /*!
       \brief Add a stdFeature to the fragment
       \param f stdFeature pointer
       \return stdFeature
    */
    stdFeature*              addStdFeature(stdFeature* f);

    /*!
       \brief Has the stdAtom got a stdFeature matching f
       \param a stdAtom pointer
       \param f stdFeature string
       \return bool
    */
    bool                     hasStdFeature(stdAtom* a, std::string f);

    /*!
       \brief Add a stdFuncGroup to the stdMolecule
       \return stdFeature
    */
    stdFuncGroup*            addStdFuncGroup();

    /*!
       \brief Add a stdFuncGroup to the stdMolecule
       \param f stdFuncGroup pointer
       \return stdFeature
    */
    stdFuncGroup*            addStdFuncGroup(stdFuncGroup* f);

    /*!
       \brief Add connection points to the fragment
       \param v connection points vector
    */
    void                     addStdConnPts(std::vector<int> &v);

    /*!
       \brief Add a stdConnTorsion to the stdMolecule
       \return stdConnTorsion
    */
    stdConnTorsion*          addStdConnTorsion();

    /*!
       \brief Add a stdRotBond to the stdMolecule
       \return stdRotBond pointer
    */
    stdRotBond*              addStdRotBond();

    /*!
       \brief Get all stdAtoms in the fragment
       \return vector of stdAtoms
    */
    std::vector<stdAtom*>    getStdAtomList();

    /*!
       \brief Get all stdBonds in the fragment
       \return vector of stdBonds
    */
    std::vector<stdBond*>    getStdBondList();

    /*!
       \brief Get all stdLoops in the fragment
       \return vector of stdLoops
    */
    std::vector<stdLoop*>    getStdLoopList();

    /*!
       \brief Get all stdAlias' in the fragment
       \return vector of stdAlias
    */
    std::vector<stdAlias*>   getStdAliasList();

    /*!
       \brief Get all stdImpropers in the fragment
       \return vector of stdImpropers
    */
    std::vector<stdImproper*>getStdImproperList();

    /*!
       \brief Get all stdRings in the fragment
       \return vector of stdRings
    */
    std::vector<stdRing*>    getStdRingList();

    /*!
       \brief Get all stdFeatures in the fragment
       \return vector of stdFeatures
    */
    std::vector<stdFeature*> getStdFeatureList();

    /*!
       \brief Get all stdFuncGroups in the fragment
       \return vector of stdFuncGroups
    */
    std::vector<stdFuncGroup*> getStdFuncGroupList();

    /*!
       \brief Get all stdConnPts in the fragment
       \return vector of stdConnPts (atom indices)
    */
    std::vector<int>         getStdConnPtsList();

    /*!
       \brief Get number of stdConnPts in the fragment
       \return number of stdConnPts
    */
    int                      getNumStdConnPtsList();

    /*!
       \brief Get all stdConnTorsion in the fragment
       \return vector of stdConnTorsion
    */
    std::vector<stdConnTorsion*> getStdConnTorList();

    /*!
       \brief Get all stdRotBonds in the fragment
       \return vector of stdRotBond
    */
    std::vector<stdRotBond*> getStdRotBondList();

    /*!
       \brief Get stdAlias by name
       \param name atom name
       \return stdAtom alias
    */
    std::string              getAlias(const std::string& name);

    /*!
       \brief Generate fragment coordinates
       \return success
    */
    int                      generateCoordinates();

    /*!
       \brief Generate fragment coordinates using known coordinates
       \param bd atom to which the first standard atom is bonded to
       \param ag atom to which the first standard atom forms an angle with
       \param tr atom to which the first standard atom forms a torsion with
       \param forward proceed through the fragment in a forward/reverse manner
       \return success
    */
    int                      generateCoordinates(vector3d* bd, vector3d* ag,
                             vector3d* tr, const int &forward = 1);

    /*!
       \brief Get fragment coordinates
    */
    std::vector<vector3d*>   getCoordinates();

    /*!
       \brief Print details of stdFrag
    */
    void                     print();

    /*!
       \brief Generate simple fingerprint
    */
    void                     generateSimpleFP();

    /*!
       \brief Get simple fingerprint
       \return simple fingerprint
    */
    std::vector<unsigned int>getSimpleFP();

    /*!
       \brief Generate adjacency matrix
    */
    int                      generateAdjMatrix();

    /*!
       \brief Generate heavy atom adjacency matrix
    */
    int                      generateHeavyAdjMatrix();

    /*!
       \brief Get adjacency matrix
       \return Adjacency matrix
    */
    int*                     getAdjMatrix();

    /*!
       \brief Get heavy atom adjacency matrix
       \return heavy atom Adjacency matrix
    */
    int*                     getHeavyAdjMatrix();

    /*!
       \brief Get adjacency matrix size
       \return Adjacency matrix size
    */
    int                      getAdjMatrixSize();

    /*!
       \brief Get heavy atom adjacency matrix size
       \return heavy atom adjacency matrix size
    */
    int                      getHeavyAdjMatrixSize();

    /*!
       \brief Generate atom symbols array
       \return success
    */
    int                      generateAtomSymbols();

    /*!
       \brief Generate heavy atom symbols array
       \return success
    */
    int                      generateHeavyAtomSymbols();

    /*!
       \brief Get atom symbols array
       \return atom symbols array
    */
    char*                    getAtomSymbols();

    /*!
       \brief Get heavy atom symbols array
       \return heavy atom symbols array
    */
    char*                    getHeavyAtomSymbols();

    /*!
       \brief Generate atom kinds array
       \return success
    */
    int                      generateAtomKinds();

    /*!
       \brief Generate heavy atom kinds array
       \return success
    */
    int                      generateHeavyAtomKinds();

    /*!
       \brief Get atom kinds array
       \return atom kinds array
    */
    int*                     getAtomKinds();

    /*!
       \brief Get heavy atom kinds array
       \return atom kinds array
    */
    int*                     getHeavyAtomKinds();

    /*!
       \brief Get charge
       \return charge
    */
    double                   getCharge();

    /*!
       \brief Get index of stdAtom
       \param a stdAtom pointer
       \return index of standard atom
    */
    int                      getStdAtomIndex(stdAtom* a);

protected:

    //! standard fragment name
    std::string              itsName;

    //! standard fragment symbol (3 Character)
    std::string              itsSymbol;

    //! standard fragment code (8 Character)
    std::string              itsCode;

    //! 1 Letter character
    std::string              itsCharacter;

    /*!
        fragment type Definitions
        - m   Independent molecule
        - s   terminal (start) chain group
        - e   terminal (end) chain group
        - t   terminal chain (start/ end) group
        - l   chain group
    */
    std::string              itsType;

    /*!
       fragment symmetry Definitions
       - C2   C2 symmetric (180 degree rotation results in equivalent geometry)
    */
    std::string              itsSymmetry;

    //! list of 3L code fragments which this fragment is a subgraph of
    std::vector<std::string> itsSubGraphs;

    //! standard atom iterator
    typedef std::vector<stdAtom*>::iterator stdAtomIterator;

    //! standard bond iterator
    typedef std::vector<stdBond*>::iterator stdBondIterator;

    //! standard loop iterator
    typedef std::vector<stdLoop*>::iterator stdLoopIterator;

    //! standard alias iterator
    typedef std::vector<stdAlias*>::iterator stdAliasIterator;

    //! standard improper iterator
    typedef std::vector<stdImproper*>::iterator stdImproperIterator;

    //! standard ring iterator
    typedef std::vector<stdRing*>::iterator stdRingIterator;

    //! standard feature iterator
    typedef std::vector<stdFeature*>::iterator stdFeatureIterator;

    //! standard functional group iterator
    typedef std::vector<stdFuncGroup*>::iterator stdFuncGroupIterator;

    //! standard connection points iterator
    typedef std::vector<int>::iterator stdConnPtsIterator;

    //! standard connTorsion iterator
    typedef std::vector<stdConnTorsion*>::iterator stdConnTorsionIterator;

    //! standard rotatable bond iterator
    typedef std::vector<stdRotBond*>::iterator stdRotBondIterator;

    //! standard atom list
    std::vector<stdAtom*>    itsStdAtomList;

    //! standard bond list
    std::vector<stdBond*>    itsStdBondList;

    //! standard loop list
    std::vector<stdLoop*>    itsStdLoopList;

    //! standard alias list
    std::vector<stdAlias*>   itsStdAliasList;

    //! standard improper list
    std::vector<stdImproper*>itsStdImproperList;

    //! standard ring list
    std::vector<stdRing*>    itsStdRingList;

    //! standard feature list
    std::vector<stdFeature*> itsStdFeatureList;

    //! standard functional group list
    std::vector<stdFuncGroup*> itsStdFuncGroupList;

    //! standard connection points list
    std::vector<int>         itsStdConnPtsList;

    //! stdConnTorsion list
    std::vector<stdConnTorsion*> itsStdConnTorsionList;

    //! stdRotBond list
    std::vector<stdRotBond*> itsStdRotBondList;

    //! standard improper list
    std::vector<vector3d*>   itsCoords;

    //! Simple fingerprint
    std::vector<unsigned int>itsSimpleFP;

    //! Fragment fingerprint
    std::vector<unsigned int>itsFragmentFP;

    //------------//
    //- POINTERS -//
    //------------//

    //! standard group which fragment belongs to
    stdGroup*                pParent;

    //! standard atom pointer
    stdAtom*                 pStdAtom;

    //! standard atom pointer
    stdAtom*                 pStdAtom1;

    //! standard atom pointer
    stdAtom*                 pStdAtom2;

    //! vector3d pointer
    vector3d*                pCoords;

    //! standard bond pointer
    stdBond*                 pStdBond;

    //! standard loop pointer
    stdLoop*                 pStdLoop;

    //! standard alias pointer
    stdAlias*                pStdAlias;

    //! standard improper pointer
    stdImproper*             pStdImproper;

    //! standard ring pointer
    stdRing*                 pStdRing;

    //! standard feature pointer
    stdFeature*              pStdFeature;

    //! standard function group pointer
    stdFuncGroup*            pStdFuncGroup;

    //! standard connection torion pointer
    stdConnTorsion*          pStdConnTorsion;

    //! standard rotatable bond pointer
    stdRotBond*              pStdRotBond;

    //! Adjacency Matrix
    int                      *adjMatrix;

    //! Adjacency Matrix size
    int                      adjMatrixSize;

    //! Atom symbols
    char                     *atomSymbols;

    //! Atom kinds
    int                      *atomKinds;

    //! Adjacency Matrix
    int                      *heavyAtomAdjMatrix;

    //! Adjacency Matrix size
    int                      heavyAtomAdjMatrixSize;

    //! Atom symbols
    char                     *heavyAtomSymbols;

    //! Atom kinds
    int                      *heavyAtomKinds;
};

} // MTKpp namespace

#endif // STDFRAG_H

