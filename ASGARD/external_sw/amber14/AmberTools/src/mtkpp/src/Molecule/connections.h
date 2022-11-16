/*!
   \file connections.h
   \brief Assigns connectivity
   \author Martin Peters

   $Date: 2010/04/29 18:59:17 $
   $Revision: 1.12 $

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

#ifndef CONNECTIONS_H
#define CONNECTIONS_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctype.h>
#include <fstream>
#include <sstream>

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
class metalCenter;
struct Bond;
struct Angle;
struct Torsion;
struct Improper;

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

class protonate;

// ============================================================
// Class : connections()
// ------------------------------------------------------------
/*! 
   \class connections
   \brief Class to assign bonds, angles, torsion and improper
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class connections
{
    friend class protonate;
public:
    /*!
       \brief connections Constructor
       \param parent collection pointer
    */
    connections();

    /*!
       \brief connections Constructor
       \param parent collection pointer
    */
    connections(collection* parent);

    //! connections Destructor
    virtual ~connections();

    /*!
       \brief Assigns all Bonds for Every Molecule in the Collection 
       \todo Figure out how to set bond types, stereo, and topology.
    */
    void assignBonds();

    /*!
       \brief Assigns all Bonds in a Molecule
       \param pMolecule molecule pointer
    */
    void assignBonds(molecule* pMolecule);

    /*!
       \brief Determine all disulfide bonds in the collection
    */
    void assignDisulfideBonds();

    /*!
       \brief Determine all disulfide bonds in a molecule
       \param pMolecule molecule pointer
    */
    void assignDisulfideBonds(molecule* pMolecule);

    /*!
       \brief Assigns all Angles for Every Molecule in the Collection 
    */
    void assignAngles();

    /*!
       \brief Assigns all Angles in a Molecule
       \param pMolecule Molecule pointer
    */
    void assignAngles(molecule* pMolecule);

    /*!
       \brief Assigns all Angles in a metal center
       \param pMetalCenter metalCenter pointer
    */
    void assignAngles(metalCenter* pMetalCenter);

    /*!
       \brief Assigns all Torsions for Every Molecule in the Collection 
    */
    void assignTorsions();

    /*!
       \brief Assigns all Torsions in a Molecule
       \param pMolecule Molecule pointer
    */
    void assignTorsions(molecule* pMolecule);

    /*!
       \brief Assigns all Impropers for Every Molecule in the Collection 
    */
    void assignImpropers();

    /*!
       \brief Assigns all Impropers in a Molecule
       \param pMolecule Molecule pointer
    */
    void assignImpropers(molecule* pMolecule);

    /*!
       \brief Convenience function to carry out assignBonds, assignAngles, assignTorsions, assignImpropers
    */
    void run();

    /*!
       \brief Convenience function to carry out assignBonds, assignAngles, assignTorsions, assignImpropers
    */
    void run(collection* pCol);

    /*!
       \brief Convenience function to carry out assignBonds, assignAngles, assignTorsions, assignImpropers
       \param pMol molecule pointer
    */
    void run(molecule* pMol);

    /*!
       \brief Convenience function to carry out assignStdBonds, assignStdAngles, assignStdTorsions, assignStdImpropers
    */
    void assignStd();

    /*!
       \brief Convenience function to carry out assignStdBonds, assignStdAngles, assignStdTorsions, assignStdImpropers
       \param pMol molecule pointer
    */
    void assignStd(molecule* pMol);

    /*!
       \brief Convenience function to carry out assignStdBonds and assignStdAngles
       \param pMol molecule pointer
    */
    void assignStdBondsAngles(molecule* pMol);

public: // FUNCTIONS

    /*!
       \brief Bond by library
       \param m molecule pointer
       \param s1 submolecule 1 pointer
       \param f1 standard residue 1 
       \param s2 submolecule 2 pointer
       \param f2 standard residue 2
    */
    void bondByLibrary(molecule* m, submolecule* s1, stdFrag* f1, submolecule* s2, stdFrag* f2);

    /*!
       \brief Bond by distance
       \param m molecule pointer 
       \param s1 submoleculs pointer
       \param s2 submolecule pointer
    */
    void bondByDistance(molecule* m, submolecule* s1, submolecule* s2);

    /*!
       \brief Does a bond exist between a and b (used by bondByDistance)
       \param d distance between a and b
       \param a atom pointer
       \param b atom pointer
    */
    virtual bool BondExists(double d, atom* a, atom* b);

    //! Assign Standard Bonds
    void assignStdBonds();

    /*!
       \brief Assign Standard Bonds
       \param pMol molecule pointer
    */
    void assignStdBonds(molecule* pMol);

    //! Assign Standard Angles
    void assignStdAngles();

    /*!
       \brief Assign Standard Angles
       \param pMol molecule pointer
    */
    void assignStdAngles(molecule* pMol);

    //! Assign Standard Torsions
    void assignStdTorsions();

    /*!
       \brief Assign Standard Torsions
       \param pMol molecule pointer
    */
    void assignStdTorsions(molecule* pMol);

    //! Assign Standard Impropers
    void assignStdImpropers();

    /*!
       \brief Assign Standard Impropers
       \param pMol molecule pointer
    */
    void assignStdImpropers(molecule* pMol);

protected: // DATA

    //! collection pointer
    collection*    pCollection;

    //! molecule pointer
    molecule*      pMolecule;

    //! submolecule pointer
    submolecule*   pSubMolecule;

    //! submolecule pointer
    submolecule*   pSubMoleculeMinus1;

    //! atom pointer
    atom*          pAtom;

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

    //! Angle pointer
    Angle*         pAngle;

    //! Torsion pointer
    Torsion*       pTorsion;

    //! Improper pointer
    Improper*      pImproper;

    //! coordinate pointer
    vector3d*      coord1;

    //! coordinate pointer
    vector3d*      coord2;

    //! standard library pointer
    stdLibrary*    pStdLibrary;

    //! standard group pointer
    stdGroup*      pStdGroup;

    //! standard fragment pointer
    stdFrag*       pStdFrag;

    //! standard fragment pointer
    stdFrag*       pStdFragMinus1;

    //! standard atom pointer
    stdAtom*       pStdAtom;

    //! standard atom pointer
    stdAtom*       pStdAtom1;

    //! standard atom pointer
    stdAtom*       pStdAtom2;

    //! standard loop pointer
    stdLoop*       pStdLoop;

    //! molecule vector
    std::vector<molecule*>        moleculeList;

    //! submolecule vector
    std::vector<submolecule*>     subMoleculeList;

    //! atom vector
    std::vector<atom*>            atomList;

    //! atom vector
    std::vector<atom*>            atomList2;

    //! standard loop vector
    std::vector<stdLoop*>         stdLoopList;

    //! parameters pointer
    parameters*    pParameters;

    //! Missing bond parameters
    int missingBondParams;

    //! Missing angle parameters
    int missingAngleParams;

    //! Missing torsion parameters
    int missingTorsionParams;

    //! Error message
    std::string errorMessage;

};

} // MTKpp namespace

#endif // CONNECTIONS_H

