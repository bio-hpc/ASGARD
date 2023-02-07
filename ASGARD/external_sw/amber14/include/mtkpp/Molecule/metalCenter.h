/*!
   \file metalCenter.h
   \brief Container for metaloprotein information
   \author Martin Peters

   Container for metalloprotein information

   $Date: 2010/04/29 18:59:17 $
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

#ifndef METALCENTER_H
#define METALCENTER_H

#include "Utils/constants.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

namespace MTKpp
{
class collection;
class molecule;
class submolecule;
class atom;
class Bond;
class Angle;
class Torsion;
class vector3d;
class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;
struct bondParam;
struct angleParam;
struct torsionParam;

// ============================================================
// Class : metalCenter()
// ------------------------------------------------------------
/*!
   \class metalCenter
   \brief Container for metal centers
*/
// ============================================================
class metalCenter
{
public:

    /*!
       \brief Constructor
    */
    metalCenter(atom* metal, std::map<std::string, double> m_d);

    /*!
       \brief Destructor
    */
    virtual ~metalCenter();

    /*!
       \brief Set name
    */
    void setName(std::string n) {itsName = n;}

    /*!
       \brief Get name
    */
    std::string getName() {return itsName;}

    /*!
       \brief Get metal, this may want to be a vector in the future
    */
    atom* getMetalAtom();

    /*!
       \brief pdb file where metal atom comes from
    */
    void setPdb(std::string f);

    /*!
       \brief Returns the value of error
    */
    void setError(int e);

    /*!
       \brief Returns the value of error
    */
    int getError();

    /*!
       \brief add bonded atom
    */
    void addAtom(atom* a);

    /*!
       \brief add bonded atom label
    */
    void addLabel(std::string s);

    /*!
       \brief Get Primary shell
    */
    std::string getPrimaryShell();

    /*!
       \brief Get Primary shell atoms
       \param vSA vector of primary shell atoms
    */
    void getPrimaryShellAtoms(std::vector<atom*>& vSA);

    /*!
       \brief Get Secondary shell
    */
    std::string getSecondaryShell();

    /*!
       \brief Get Secondary shell atoms
       \param vSA vector of secondary shell atoms
    */
    void getSecondaryShellAtoms(std::vector<atom*>& vSA);

    /*!
       \brief Get most likely geometry
    */
    std::string getGeometry();

    /*!
       \brief Get most likely geometry rmsd
    */
    double getGeometryRMS();

    /*!
       \brief Check if a residue is bidentate
    */
    void checkBidentate(int at1, int at2);

    /*!
       \brief Assign bond type between metal center and bonded atoms
    */
    void assignBondType();

    /*!
       \brief Write metal cluster pdb file
    */
    void writeEnvironment(std::string file_name);

    /*!
       \brief Assign coordination state
    */
    void assignGeometry();

    /*!
       \brief Print metal center info
    */
    void print(std::ostream& os);

    /*!
       \brief Return metal center info
    */
    std::string getInfo();

    /*!
       \brief Test if 4/5 coordinate is tet, sqp, or tnb
    */
    void test4Coord(std::vector<atom*> atoms);

    /*!
       \brief Determine the tbp_rms and ttp_rms
    */
    void test5Coord(std::vector<atom*> atoms);

    /*!
       \brief Determine the oct_rms
    */
    void test6Coord(std::vector<atom*> atoms);

    /*!
       \brief Set stdGroup pointer
       \param pStdG stdGroup pointer
    */
    void setStdGroup(stdGroup* pStdG);

    /*!
       \brief Get stdGroup pointer
       \return stdGroup pointer
    */
    stdGroup* getStdGroup();

    /*!
       \brief Add Bond to metal centers
       \param at1 atom pointer
       \param at2 atom pointer
       \param type integer value designating single, double, triple etc., bond type
       \param stereo integer value designating bond stereo
       \param topology integer value assigning ring membership or not
       \param size distance between at1 and at2
       \return Bond pointer
    */
    Bond*                    addBond(atom* at1, atom* at2,
                                     const int& type = 1,
                                     const int& stereo = 0,
                                     const int& topology = 0,
                                     const double& size = 0.0);

    /*!
       \brief Does metal center contain a certain bond
       \param at1 atom pointer
       \param at2 atom pointer
       \return boolean
    */
    bool                     hasBond(atom* at1, atom* at2);

    /*!
       \brief Get Bond map
       \return bond map <int Bond*>
    */
    std::map<int, Bond*>     getBondMap();

    /*!
       \brief Get number of bonds
       \return number of bonds
    */
    int                      numBonds();

    /*!
       \brief Add Angle to metal center
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param ang Angle between at1, at2, and at3
       \return Angle pointer
    */
    Angle*                   addAngle(atom* at1, atom* at2, atom* at3,
                                      const double& ang = 0.0);

    /*!
       \brief Does metal center contain a certain angle
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \return boolean
    */
    bool                     hasAngle(atom* at1,atom* at2,atom* at3);

    /*!
       \brief Get Angle map
       \return angle map <ULONG_KING Angle*>
    */
    std::map<ULONG_KIND, Angle*>    getAngleMap();

    /*!
       \brief Get number of angles
       \return number of angles
    */
    int                      numAngles();

    /*!
       \brief Add Torsion to metal center
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \param tor Torsion between at1, at2, at3, and at4
       \return Torsion pointer
    */
    Torsion*                 addTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4,
                                        const double& tor);

    /*!
       \brief Does metal center contain a certain torsion
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \return boolean
    */
    bool                     hasTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4);

    /*!
       \brief Get Torsion map
       \return torsion map <ULONG_KIND Torsion*>
    */
    std::map<ULONG_KIND, Torsion*>     getTorsionMap();

    /*!
       \brief Get number of torsions
       \return number of torsions
    */
    int                                numTorsions();

protected:

    //! metal center name
    std::string              itsName;

    //! metal atom
    atom*                    metalAtom;

    //! molecule pointer
    molecule*                pMolecule;

    //! collection pointer
    collection*              pCollection;

    //! stdLibrary pointer
    stdLibrary*              pStdLibrary;

    //! metal atom coordinations
    vector3d*                metalCoords;

    //! metal atom's bonded atoms
    std::vector<atom*>       itsBondedAtoms;

    //! metal:residue label
    std::vector<std::string> itsLabels;

    //! bonded atoms residue names
    std::vector<std::string> itsResNames;

    //! bonded atoms residue numbers
    std::vector<int>         itsResNums;

    //! metal:residue label = double map
    std::map<std::string, double> metal_donors;

    //! bonded atom coordination types
    std::vector<std::string> coordinationType;

    //! bonded atoms bools
    std::vector<bool>        assigned;

    //! storage for primary shell info
    std::string              primaryShell;

    //! storage for secondary shell info
    std::string              secondaryShell;

    //! pdb file where metal atom comes from
    std::string              pdbFile;

    /*!
       \brief Most likely geometry
       - tetrahedral (tet) 109.5
       - tetrahedral with a weak bonded atom (tnb)
       - square planar (sqp) 127.3
       - trigonal bipyramidal (tbp) 111.4
       - tetragonal pyramid (ttp)
       - octahedral (oct)
    */
    std::string              geometry1;

    /*!
       \brief Second most likely geometry
       - tetrahedral (tet) 109.5
       - tetrahedral with a weak bonded atom (tnb)
       - square planar (sqp) 127.3
       - trigonal bipyramidal (tbp) 111.4
       - tetragonal pyramid (ttp)
       - octahedral (oct)
    */
    std::string              geometry2;

    /*!
       RMSD for geometry1 from ideal geometry
    */
    double                   geometry1RMS;

    /*!
       RMSD for geometry2 from ideal geometry
    */
    double                   geometry2RMS;

    //! Error
    int                      error;

    //! stdGroup pointer
    stdGroup*                pStdGroup;

    //! string:double map iterator
    typedef std::map<std::string, double>::iterator strDbMapIterator;

    //! string:int map iterator
    typedef std::map<std::string, int>::iterator strIntMapIterator;

    //! string:string map iterator
    typedef std::map<std::string, std::string>::iterator strStrMapIterator;

    //! int vector iterator
    typedef std::vector<unsigned int>::iterator intVectorIterator;

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator BondMapIterator;

    //! Angle map iterator
    typedef std::map<ULONG_KIND, Angle*>::iterator AngleMapIterator;

    //! Torsion map iterator
    typedef std::map<ULONG_KIND, Torsion*>::iterator TorsionMapIterator;

    //! Bond pointer
    Bond*                    pBond;

    //! Angle pointer
    Angle*                   pAngle;

    //! Torsion pointer
    Torsion*                 pTorsion;

    //! Bond parameter pointer
    bondParam*               pBondParam;

    //! Angle parameter pointer
    angleParam*              pAngleParam;

    //! Torsion parameter pointer
    torsionParam*            pTorsionParam;

    //! Bond map
    std::map<int, Bond*>     itsBondMap;

    //! Angle map
    std::map<ULONG_KIND, Angle*>    itsAngleMap;

    //! Torsion map
    std::map<ULONG_KIND, Torsion*> itsTorsionMap;
};

// ============================================================
// Class : metalGroup()
// ------------------------------------------------------------
/*! 
   \class metalGroup
   \brief Container for metal centers
   \author Martin Peters
*/
// ============================================================
class metalGroup
{
public:

    /*!
       \brief Constructor
    */
    metalGroup();

    /*!
       \brief Destructor
    */
    virtual ~metalGroup();

    /*!
       \brief Add metal center
    */
    void addMetalCenter(metalCenter* m);

    /*!
       \brief Add metal center
    */
    std::vector<metalCenter*> getMetalCenters();

private:
    //! Metal centers
    std::vector<metalCenter*> metalCenters;
};

} // MTKpp namespace

#endif // METALCENTER_H
