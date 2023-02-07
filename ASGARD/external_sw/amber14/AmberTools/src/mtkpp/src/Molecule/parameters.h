/*!
   \file parameters.h
   \brief Container for parameter information
   \author Martin Peters

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.14 $

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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include "Utils/constants.h"

namespace MTKpp
{

class elements;

// ============================================================
// Struct : atomType
// ------------------------------------------------------------
// 
// ============================================================
struct atomType
{
    //! atom type name
    std::string name;

    //! element symbol
    std::string element;

    //! atomic number
    int atNum;

    //! mass of atom type
    double mass;

    //! hybridization
    std::string hybridization;

    //! description of atom type
    std::string description;

    //! van der Waals radius
    double rvalue;

    //! Well Depth
    double evalue;

    //! Atomic polarizability in A**3
    double atomPolarizability;

    //! group name
    std::string groupName;

    //! Optimize
    bool optimize;
};

// ============================================================
// Struct : bondParam
// ------------------------------------------------------------
// 
// ============================================================
struct bondParam
{
    //! Atom Type of atom 1
    std::string atomType1;

    //! Atom Type of atom 2
    std::string atomType2;

    //! Bond Force Constant
    double keq;

    //! Equilibrium Bond Length
    double req;

    //! group name
    std::string groupName;

    //! Optimize
    bool optimize;
};

// ============================================================
// Struct : angleParam
// ------------------------------------------------------------
// 
// ============================================================
struct angleParam
{
    //! Atom Type of atom 1
    std::string atomType1;

    //! Atom Type of atom 2
    std::string atomType2;

    //! Atom Type of atom 3
    std::string atomType3;

    //! Angle Force Constant
    double keq;

    //! Equilibrium Angle Size
    double req;

    //! group name
    std::string groupName;

    //! Optimize
    bool optimize;
};

// ============================================================
// Struct : torsionParam
// ------------------------------------------------------------
// 
// ============================================================
struct torsionParam
{
    //! Atom Type of atom 1
    std::string atomType1;

    //! Atom Type of atom 2
    std::string atomType2;

    //! Atom Type of atom 3
    std::string atomType3;

    //! Atom Type of atom 4
    std::string atomType4;

    //! The periodicity of the torsion
    double Nt;

    //! Magnitude of torsion in kcal/mol
    double Vn;

    //! Phase offset in deg
    double gamma;

    //! No. of bond paths
    int npth;

    //! group name
    std::string groupName;
}; 

// ============================================================
// Struct : improperParam
// ------------------------------------------------------------
// 
// ============================================================
struct improperParam
{
    //! Atom Type of atom 1
    std::string atomType1;

    //! Atom Type of atom 2
    std::string atomType2;

    //! Atom Type of atom 3
    std::string atomType3;

    //! Atom Type of atom 4
    std::string atomType4;

    //!
    double Nt;

    //!
    double Vn;

    //!
    double gamma;

    //! group name
    std::string groupName;
};

// ============================================================
// Struct : hBondParam
// ------------------------------------------------------------
// 
// ============================================================
struct hBondParam
{
    //! Atom Type of atom 1
    std::string atomType1;

    //! Atom Type of atom 2
    std::string atomType2;

    //!
    double p10;

    //!
    double p12;

    //! group name
    std::string groupName;
};

// ============================================================
// Struct : equivalentAtoms
// ------------------------------------------------------------
// 
// ============================================================
struct equivalentAtomsParam
{
    //! Atom Type of original
    std::string original;

    //! List of equivalent Atom Types
    std::vector<std::string> itsEquivalentList;

    //! group name
    std::string groupName;
};

// ============================================================
// Struct : LJ612SE
// ------------------------------------------------------------
// 
// ============================================================
struct LJ612SE
{
    //! Atom Type of atom 1
    std::string atomType1;

    //! Atom Type of atom 2
    std::string atomType2;

    //!
    double sigma;

    //!
    double epsilon;
};

// ============================================================
// Class : parameters
// ------------------------------------------------------------
// Class parameters
// ============================================================
class parameters
{
protected:

    /*!
       \brief parameters Constructor
       \param p elements pointer
    */
    parameters(elements* p);

    /*!
       \brief parameters Destructor
    */
    virtual ~parameters();

public:
    /*!
       \brief Returns the singleton for the parameter.  It will be used by all collections
       \param p elements pointer
     */
    static parameters* getInstance(elements* p);

    /*!
       \brief print atom types
    */
    void printAtomTypes();

    /*!
       \brief Add atom type
    */
    atomType*                addAtomType();

    /*!
       \brief Add atom type
       \param a atom type pointer
       \param n new atom type name
       \param g new group name
       \return atom type pointer
    */
    atomType*                addAtomType(atomType* a, std::string n,
                                                      std::string g);

    /*!
       \brief Add bond parameter
       \return bond parameter pointer
    */
    bondParam*               addBondParam();

    /*!
       \brief Add angle parameter
       \return angle parameter pointer
    */
    angleParam*              addAngleParam();

    /*!
       \brief Add torsion parameter
       \return torsion parameter pointer
    */
    torsionParam*            addTorsionParam();

    /*!
       \brief Add improper parameter
       \return improper parameters pointer
    */
    improperParam*           addImproperParam();

    /*!
       \brief Add H-Bond parameter
       \return H-Bond parameter pointer
    */
    hBondParam*              addHBondParam();

    /*!
       \brief Add equivalent atoms parameter
       \return equivalent atoms parameter pointer
    */
    equivalentAtomsParam*    addEquivalentAtomsParam();

    /*!
       \brief Calculate sigma and epsilon values for each pair of atom types
    */
    void                     calculateSigmaEpsilon();

    /*!
       \brief Set atomic number of atomType
       \param a atomType pointer
    */
    void                     setAtomNumber(atomType* a);

    /*!
       \brief Get atomic number using atom name
       \param at atom name
       \return atomic number
    */
    int                      getAtomicNum(const std::string& at);

    /*!
       \brief Has atom type using atom name
       \param a atom name
       \param g group name
       \return boolean
    */
    bool                     hasAtomType(const std::string& a, const std::string& g);

    /*!
       \brief Get atom type using atom name
       \param a atom name
       \return atomType pointer
    */
    atomType*                getAtomType(const std::string& a);

    /*!
       \brief Get atom types
       \return atomType vector
    */
    std::vector<atomType*>   getAtomTypes();

    /*!
       \brief Get bond parameters
       \return bondParam vector
    */
    std::vector<bondParam*>  getBondParams();

    /*!
       \brief Get angle parameters
       \return angleParams vector
    */
    std::vector<angleParam*> getAngleParams();

    /*!
       \brief Get torsion parameters
       \return torsionParams vector
    */
    std::vector<torsionParam*>    getTorsionParams();

    /*!
       \brief Get improper parameters
       \return torsionParams vector
    */
    std::vector<improperParam*>   getImproperParams();

    /*!
       \brief Get element symbol using atom type
       \param a atom type
       \return atom symbol
    */
    std::string              getAtomTypeSymbol(const std::string& a);

    /*!
       \brief Get atom hybridization using atom type
       \param a atom type
       \return atom hybridization
    */
    std::string              getAtomTypeHybridization(const std::string& a);

    /*!
       \brief Get bond parameter
       \param at1 atom type 1
       \param at2 atom type 2
    */
    bondParam*               getBondParam(const std::string& at1,
                                          const std::string& at2);

    /*!
       \brief Get bond parameter
       \param at1 atom type 1
       \param at2 atom type 2
       \param group name
    */
    bool                     hasBondParam(const std::string& at1,
                                          const std::string& at2,
                                          const std::string& g);

    /*!
       \brief Get angle parameter
       \param at1 atom type 1
       \param at2 atom type 2
       \param at3 atom type 3
    */
    angleParam*              getAngleParam(const std::string& at1,
                             const std::string& at2, const std::string& at3);

    /*!
       \brief Get angle parameter
       \param at1 atom type 1
       \param at2 atom type 2
       \param at3 atom type 3
       \param group name
    */
    bool                     hasAngleParam(const std::string& at1,
                             const std::string& at2, const std::string& at3,
                             const std::string& g);

    /*!
       \brief Get torsion parameter vector
       \param at1 atom type 1
       \param at2 atom type 2
       \param at3 atom type 3
       \param at4 atom type 4
    */
    std::vector<torsionParam*> getTorsionParamList(
                               const std::string& at1, const std::string& at2,
                               const std::string& at3, const std::string& at4);

    /*!
       \brief Remove Proline torsion
       \param t torsion parameters
    */
    void removeProlineTorsion(std::vector<torsionParam*>& t);

    /*!
       \brief Get improper parameter vector
       \param at1 atom type 1
       \param at2 atom type 2
       \param at3 atom type 3
       \param at4 atom type 4
       \param l order
    */
    std::vector<improperParam*> getImproperParamList(
                                const std::string& at1, const std::string& at2,
                                const std::string& at3, const std::string& at4,
                                std::vector<std::vector<int> >& l);

    /*!
       \brief
    */
    std::vector<equivalentAtomsParam*> getEquivalentAtomList();

    /*!
       \brief
    */
    LJ612SE* getLJ612SE(const std::string&, const std::string&);

protected:

    //! elements pointer
    elements*                          pElements;
    //! atom type vector
    std::vector<atomType*>             itsTypeList;

    //! bond parameter vector
    std::vector<bondParam*>            itsBondList;

    //! angle parameter vector
    std::vector<angleParam*>           itsAngleList;

    //! torsion parameter vector
    std::vector<torsionParam*>         itsTorsionList;

    //! improper parameter vector
    std::vector<improperParam*>        itsImproperList;

    //! H-Bond parameter vector
    std::vector<hBondParam*>           itsHBondList;

    //! equivalent atom parameter vector
    std::vector<equivalentAtomsParam*> itsEquivalentAtomsList;

    //! L-J 6-12 parameter vector
    std::vector<LJ612SE*>              itsLJ612SEList;

    //! atom type iterator
    typedef std::vector<atomType*>::iterator atomTypeIterator;

    //! bond parameter iterator
    typedef std::vector<bondParam*>::iterator bondParamIterator;

    //! angle parameter iterator
    typedef std::vector<angleParam*>::iterator angleParamIterator;

    //! torsion parameter iterator
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    //! improper parameter iterator
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    //! H-Bond parameter iterator
    typedef std::vector<hBondParam*>::iterator hBondParamIterator;

    //! L-J 6-12 parameter iterator
    typedef std::vector<LJ612SE*>::iterator LJ612SEIterator;

    //! atom type pointer
    atomType*                pAtomType;

    //! bond parameter pointer
    bondParam*               pBondParam;

    //! angle parameter pointer
    angleParam*              pAngleParam;

    //! torsion parameter pointer
    torsionParam*            pTorsionParam;

    //! improper parameter pointer
    improperParam*           pImproperParam;

    //! H-Bond parameter pointer
    hBondParam*              pHBondParam;

    //! equivalent atom pointer
    equivalentAtomsParam*    pEquivalentAtomsParam;

    //! L-J 6-12 parameter pointer
    LJ612SE*                 pLJ612SE;

    //! L-J 6-12 parameter boolean
    bool                     bLJ612SE;
};

} // MTKpp namespace

#endif // PARAMETERS_H
