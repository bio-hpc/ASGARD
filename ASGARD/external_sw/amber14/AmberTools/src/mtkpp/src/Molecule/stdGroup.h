/*!
   \file stdGroup.h
   \brief Container for standard fragments
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
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

#ifndef STDGROUP_h
#define STDGROUP_h

#include <iostream>
#include <string>
#include <vector>

namespace MTKpp
{

class stdLibrary;
class stdFrag;
class molecule;

// ============================================================
// Class : stdGroup()
// ------------------------------------------------------------
// Class stdGroup - Container for stdFrag.
// ============================================================

class stdGroup
{
public:

    /*!
       \brief stdGroup Constructor
       \param parent stdLibrary pointer
    */
    stdGroup(stdLibrary *parent = 0);

    //! stdGroup Destructor
    virtual ~stdGroup();

    /*!
       \brief Get stdLibrary which stdGroup is apart of
       \return stdLibrary pointer
    */
    stdLibrary*              getParent();

    /*!
       \brief Set name of stdGroup
       \param name stdGroup name
    */
    void                     setName(const std::string& name);

    /*!
       \brief Get name of stdGroup
       \return stdGroup name
    */
    std::string              getName();

    /*!
       \brief Set info string of stdGroup
       \param info stdGroup info
    */
    void                     setInfo(const std::string& info);

    /*!
       \brief Get info of stdGroup
       \return stdGroup info
    */
    std::string              getInfo();

    /*!
       \brief Add stdFrag to stdGroup
       \return stdFrag
    */
    virtual stdFrag*         addStdFrag();

    /*!
       \brief Add stdFrag to stdGroup
       \param f stdFrag pointer
       \return stdFrag
    */
    virtual stdFrag*         addStdFrag(stdFrag* f);

    /*!
       \brief Add stdFrag to stdGroup from another stdGroup
       \param s stdFrag pointer
       \return stdFrag
    */
    void                     appendStdFrag(stdFrag* s);

    /*!
       \brief Get stdFrag by name
       \param name stdFrag name
    */
    virtual stdFrag*         getStdFrag(const std::string& name);

    /*!
       \brief Has stdGroup a certain stdFrag
       \param name stdFrag name
       \return boolean
    */
    bool                     hasStdFrag(const std::string& name);

    /*!
       \brief Get list of stdFrags
       \return list of stdFrags
    */
    std::vector<stdFrag*>    getStdFragList();

    /*!
       \brief generate simple fingerprints
    */
    void                     generateSimpleFP();

    /*!
       \brief generate adjacency matrices
    */
    void                     generateAdjMatrices();

    /*!
       \brief Has stdMolecule
       \return boolean
    */
    bool                     hasStdMolecule();

    /*!
       \brief Set stdMolecule
       \param m pointer
    */
    void                     setStdMolecule(molecule* m);

    /*!
       \brief Get stdMolecule
       \return molecule pointer
    */
    molecule*                getStdMolecule();

    /*!
       \brief Get getCharge
       \return molecule charge
    */
    double                   getCharge();

protected:
    //! stdLibrary pointer
    stdLibrary*              pParent;

    //! stdGroup name
    std::string              itsName;

    //! stdGroup information
    std::string              itsInfo;

    //! stdFrag iterator
    typedef std::vector<stdFrag*>::iterator stdFragIterator;

    //! stdFrag list
    std::vector<stdFrag*>    itsStdFragList;

    //! stdFrag pointer
    stdFrag*                 pStdFrag;

    //!
    molecule*                pMolecule;
};

} // MTKpp namespace

#endif // STDGROUP_H

