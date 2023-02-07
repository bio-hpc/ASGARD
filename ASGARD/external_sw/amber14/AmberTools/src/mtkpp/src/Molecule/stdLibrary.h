/*!
   \file stdLibrary.h
   \brief Container for standard groups
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
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

#ifndef STDLIBRARY_h
#define STDLIBRARY_h

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace MTKpp
{

class stdGroup;
class stdFrag;

// ============================================================
// Class : stdLibrary()
// ------------------------------------------------------------
// Class stdLibrary - Container for stdGroup.
// ============================================================

class stdLibrary
{
protected:
    /*!
       \brief stdLibrary Constructor
    */
    stdLibrary();

    /*!
       \brief Standard Library Destructor
    */
    virtual ~stdLibrary();

public:
    /*
     * Returns the singleton for the stdLibrary.  It will be used by all collections
     *
     */
    static stdLibrary* getInstance();

    /*!
       \brief Set name of library
       \param n name of library
    */
    void                     setName(const std::string& n);

    /*!
       \brief Get name of library
       \return name of library
    */
    std::string              getName();

    /*!
       \brief Add a stdGroup to the library
       \return stdGroup pointer
    */
    virtual stdGroup*        addStdGroup();

    /*!
       \brief Get stdGroup from library by name
       \param n stdGroup name
       \return stdGroup pointer
    */
    virtual stdGroup*        getStdGroup(const std::string& n);

    /*!
       \brief Get stdGroups from library
       \return vector of stdGroup pointers
    */
    std::vector<stdGroup*>   getStdGroupList();

    /*!
       \brief Get stdFrag from library by name
       \param fragName stdFrag Name
       \return stdFrag pointer
    */
    stdFrag*                 getStdFrag(const std::string& fragName);

    /*!
       \brief Get stdFrag from library by fragment name and group name
       \param fragName stdFrag Name
       \param grpName stdGroup name or stdFrag type
       \return stdFrag pointer
    */
    stdFrag*                 getStdFrag(const std::string& fragName,
                                        const std::string& grpName);

    /*!
       \brief Get number of stdFrags in the library
       \return number of stdFrags
    */
    int                      getNumberStdFrag();

    /*!
       \brief Set one letter code of amino acid or base
       \param lll 3 letter code
       \param l 1 letter code
    */
    void                     setL(std::string lll, std::string l);

    /*!
       \brief Get one letter code of amino acid or base
       \param lll 3 leter code
       \return 1 letter code
    */
    std::string              getL(std::string lll);

    /*!
       \brief generate simple fingerprints
    */
    void                     generateSimpleFP();

    /*!
       \brief generate adjacency matrices
    */
    void                     generateAdjMatrices();

    /*!
       \brief generate atom kind matrices
    */
    void                     generateAtomKinds();

    bool hasSimpleFPGenerated() {return this->bSimpleFPGenerated;}
    bool hasAdjMatricesGenerated() {return this->bAdjMatricesGenerated;}
    bool hasAtomKindsGenerated() {return this->bAtomKindsGenerated;}

protected:
    //! stdLibrary name
    std::string              itsName;

    //! stdGroup iterator
    typedef std::vector<stdGroup*>::iterator stdGroupIterator;

    //! stdGroup list
    std::vector<stdGroup*>   itsStdGroupList;

    //! stdGroup pointer
    stdGroup*                pStdGroup;

    //! stdFrag pointer
    stdFrag*                 pStdFrag;

    //! 3L code to 1L code for proteins and dna
    std::map<std::string, std::string>   LLL2L;

    //!
    bool bSimpleFPGenerated;

    //!
    bool bAdjMatricesGenerated;

    //!
    bool bAtomKindsGenerated;
};

} // MTKpp namespace

#endif // STDLIBRARY_H

