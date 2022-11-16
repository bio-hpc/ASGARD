/*!
   \file element.h
   \brief Container for element information
   \author Martin Peters

   $Date: 2010/03/29 20:43:22 $
   $Revision: 1.9 $

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

#ifndef ELEMENT_h
#define ELEMENT_h

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "Utils/constants.h"

namespace MTKpp
{

// ============================================================
// Struct : element
// ------------------------------------------------------------
/*! 
   \struct element
   \brief element information
*/
// ============================================================
struct element
{
    //! Atomic number
    int number;

    //! Symbol
    std::string symbol;

    //! Name
    std::string name;

    //! Mass
    double mass;

    //! Periodic table group
    int group;

    //! Periodic table period
    int period;

    //! Color
    double red;

    //! Color
    double green;

    //! Color
    double blue;

    //! Standard valence
    int valence;

    //! Standard full shell
    int filledShell;

    //! Covalent radius
    double covalentRadius;

    //! van der Waals radius
    double vdWRadius;

    //! Pauling electronegativity value
    double paulingEN;

    //! Semiempirical Hamiltonians available
    std::vector<std::string> seHams;
};

// ============================================================
// Class : elements()
// ------------------------------------------------------------
/*! 
   \class elements
   \brief Container for element information
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class elements
{
public:
    //! Class constructor
    elements();

    //! Class destructor
    virtual ~elements();

    /*!
       \brief Add element
       \return element pointer
    */
    element* addElement();

    /*!
       \brief Set element name
       \param s element name 
    */
    void     setElementName(const std::string s);

    /*!
       \brief Get element by name
       \param s element symbol
       \return element pointer
    */
    element* getElement(const std::string s);

    /*!
       \brief Returns atomic number
       \param s element symbole
       \return atomic number
    */
    int getElementNumber(const std::string s);

    /*!
       \brief Get element name
       \param s element symbol
       \return element name
    */
    std::string getElementName(const std::string s);

    /*!
       \brief Get element mass
       \param s element symbol
       \return element mass
    */
    double getElementMass(const std::string s);

    /*!
       \brief Do we have parameters of this element in a certain Hamiltonian
       \param e element
       \param h Hamiltonian
       \return yes/no
    */
    bool hasSEHamiltonian(const std::string e, const std::string h);

protected:
    //! element map
    std::map<std::string, element*> itsElementMap;

    //! element iterator
    typedef std::map<std::string, element*>::iterator ElementMapIterator;

    //! element pointer
    element* pElement;

    //! string iterator
    typedef std::vector<std::string>::iterator strIterator;
};

} // MTKpp namespace

#endif // ELEMENT_H
