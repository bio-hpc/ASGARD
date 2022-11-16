/*!
   \file paramParser.h
   \brief Parses parameter xml files using xercesc
   \author Martin Peters

   $Date: 2010/08/11 21:11:00 $
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

#ifndef PARAMPARSER_H
#define PARAMPARSER_H

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "baseParser.h"

#ifdef USE_QT
#include <QtCore>
#include <QtXml>
#elif defined(USE_XERCES)
#include "domErrorHandler.h"
#else // USE_TINYXML
#include <sstream>
#include "tinyxml/tinyxml.h"
#endif

namespace MTKpp
{

class parameters;
struct atomType;
struct bondParam;
struct angleParam;
struct torsionParam;
struct improperParam;
struct hBondParam;
struct equivalentAtomsParam;

// ============================================================
// Class : paramParser()
// ------------------------------------------------------------
/*!
   \class paramParser
   \brief Reads parameter xml files
   \author Martin Peters
*/
// ============================================================

#ifdef USE_QT
class paramParser : public QObject, public baseParser
{
    Q_OBJECT

#else // USE_XERCES or USE_TINYXML
class paramParser : public baseParser
{
#endif

public:

    /*!
       \brief paramParser Constructor
       \param c parameters pointer
    */
    paramParser(parameters* c);

    //! paramParser Destructor
    virtual ~paramParser();

    /*!
       \brief Reads parameter xml files
       \param fileName parameter xml file
       \return boolean
    */
    int Read(std::string fileName);

    /*!
       \brief Write standard library xml file
       \param i parameter xml file
       \param g group name
    */
    int Write(std::string i, std::string g);

protected:
#ifdef USE_QT

    /*!
       \brief Read types
       \param d dom node
    */
    void typeFiller(QDomNode d);

    /*!
       \brief Read bonds
       \param d dom node
    */
    void bondFiller(QDomNode d);

    /*!
       \brief Read angles
       \param d dom node
    */
    void angleFiller(QDomNode d);

    /*!
       \brief Read torsions
       \param d dom node
    */
    void torsionFiller(QDomNode d);

    /*!
       \brief Read impropers
       \param d dom node
    */
    void improperFiller(QDomNode d);

    /*!
       \brief Read hbonds
       \param d dom node
    */
    void hBondFiller(QDomNode d);

    /*!
       \brief Read equivalent atoms
       \param d dom node
    */
    void equivalAtomFiller(QDomNode d);
#endif // USE_QT

#ifdef USE_XERCES

    /*!
       \brief Read parameters
       \param d dom node
    */
    void paramFiller(DOMNode* d);

    /*!
       \brief Read types
       \param d dom node
    */
    void typeFiller(DOMNode* d);

    /*!
       \brief Read bonds
       \param d dom node
    */
    void bondFiller(DOMNode* d);

    /*!
       \brief Read angles
       \param d dom node
    */
    void angleFiller(DOMNode* d);

    /*!
       \brief Read torsions
       \param d dom node
    */
    void torsionFiller(DOMNode* d);

    /*!
       \brief Read impropers
       \param d dom node
    */
    void improperFiller(DOMNode* d);

    /*!
       \brief Read hbonds
       \param d dom node
    */
    void hBondFiller(DOMNode* d);

    /*!
       \brief Read equivalent atoms
       \param d dom node
    */
    void equivalAtomFiller(DOMNode* d);
#endif // USE_XERCES

    /*!
       \brief Update equivalent atoms
    */
    void updateEquivalentAtoms();

protected:

    //! group name
    std::string              groupName;

    //! parameters pointer
    parameters*              pParameters;

    //! atomType pointer
    atomType*                pAtomType;

    //! bondParam pointer
    bondParam*               pBondParam;

    //! angleParam pointer
    angleParam*              pAngleParam;

    //! torsionParam pointer
    torsionParam*            pTorsionParam;

    //! improperParam pointer
    improperParam*           pImproperParam;

    //! hBondParam pointer
    hBondParam*              pHBondParam;

    //! equivalentAtomsParam pointer
    equivalentAtomsParam*    pEquivalentAtomsParam;
};

} // MTKpp namespace

#endif // PARAMPARSER_H
