/*!
   \file elementParser.h
   \brief Parses elements xml file using XERCES-C or Trolltech's Qt
   \author Martin Peters

   $Date: 2010/08/11 21:11:00 $
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

#ifndef ELEMENTPARSER_H
#define ELEMENTPARSER_H

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
#include <QObject>
#elif defined(USE_XERCES)
#include "domErrorHandler.h"
#else // USE_TINYXML
#include <sstream>
#include "tinyxml/tinyxml.h"
#endif

namespace MTKpp
{

class elements;
struct element;

// ============================================================
// Class : elementParser()
// ------------------------------------------------------------
/*!
   \class elementParser
   \brief Reads element xml file
   \author Martin Peters
*/
// ============================================================
#ifdef USE_QT
class elementParser : public QObject, public baseParser
{
    Q_OBJECT

#else // USE_XERCES or USE_TINYXML
class elementParser : public baseParser
{
#endif

public:

    /*!
       \brief elementParser Constructor
    */
     elementParser(elements *e);

    //! elementParser Destructor
     virtual ~elementParser();

    /*!
       \brief Read element xml file
       \param i Input file
       \return boolean
    */
     int Read(std::string i);

#ifdef USE_XERCES
protected:
    /*!
       \brief Read elements section
       \param d dom node
     */
     void elementsFiller(DOMNode* d);

    /*!
       \brief Reads element section
       \param d dom node
     */
     void elementFiller(DOMNode* d);
#endif // USE_XERCES

protected:

    //! elements pointer
    elements*          pElements;

    //! element pointer
    element*           pElement;
};

} // MTKpp namespace

#endif // ELEMENTPARSER_H
