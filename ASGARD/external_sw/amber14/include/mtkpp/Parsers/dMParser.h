/*!
   \file dMParser.h
   \brief sheet/tables parser
   \author Martin B. Peters

   $Date: 2010/08/11 21:11:00 $
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

#ifndef DMPARSER_H
#define DMPARSER_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include "baseParser.h"

#include "Statistics/table.h"

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

class sheet;

// ============================================================
// Class : dMParser()
// ------------------------------------------------------------
/*!
   \class dMParser
   \brief sheet/table parser
   \author Martin Peters
*/
// ============================================================
#ifdef USE_QT
class dMParser : public QObject, public baseParser
{
    Q_OBJECT
#else // USE_XERCES or USE_TINYMXL
class dMParser : public baseParser
{
#endif

public:
    /*!
       \brief dMParser Constructor
    */
    dMParser();

    //! dMParser Destructor
    virtual ~dMParser();

    /*!
       \brief Read xml file
       \param s sheet pointer
       \param f filename
       \return success
    */
    int read(sheet* s, std::string f);

    /*!
       \brief Import txt file
         Line
         1st tableName nRows nColumns
         2nd Column Labels
         3rd Row Label data ...
         ...
         Nth
       \param s sheet pointer
       \param f filename
       \return success
    */
    int import(sheet* s, std::string f);

    /*!
       \brief Write xml file
       \param s sheet pointer
       \param f filename
       \return success
    */
    int write(sheet* s, std::string f, bool bComments = false);

#ifdef USE_QT
protected:

    /*!
       \brief table filler
       \param d dom node
    */
    void tableFiller(QDomNode d);

    /*!
       \brief row and column filler
       \param d dom node
    */
    void rowColFiller(QDomNode d);

    /*!
       \brief cell filler
       \param d dom node
       \param i int
    */
    void cellFiller(QDomNode d, int i);

    /*!
       \brief integer cell filler
       \param d dom node
       \param i int
    */
    void cellIntFiller(QDomNode d, int i);

    /*!
       \brief Write group
       \param doc dom document
       \param t table
    */
    void writeTable(QDomDocument doc, table<double>* t);

    /*!
       \brief Write group
       \param doc dom document
       \param t table
    */
    void writeIntTable(QDomDocument doc, table<int>* t);

#endif // USE_QT

#ifdef USE_XERCES
protected:

    /*!
       \brief 
       \param d dom node
    */
    void sheetFiller(DOMNode* d);

    /*!
       \brief 
       \param d dom node
    */
    void tableFiller(DOMNode* d);

    /*!
       \brief 
       \param d dom node
    */
    void rowColFiller(DOMNode* d);


    /*!
       \brief Cell filler
       \param d dom node
       \param i cell index
    */
    void cellFiller(DOMNode* d, int i);

    /*!
       \brief Cell filler
       \param d dom node
       \param i cell index
    */
    void cellIntFiller(DOMNode* d, int i);

    /*!
       \brief Write table
       \param doc dom document
       \param t table
    */
    void writeTable(XERCES_CPP_NAMESPACE::DOMDocument* doc, table<double>* t);

    /*!
       \brief Write table
       \param doc dom document
       \param t table
    */
    void writeIntTable(XERCES_CPP_NAMESPACE::DOMDocument* doc, table<int>* t);

#endif // USE_XERCES

#ifdef USE_TINYXML
protected:

    /*!
       \brief Write table
       \param root dom element
       \param t table
    */
    void writeTable(TiXmlElement* root, table<double>* t);

    /*!
       \brief Write table
       \param root dom element
       \param t table
    */
    void writeIntTable(TiXmlElement* root, table<int>* t);

#endif // USE_TINYXML

protected:
    //! sheet pointer
    sheet* mySheet;

    //! table containing doubles
    table<double>* myDoubleTable;

    //! table containing integers
    table<int>* myIntTable;

    //! Current table type
    std::string currentType;
};

} // MTKpp namespace

#endif // DMPARSER_H
