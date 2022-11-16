/*!
   \file stdLibParser.h
   \brief Parses standard library xml files using xercesc
   \author Martin Peters

   $Date: 2010/08/11 21:06:59 $
   $Revision: 1.18 $

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

#ifndef STDLIBPARSER_H
#define STDLIBPARSER_H

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

class collection;
class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;
struct stdBond;
struct stdImproper;
struct stdLoop;
struct stdAlias;
struct stdRing;
struct stdFeature;
struct stdFuncGroup;
class parameters;

#ifdef USE_QT
class stdLibParser : public QObject, public baseParser
{
    Q_OBJECT
#else // USE_XERCES or USE_TINYXML
class stdLibParser : public baseParser
{
#endif

public:

    /*!
       \brief baseParser Constructor
       \param c collection pointer
       \param s stdLibrary pointer
       \param p parameters pointer
    */
    stdLibParser(collection* c, stdLibrary* s, parameters* p);

    /*!
       \brief baseParser Constructor
       \param s stdLibrary pointer
       \param p parameters pointer
    */
    stdLibParser(stdLibrary* s, parameters* p);

    /*!
       \brief baseParser Constructor
       \param s stdLibrary pointer
    */
    stdLibParser(stdLibrary* s);

    //! baseParser Destructor
    virtual ~stdLibParser();

    /*!
       \brief Read standard library xml file
       \param i Standard library file
    */
    void Read(std::string i);

    /*!
       \brief Write standard library xml file
       \param i Standard library file
    */
    void Write(std::string i);

    /*!
       \brief Write standard library xml file
       \param i Standard library file
       \param g group name
    */
    void Write(std::string i, std::string g);

    /*!
       \brief Write standard library xml file
       \param i Standard library file
       \param g group name
       \param f fragment name
    */
    void Write(std::string i, std::string g, std::string f);

#ifdef USE_TINYXML
protected:

    /*!
       \brief Standard group writer
       \param root TiXmlElement object
    */
    void writeGroup(TiXmlElement* root, stdGroup* stdGp);

    /*!
       \brief Standard group writer
       \param root TiXmlElement object
       \param fragName standard fragment name
    */
    void writeGroup(TiXmlElement* root, stdGroup* stdGp, std::string fragName);

    /*!
       \brief Standard fragment writer
       \param root TiXmlElement object
       \param stdFg standard fragment object
    */
    void writeFrag(TiXmlElement* root, stdFrag* stdFg);

#endif // USE_TINYXML

#ifdef USE_QT
protected:

    /*!
       \brief Standard group reader
       \param d QDomNode object
    */
    void groupFiller(QDomNode d);

    /*!
       \brief Standard fragment reader
       \param d QDomNode object
    */
    void fragmentFiller(QDomNode d);

    /*!
       \brief Standard atom reader
       \param d QDomNode object
    */
    void atomFiller(QDomNode d);

    /*!
       \brief Standard alias reader
       \param d QDomNode object
    */
    void aliasFiller(QDomNode d);

    /*!
       \brief Standard improper reader
       \param d QDomNode object
    */
    void improperFiller(QDomNode d);

    /*!
       \brief Standard loop reader
       \param d QDomNode object
    */
    void loopFiller(QDomNode d);

    /*!
       \brief Standard ring reader
       \param d QDomNode object
    */
    void ringFiller(QDomNode d);

    /*!
       \brief Standard feature reader
       \param d QDomNode object
    */
    void featureFiller(QDomNode d);

    /*!
       \brief Standard functional group reader
       \param d QDomNode object
    */
    void funcGroupFiller(QDomNode d);

    /*!
       \brief Connection Points reader
       \param d QDomNode object
    */
    void connPtsFiller(QDomNode d);

    /*!
       \brief Connection Torsions reader
       \param d QDomNode object
    */
    void connTorFiller(QDomNode d);

    /*!
       \brief Rotatable Bond reader
       \param d QDomNode object
    */
    void rotBondFiller(QDomNode d);

    /*!
       \brief Write group
       \param doc QDomDocument object
       \param stdGp standard group
    */
    void writeGroup(QDomDocument doc, stdGroup* stdGp);

    /*!
       \brief Write group
       \param doc QDomDocument object
       \param stdGp standard group
       \param fragName standard fragment name
    */
    void writeGroup(QDomDocument doc, stdGroup* stdGp, std::string fragName);

    /*!
       \brief Write fragment
       \param doc QDomDocument object
       \param groupElem dom element
       \param stdFg standard fragment
    */
    void writeFrag(QDomDocument doc, QDomElement groupElem, stdFrag* stdFg);

    /*!
       \brief Write atom
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdFg standard fragment
       \param stdAt standard atom
    */
    void writeAtom(QDomDocument doc, QDomElement fragElem, stdFrag* stdFg, stdAtom* stdAt);

    /*!
       \brief Write loop
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdLp standard loop
    */
    void writeLoop(QDomDocument doc, QDomElement fragElem, stdLoop* stdLp);

    /*!
       \brief Write alias
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdAl standard alias
    */
    void writeAlias(QDomDocument doc, QDomElement fragElem, stdAlias* stdAl);

    /*!
       \brief Write improper
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdIm standard improper
    */
    void writeImproper(QDomDocument doc, QDomElement fragElem, stdImproper* stdIm);

    /*!
       \brief Write ring
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdRg standard ring
    */
    void writeRing(QDomDocument doc, QDomElement fragElem, stdRing* stdRg);

    /*!
       \brief Write feature
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdFt standard feature
    */
    void writeFeature(QDomDocument doc, QDomElement fragElem, stdFeature* stdFt);

    /*!
       \brief Write funcGroup
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdFg standard functional group
    */
    void writeFuncGroup(QDomDocument doc, QDomElement fragElem, stdFuncGroup* stdFg);

    /*!
       \brief Write connection points
       \param doc QDomDocument object
       \param fragElem dom element
       \param stdFg standard fragment
    */
    void writeConnPts(QDomDocument doc, QDomElement fragElem, stdFrag* stdFg);

#endif // USE_QT

#ifdef USE_XERCES
protected:

    /*!
       \brief Standard library reader
       \param d dom node
    */
    void stdLibFiller(DOMNode* d);

    /*!
       \brief Standard group reader
       \param d dom node
    */
    void groupFiller(DOMNode* d);

    /*!
       \brief Standard fragment reader
       \param d dom node
    */
    void fragmentFiller(DOMNode* d);

    /*!
       \brief Standard atom reader
       \param d dom node
       \todo electron to kcal for MM calculations
    */
    void atomFiller(DOMNode* d);

    /*!
       \brief Standard alias reader
       \param d dom node
    */
    void aliasFiller(DOMNode* d);

    /*!
       \brief Standard improper reader
       \param d dom node
    */
    void improperFiller(DOMNode* d);

    /*!
       \brief Standard loop reader
       \param d dom node
    */
    void loopFiller(DOMNode* d);

    /*!
       \brief Standard ring reader
       \param d dom node
    */
    void ringFiller(DOMNode* d);

    /*!
       \brief Standard feature reader
       \param d dom node
    */
    void featureFiller(DOMNode* d);

    /*!
       \brief Standard functional group reader
       \param d dom node
    */
    void funcGroupFiller(DOMNode* d);

    /*!
       \brief Connection Points reader
       \param d dom node
    */
    void connPtsFiller(DOMNode* d);

    /*!
       \brief Connection Torsions reader
       \param d dom node
    */
    void connTorFiller(DOMNode* d);

    /*!
       \brief Rotatable Bond reader
       \param d dom node
    */
    void rotBondFiller(DOMNode* d);

    /*!
       \brief Write group
       \param doc dom document
       \param stdGp standard group
    */
    void writeGroup(XERCES_CPP_NAMESPACE::DOMDocument* doc, stdGroup* stdGp);

    /*!
       \brief Write group
       \param doc dom document
       \param stdGp standard group
       \param fragName standard fragment name
    */
    void writeGroup(XERCES_CPP_NAMESPACE::DOMDocument* doc, stdGroup* stdGp, std::string fragName);

    /*!
       \brief Write fragment
       \param doc dom document
       \param groupElem dom element
       \param stdFg standard fragment
    */
    void writeFrag(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* groupElem, stdFrag* stdFg);

    /*!
       \brief Write atom
       \param doc dom document
       \param fragElem dom element
       \param stdFg standard fragment
       \param stdAt standard atom
    */
    void writeAtom(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFrag* stdFg, stdAtom* stdAt);

    /*!
       \brief Write loop
       \param doc dom document
       \param fragElem dom element
       \param stdLp standard loop
    */
    void writeLoop(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdLoop* stdLp);

    /*!
       \brief Write alias
       \param doc dom document
       \param fragElem dom element
       \param stdAl standard alias
    */
    void writeAlias(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdAlias* stdAl);

    /*!
       \brief Write improper
       \param doc dom document
       \param fragElem dom element
       \param stdIm standard improper
    */
    void writeImproper(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdImproper* stdIm);

    /*!
       \brief Write ring
       \param doc dom document
       \param fragElem dom element
       \param stdRg standard ring
    */
    void writeRing(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdRing* stdRg);

    /*!
       \brief Write feature
       \param doc dom document
       \param fragElem dom element
       \param stdFt standard feature
    */
    void writeFeature(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFeature* stdFt);

    /*!
       \brief Write funcGroup
       \param doc dom document
       \param fragElem dom element
       \param stdFg standard functional group
    */
    void writeFuncGroup(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFuncGroup* stdFg);

    /*!
       \brief Write connection points
       \param doc dom document
       \param fragElem dom element
       \param stdFg standard fragment
    */
    void writeConnPts(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFrag* stdFg);
#endif // USE_XERCES

protected:

    //! collection pointer
    collection*         pCollection;

    //! stdLibrary pointer
    stdLibrary*         pStdLibrary;

    //! parameters pointer
    parameters*         pParameters;

    //! stdGroup pointer
    stdGroup*           pStdGroup;

    //! stdFrag pointer
    stdFrag*            pStdFrag;

    //! stdAtom pointer
    stdAtom*            pStdAtom;

    //! stdBond pointer
    stdBond*            pStdBond;

    //! stdImproper pointer
    stdImproper*        pStdImproper;

    //! stdLoop pointer
    stdLoop*            pStdLoop;

    //! stdAlias pointer
    stdAlias*           pStdAlias;

    //! stdRing pointer
    stdRing*            pStdRing;

    //! stdFeature pointer
    stdFeature*         pStdFeature;

    //! stdFuncGroup pointer
    stdFuncGroup*       pStdFuncGroup;
};

} // MTKpp namespace

#endif // STDLIBPARSER_H
