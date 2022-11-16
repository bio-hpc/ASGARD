/*!
   \file mtkppParser.h
   \brief Parses MTK++ State files
   \author Martin Peters

   Reads and writes MTK++ State files

   $Date: 2010/08/19 13:48:48 $
   $Revision: 1.10 $

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

#ifndef MTKPPPARSER_H
#define MTKPPPARSER_H

#include <vector>
#include <map>

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "Utils/constants.h"
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
class molecule;
class submolecule;
class atom;
class vector3d;
struct Bond;
struct Angle;
struct Torsion;
struct Improper;
class metalCenter;

// ============================================================
// Class : mtkppParser()
// ------------------------------------------------------------
/*!
   \class mtkppParser
   \brief Reads and writes MTK++ State xml files
   \author Martin Peters
*/
// ============================================================

#ifdef USE_QT
class mtkppParser : public QObject, public baseParser
{
    Q_OBJECT

#else // USE_XERCES and USE_TINYXML
class mtkppParser : public baseParser
{
#endif

friend class stdLibParser;

public:
    /*!
       \brief mtkppParser Constructor
    */
    mtkppParser();

    //! mtkppParser Destructor
    virtual ~mtkppParser();

    /*!
       \brief Read MTK++ xml file
       \param i input file
       \param c collection pointer
    */
    void           Read(const std::string &i, collection* c);

    /*!
       \brief Write MTK++ xml file
       \param o output file
       \param c collection pointer
    */
    void           Write(const std::string &o, collection* c);


#ifdef USE_QT

    /*!
       \brief Read molecule
       \param node dom node
       \param mol molecule
    */
    void readMolecule(QDomNode node, molecule* mol);

protected: // functions

    /*!
       \brief Read submolecule
       \param node dom node
       \param mol molecule
    */
    void readSubmolecule(QDomNode node, molecule* mol);

    /*!
       \brief Read atom
       \param node dom node
       \param mol molecule
    */
    void readAtom(QDomNode node, submolecule* subMol);

    /*!
       \brief Read bond
       \param node dom node
       \param mol molecule
    */
    void readBond(QDomNode node, molecule* mol);

    /*!
       \brief Read angle
       \param node dom node
       \param mol molecule
    */
    void readAngle(QDomNode node, molecule* mol);

    /*!
       \brief Read angle
       \param node dom node
       \param mol molecule
    */
    void readTorsion(QDomNode node, molecule* mol);

    /*!
       \brief Write molecule
       \param doc dom document
       \param mol molecule
    */
    void writeMolecule(QDomDocument doc, molecule* mol);

    /*!
       \brief Write molecule
       \param doc dom document
       \param groupElem dom element
       \param mol molecule
    */
    void writeMolecule(QDomDocument doc, QDomElement groupElem, molecule* mol);

    /*!
       \brief Write molecule
       \param doc dom document
       \param moleculeElem dom element
       \param mol molecule
    */
    void writeSubMolecule(QDomDocument doc, QDomElement moleculeElem, submolecule* smol);

    /*!
       \brief Write atom
       \param doc dom document
       \param smolElem dom element
       \param smol submolecule pointer
       \param at atom pointer
    */
    void writeAtom(QDomDocument doc, QDomElement smolElem, submolecule* smol, atom* at);

    /*!
       \brief Write bonds
       \param doc dom document
       \param moleculeElem molecule element
       \param bonds bond map
    */
    void writeBonds(QDomDocument doc, QDomElement moleculeElem, std::map<int, Bond*> bonds);

    /*!
       \brief Write angles
       \param doc dom document
       \param moleculeElem molecule element
       \param angles angle map
    */
    void writeAngles(QDomDocument doc, QDomElement moleculeElem, std::map<ULONG_KIND, Angle*> angles);

    /*!
       \brief Write angles
       \param doc dom document
       \param moleculeElem molecule element
       \param torsions torsion map
    */
    void writeTorsions(QDomDocument doc, QDomElement moleculeElem, std::map<ULONG_KIND, Torsion*> torsions);

#elif defined(USE_XERCES)

    /*!
       \brief Read molecule
       \param node dom node
       \param mol molecule
    */
    void readMolecule(DOMNode* molNode, molecule* mol);

protected: // functions

    // todo
    /*!
       \brief Read submolecule
       \param smolNode dom node
       \param mol molecule
    */
    void readSubmolecule(DOMNode* smolNode, molecule* mol);

    /*!
       \brief Read atom
       \param atomNode dom node
       \param mol molecule
    */
    void readAtom(DOMNode* atomNode, submolecule* subMol);

    /*!
       \brief Read bond
       \param bondNode dom node
       \param mol molecule
    */
    void readBond(DOMNode* bondNode, molecule* mol);

    /*!
       \brief Read angle
       \param angleNode dom node
       \param mol molecule
    */
    void readAngle(DOMNode* angleNode, molecule* mol);

    /*!
       \brief Read angle
       \param torNode dom node
       \param mol molecule
    */
    void readTorsion(DOMNode* torNode, molecule* mol);

    /*!
       \brief Write molecule
       \param doc dom document
       \param mol molecule
    */
    void writeMolecule(XERCES_CPP_NAMESPACE::DOMDocument* doc, molecule* mol);

    /*!
       \brief Write submolecule
       \param doc dom document
       \param moleculeElem molecule element
       \param smol submolecule pointer
    */
    void writeSubMolecule(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, submolecule* smol);

    /*!
       \brief Write atom
       \param doc dom document
       \param smolElem dom element
       \param smol submolecule pointer
       \param at atom pointer
    */
    void writeAtom(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* smolElem, submolecule* smol, atom* at);

    /*!
       \brief Write metal center
       \param doc dom document
       \param pMetalCenter metal center pointer
    */
    void writeMetalCenter(XERCES_CPP_NAMESPACE::DOMDocument* doc, metalCenter* pMetalCenter);

    /*!
       \brief Write bonds
       \param doc dom document
       \param moleculeElem molecule element
       \param bonds bond map
    */
    void writeBonds(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, std::map<int, Bond*> bonds);

    /*!
       \brief Write angles
       \param doc dom document
       \param moleculeElem molecule element
       \param angles angle map
    */
    void writeAngles(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, std::map<ULONG_KIND, Angle*> angles);

    /*!
       \brief Write angles
       \param doc dom document
       \param moleculeElem molecule element
       \param torsions torsion map
    */
    void writeTorsions(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, std::map<ULONG_KIND, Torsion*> torsions);

#else // USE_TINYXML
    /*!
       \brief Read molecule
       \param node dom node
       \param mol molecule
    */
    void readMolecule(TiXmlNode* molNode, molecule* mol);

protected: // functions

    /*!
       \brief Read submolecule
       \param node dom node
       \param mol molecule
    */
    void readSubmolecule(TiXmlNode* smolNode , molecule* mol);

    /*!
       \brief Read atom
       \param node dom node
       \param mol molecule
    */
    void readAtom(TiXmlNode* atomNode, submolecule* subMol);

    /*!
       \brief Read bond
       \param node dom node
       \param mol molecule
    */
    void readBond(TiXmlNode* bondNode, molecule* mol);

    /*!
       \brief Read angle
       \param node dom node
       \param mol molecule
    */
    void readAngle(TiXmlNode* angleNode, molecule* mol);

    /*!
       \brief Read angle
       \param node dom node
       \param mol molecule
    */
    void readTorsion(TiXmlNode* torNode, molecule* mol);

    /*!
       \brief Write molecule
       \param root TiXmlElement object
       \param mol molecule
    */
    void writeMolecule(TiXmlElement* root, molecule* mol);

    /*!
       \brief Write submolecule
       \param root TiXmlElement object
       \param mol molecule
    */
    void writeSubMolecule(TiXmlElement* root, submolecule* smol);

    /*!
       \brief Write atom
       \param root TiXmlElement object
       \param smol submolecule pointer
       \param at atom pointer
    */
    void writeAtom(TiXmlElement* root, submolecule* smol, atom* at);

    /*!
       \brief Write metal center
       \param doc dom document
       \param pMetalCenter metal center pointer
    */
    void writeMetalCenter(TiXmlElement* root, metalCenter* pMetalCenter);

    /*!
       \brief Write bonds
       \param root TiXmlElement object
       \param bonds bond map
    */
    void writeBonds(TiXmlElement* root, std::map<int, Bond*> bonds);

    /*!
       \brief Write angles
       \param root TiXmlElement object
       \param angles angle map
    */
    void writeAngles(TiXmlElement* root, std::map<ULONG_KIND, Angle*> angles);

    /*!
       \brief Write angles
       \param root TiXmlElement object
       \param torsions torsion map
    */
    void writeTorsions(TiXmlElement* root, std::map<ULONG_KIND, Torsion*> torsions);

#endif // USE_QT/USE_XERCES

protected: // data

    //////////////////////////////////////////
    // - pointers to the backbone objects - //
    //////////////////////////////////////////

    //! collection pointer
    collection*    pCollection;

    //! atom pointer
    atom*          pAtom;

    //! atom pointer
    atom*          pAtom1;

    //! atom pointer
    atom*          pAtom2;

    //! Bond pointer
    Bond*          pBond;

    //! Angle pointer
    Angle*         pAngle;

    //! Torsion pointer
    Torsion*       pTorsion;

    //! Improper pointer
    Improper*      pImproper;

    //! coordinate pointer
    vector3d*      pCoord1;

    //! coordinate pointer
    vector3d*      pCoord2;

    //! molecule iterator
    typedef std::vector<molecule*>::iterator          moleculeIterator;

    //! submolecule iterator
    typedef std::vector<submolecule*>::iterator       submoleculeIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator              atomIterator;

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator            BondMapIterator;

    //! Angle map iterator
    typedef std::map<ULONG_KIND, Angle*>::iterator    AngleMapIterator;

    //! Torsion map iterator
    typedef std::map<ULONG_KIND, Torsion*>::iterator  TorsionMapIterator;

    //! Improper map iterator
    typedef std::map<int, Improper*>::iterator        ImproperMapIterator;

    //! property map iterator
    typedef std::map<std::string, double>::iterator   PropertyMapIterator;

    //! property map iterator
    typedef std::map<std::string, int>::iterator      intPropertyMapIterator;
};

} // MTKpp namespace

#endif // MTKPPPARSER_H
