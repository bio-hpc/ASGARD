/*!
   \file mtkppParser.cpp
   \brief Parses MTK++ State xml files using xercesc
   \author Martin Peters

   $Date: 2010/08/19 13:48:48 $
   $Revision: 1.8 $

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

#include <sstream>

#include "mtkppParser.h"
#include "StringManip.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/angle.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"

#include "Molecule/stdFrag.h"

#include "Molecule/metalCenter.h"

#ifdef USE_XERCES
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#endif // USE_XERCES

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : mtkppParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
mtkppParser::mtkppParser() {}

// ============================================================
// Function : mtkppParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
mtkppParser::~mtkppParser() {}

    // ---------------------------- //
    // -    Qt READ FUNCTIONS     - //
    // ---------------------------- //

#ifdef USE_QT
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses MTK++ State xml files using Qt
// ------------------------------------------------------------
void mtkppParser::Read(const std::string &fileName, collection* c)
{
    std::cout << " IMPLEMENT mtkppParser::Read using Qt " << std::endl;

    pCollection = c;
}

// ============================================================
// Function : readMolecule
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void mtkppParser::readMolecule(QDomNode moleculeNode, molecule* pMolecule)
{
    pCollection = pMolecule->getParent();

    QDomElement moleculeElement = moleculeNode.toElement();

    if (!moleculeElement.isNull()) {
      if (moleculeElement.hasAttribute("identity")) {
        QString att = moleculeElement.attribute("identity");
        pMolecule->setName(att.toStdString());
      }

      QDomNode node = moleculeNode.firstChild();
      while (!node.isNull()) {
        if (node.nodeName() == "submolecule" and node.hasAttributes()) {
          readSubmolecule(node, pMolecule);
        }

        if (node.nodeName() == "bond" and node.hasAttributes()) {
          readBond(node, pMolecule);
        }

        if (node.nodeName() == "angle" and node.hasAttributes()) {
          readAngle(node, pMolecule);
        }

        if (node.nodeName() == "torsion" and node.hasAttributes()) {
          readTorsion(node, pMolecule);
        }

        node = node.nextSibling();
      }
    }
}

// ============================================================
// Function : readSubmolecule
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void mtkppParser::readSubmolecule(QDomNode submoleculeNode, molecule* pMolecule)
{
    QDomElement submoleculeElement = submoleculeNode.toElement();
    if (!submoleculeElement.isNull()) {
      submolecule* pSubmolecule = pMolecule->addSubMolecule();

      if (submoleculeElement.hasAttribute("name")) {
        QString att = submoleculeElement.attribute("name");
        pSubmolecule->setName(att.toStdString());
      }

      QDomNode node = submoleculeNode.firstChild();
      while (!node.isNull()) {
        if (node.nodeName() == "atom" and node.hasAttributes()) {
          readAtom(node, pSubmolecule);
        }
        if (node.nodeName() == "std" and node.hasAttributes()) {
          // todo
          // <std code="WT1" type="m" />

/*

    if (pSubMolecule->hasStdFrag()) {
      QDomElement stdFragElem = doc.createElement("std");
      subMolElem.appendChild(stdFragElem);

      stdFrag* pStdFrag = pSubMolecule->getStdFrag();
      stdFragElem.setAttribute("code", QString::fromStdString(pStdFrag->getSymbol()));
      stdFragElem.setAttribute("type", QString::fromStdString(pStdFrag->getType()));
    }

*/
        }
        node = node.nextSibling();
      }
    }
}

// ============================================================
// Function : readAtom
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void mtkppParser::readAtom(QDomNode atomNode, submolecule* pSubmolecule)
{
/*
    <atom hybridization="3" element="N" formalCharge="0" type="0" index="1" name=" N  " valence="0" >
     <coord x="40.506" y="-0.377" z="30.375" />
    </atom>
*/
    QDomElement atomElement = atomNode.toElement();
    if (!atomElement.isNull()) {
      atom* pAtom = pSubmolecule->addAtom();
      QString att;
      if (atomElement.hasAttribute("name")) {
        att = atomElement.attribute("name");
        pAtom->setName(att.toStdString());
      }

      if (atomElement.hasAttribute("index")) {
        att = atomElement.attribute("index");
        pAtom->setIndex(att.toInt());
      }

      if (atomElement.hasAttribute("element")) {
        att = atomElement.attribute("element");
        element* ele = pCollection->pElements->getElement(att.toStdString());
        if (ele) {
          pAtom->setElement(ele);
        }
      }

      if (atomElement.hasAttribute("valence")) {
        att = atomElement.attribute("valence");
        pAtom->setValence(att.toInt());
      }

      if (atomElement.hasAttribute("hybridization")) {
        att = atomElement.attribute("hybridization");
        pAtom->setHybridization(att.toInt());
      }

      if (atomElement.hasAttribute("type")) {
        att = atomElement.attribute("type");
        pAtom->setType(att.toInt());
      }

      if (atomElement.hasAttribute("formalCharge")) {
        att = atomElement.attribute("formalCharge");
        pAtom->setFormalCharge(att.toInt());
      }

      double xCoord = 0.0;
      double yCoord = 0.0;
      double zCoord = 0.0;
      QDomNode node = atomNode.firstChild();
      while (!node.isNull()) {
        if (node.nodeName() == "coord" and node.hasAttributes()) {
          QDomElement coordElement = node.toElement();
          if (!coordElement.isNull()) {

            if (coordElement.hasAttribute("x")) {
              att = coordElement.attribute("x");
              xCoord = att.toDouble();
            }
            if (coordElement.hasAttribute("y")) {
              att = coordElement.attribute("y");
              yCoord = att.toDouble();
            }
            if (coordElement.hasAttribute("z")) {
              att = coordElement.attribute("z");
              zCoord = att.toDouble();
            }
            pAtom->setCoords(xCoord, yCoord, zCoord);
          }
        }

        node = node.nextSibling();
      }
    }
}

// ============================================================
// Function : readBond
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void mtkppParser::readBond(QDomNode bondNode, molecule* pMolecule)
{
    // <bond at1="1" at2="2" kind="0" rotatable="1" stereo="0" type="1"  topology="0" />

    QDomElement bondElement = bondNode.toElement();
    if (!bondElement.isNull()) {
      QString att;
      int at1 = 0;
      int at2 = 0;
      int kind = 0;
      int rotatable = 0;
      int stereo = 0;
      int type = 0;
      int topology= 0;

      if (bondElement.hasAttribute("at1")) {
        att = bondElement.attribute("at1");
        at1 = att.toInt();
      }

      if (bondElement.hasAttribute("at2")) {
        att = bondElement.attribute("at2");
        at2 = att.toInt();
      }

      if (bondElement.hasAttribute("kind")) {
        att = bondElement.attribute("kind");
        kind = att.toInt();
      }

      if (bondElement.hasAttribute("rotatable")) {
        att = bondElement.attribute("rotatable");
        rotatable = att.toInt();
      }

      if (bondElement.hasAttribute("stereo")) {
        att = bondElement.attribute("stereo");
        stereo = att.toInt();
      }

      if (bondElement.hasAttribute("type")) {
        att = bondElement.attribute("type");
        type = att.toInt();
      }

      if (bondElement.hasAttribute("topology")) {
        att = bondElement.attribute("topology");
        topology = att.toInt();
      }

      atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
      atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);

      if (pAtom1 and pAtom2) {
        pMolecule->addBond(pAtom1, pAtom2, type, stereo, topology, 0.0);
      }
    }
}

// ============================================================
// Function : readAngle
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void mtkppParser::readAngle(QDomNode angleNode, molecule* pMolecule)
{
    // <angle at1="57" at2="58" at3="68" />
    QDomElement angleElement = angleNode.toElement();
    if (!angleElement.isNull()) {
      QString att;
      int at1 = 0;
      int at2 = 0;
      int at3 = 0;

      if (angleElement.hasAttribute("at1")) {
        att = angleElement.attribute("at1");
        at1 = att.toInt();
      }

      if (angleElement.hasAttribute("at2")) {
        att = angleElement.attribute("at2");
        at2 = att.toInt();
      }

      if (angleElement.hasAttribute("at3")) {
        att = angleElement.attribute("at3");
        at3 = att.toInt();
      }

      atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
      atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);
      atom* pAtom3 = pMolecule->getAtom(at3, 1, 0, 0);

      if (pAtom1 and pAtom2 and pAtom3) {
        pMolecule->addAngle(pAtom1, pAtom2, pAtom3, 0.0);
      }
    }
}

// ============================================================
// Function : readTorsion
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void mtkppParser::readTorsion(QDomNode torNode, molecule* pMolecule)
{
    // <torsion at1="48" at2="49" at3="74" at4="66" />
    QDomElement torElement = torNode.toElement();
    if (!torElement.isNull()) {
      QString att;
      int at1 = 0;
      int at2 = 0;
      int at3 = 0;
      int at4 = 0;

      if (torElement.hasAttribute("at1")) {
        att = torElement.attribute("at1");
        at1 = att.toInt();
      }

      if (torElement.hasAttribute("at2")) {
        att = torElement.attribute("at2");
        at2 = att.toInt();
      }

      if (torElement.hasAttribute("at3")) {
        att = torElement.attribute("at3");
        at3 = att.toInt();
      }

      if (torElement.hasAttribute("at4")) {
        att = torElement.attribute("at4");
        at4 = att.toInt();
      }

      atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
      atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);
      atom* pAtom3 = pMolecule->getAtom(at3, 1, 0, 0);
      atom* pAtom4 = pMolecule->getAtom(at4, 1, 0, 0);

      if (pAtom1 and pAtom2 and pAtom3 and pAtom4) {
        pMolecule->addTorsion(pAtom1, pAtom2, pAtom3, pAtom4, 0.0);
      }
    }
}

#endif // USE_QT

    // ---------------------------------- //
    // -     TINYXML READ FUNCTIONS     - //
    // ---------------------------------- //

#ifdef USE_TINYXML
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::Read(const std::string &fileName, collection* pCollection)
{
     std::cout << " IMPLEMENT mtkppParser::Read using tinyxml " << fileName << std::endl;
}

// ============================================================
// Function : readMolecule
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::readMolecule(TiXmlNode* molNode, molecule* pMolecule)
{
    pCollection = pMolecule->getParent();

    TiXmlElement* molElement = molNode->ToElement();

    if (molElement->Attribute("identity")) {
      pMolecule->setName(molElement->Attribute("identity"));
    }

    for (TiXmlNode* smolNode = molNode->FirstChild("submolecule"); smolNode; smolNode = smolNode->NextSibling("submolecule")) {
      readSubmolecule(smolNode, pMolecule);
    }

    for (TiXmlNode* bondNode = molNode->FirstChild("bond"); bondNode; bondNode = bondNode->NextSibling("bond")) {
      readBond(bondNode, pMolecule);
    }

    for (TiXmlNode* angleNode = molNode->FirstChild("angle"); angleNode; angleNode = angleNode->NextSibling("angle")) {
      readAngle(angleNode, pMolecule);
    }

    for (TiXmlNode* torNode = molNode->FirstChild("torsion"); torNode; torNode = torNode->NextSibling("torsion")) {
      readTorsion(torNode, pMolecule);
    }
}

// ============================================================
// Function : readSubmolecule
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::readSubmolecule(TiXmlNode* smolNode, molecule* pMolecule)
{
    TiXmlElement* smolElement = smolNode->ToElement();
    if (smolElement) {
      submolecule* pSubmolecule = pMolecule->addSubMolecule();

      if (smolElement->Attribute("name")) {
        pSubmolecule->setName(smolElement->Attribute("name"));
      }

      for (TiXmlNode* atomNode = smolNode->FirstChild("atom"); atomNode; atomNode = atomNode->NextSibling("atom")) {
        readAtom(atomNode, pSubmolecule);
      }

      // todo
      for (TiXmlNode* stdNode = smolNode->FirstChild("std"); stdNode; stdNode = stdNode->NextSibling("std")) {
        /*
        if (pSubMolecule->hasStdFrag()) {
          QDomElement stdFragElem = doc.createElement("std");
          subMolElem.appendChild(stdFragElem);

          stdFrag* pStdFrag = pSubMolecule->getStdFrag();
          stdFragElem.setAttribute("code", QString::fromStdString(pStdFrag->getSymbol()));
          stdFragElem.setAttribute("type", QString::fromStdString(pStdFrag->getType()));
        }
        */
      }
    }
}

// ============================================================
// Function : readAtom
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::readAtom(TiXmlNode* atomNode, submolecule* pSubmolecule)
{
    TiXmlElement* atomElement = atomNode->ToElement();
    if (atomElement) {
      atom* pAtom = pSubmolecule->addAtom();
      std::string att;

      if (atomElement->Attribute("name")) {
        pAtom->setName(atomElement->Attribute("name"));
      }

      if (atomElement->Attribute("index")) {
        pAtom->setIndex(string2Int(atomElement->Attribute("index")));
      }

      if (atomElement->Attribute("element")) {
        att = atomElement->Attribute("element");
        element* ele = pCollection->pElements->getElement(att);
        if (ele) {
          pAtom->setElement(ele);
        }
      }

      if (atomElement->Attribute("valence")) {
        pAtom->setValence(string2Int(atomElement->Attribute("valence")));
      }

      if (atomElement->Attribute("hybridization")) {
        pAtom->setHybridization(string2Int(atomElement->Attribute("hybridization")));
      }

      if (atomElement->Attribute("type")) {
        pAtom->setType(string2Int(atomElement->Attribute("type")));
      }

      if (atomElement->Attribute("formalCharge")) {
        pAtom->setFormalCharge(string2Int(atomElement->Attribute("formalCharge")));
      }

      double xCoord = 0.0;
      double yCoord = 0.0;
      double zCoord = 0.0;
      for (TiXmlNode* coordNode = atomNode->FirstChild("coord"); coordNode; coordNode = coordNode->NextSibling("coord")) {
        TiXmlElement* coordElement = coordNode->ToElement();
        if (coordElement->Attribute("x")) {
          xCoord = string2Double(coordElement->Attribute("x"));
        }

        if (coordElement->Attribute("y")) {
          yCoord = string2Double(coordElement->Attribute("y"));
        }

        if (coordElement->Attribute("z")) {
          zCoord = string2Double(coordElement->Attribute("z"));
        }
        pAtom->setCoords(xCoord, yCoord, zCoord);
      }
    }
}

// ============================================================
// Function : readBond
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::readBond(TiXmlNode* bondNode, molecule* pMolecule)
{
    TiXmlElement* bondElement = bondNode->ToElement();
    if (bondElement) {

      std::string att;
      int at1 = 0;
      int at2 = 0;
      int kind = 0;
      int rotatable = 0;
      int stereo = 0;
      int type = 0;
      int topology= 0;

      if (bondElement->Attribute("at1")) {
        at1 = string2Int(bondElement->Attribute("at1"));
      }

      if (bondElement->Attribute("at2")) {
        at2 = string2Int(bondElement->Attribute("at2"));
      }

      if (bondElement->Attribute("kind")) {
        kind = string2Int(bondElement->Attribute("kind"));
      }

      if (bondElement->Attribute("rotatable")) {
        rotatable = string2Int(bondElement->Attribute("rotatable"));
      }

      if (bondElement->Attribute("stereo")) {
        stereo = string2Int(bondElement->Attribute("stereo"));
      }

      if (bondElement->Attribute("type")) {
        type = string2Int(bondElement->Attribute("type"));
      }

      if (bondElement->Attribute("topology")) {
        topology = string2Int(bondElement->Attribute("topology"));
      }

      atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
      atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);

      if (pAtom1 and pAtom2) {
        pMolecule->addBond(pAtom1, pAtom2, type, stereo, topology, 0.0);
      }
    }
}

// ============================================================
// Function : readAngle
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::readAngle(TiXmlNode* angleNode, molecule* pMolecule)
{
    TiXmlElement* angleElement = angleNode->ToElement();
    if (angleElement) {
      std::string att;
      int at1 = 0;
      int at2 = 0;
      int at3 = 0;

      if (angleElement->Attribute("at1")) {
        at1 = string2Int(angleElement->Attribute("at1"));
      }

      if (angleElement->Attribute("at2")) {
        at2 = string2Int(angleElement->Attribute("at2"));
      }

      if (angleElement->Attribute("at3")) {
        at3 = string2Int(angleElement->Attribute("at3"));
      }

      atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
      atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);
      atom* pAtom3 = pMolecule->getAtom(at3, 1, 0, 0);

      if (pAtom1 and pAtom2 and pAtom3) {
        pMolecule->addAngle(pAtom1, pAtom2, pAtom3, 0.0);
      }
    }
}

// ============================================================
// Function : readTorsion
// ------------------------------------------------------------
// Parses MTK++ State xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::readTorsion(TiXmlNode* torsionNode, molecule* pMolecule)
{
    TiXmlElement* torsionElement = torsionNode->ToElement();
    if (torsionElement) {
      std::string att;
      int at1 = 0;
      int at2 = 0;
      int at3 = 0;
      int at4 = 0;

      if (torsionElement->Attribute("at1")) {
        at1 = string2Int(torsionElement->Attribute("at1"));
      }

      if (torsionElement->Attribute("at2")) {
        at2 = string2Int(torsionElement->Attribute("at2"));
      }

      if (torsionElement->Attribute("at3")) {
        at3 = string2Int(torsionElement->Attribute("at3"));
      }

      if (torsionElement->Attribute("at4")) {
        at4 = string2Int(torsionElement->Attribute("at4"));
      }

      atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
      atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);
      atom* pAtom3 = pMolecule->getAtom(at3, 1, 0, 0);
      atom* pAtom4 = pMolecule->getAtom(at4, 1, 0, 0);

      if (pAtom1 and pAtom2 and pAtom3 and pAtom4) {
        pMolecule->addTorsion(pAtom1, pAtom2, pAtom3, pAtom4, 0.0);
      }
    }
}

#endif // USE_TINYXML

    // ---------------------------------- //
    // -    XERCES-C READ FUNCTIONS     - //
    // ---------------------------------- //

#ifdef USE_XERCES
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::Read(const std::string &fileName, collection* pCollection)
{
    std::cout << " IMPLEMENT mtkppParser::Read using xerces-c" << std::endl;
}

// ============================================================
// Function : readMolecule
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::readMolecule(DOMNode* e, molecule* pMolecule)
{
    pCollection = pMolecule->getParent();

    if (e->hasAttributes()) {
      DOMNamedNodeMap *pAttributes = e->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());

        if (std::string(name) == std::string("identity")) {
          char *cMolName = XMLString::transcode(pAttributeNode->getNodeValue());
          pMolecule->setName(std::string(cMolName));
          delete cMolName;
        }

      }
    }

    DOMNode *a;
    for (a = e->getFirstChild(); a != 0; a = a->getNextSibling()) {
      if (a->getNodeType() == 1) {
        if ((std::string)(XC(a->getNodeName())) == "submolecule") {
          readSubmolecule(a, pMolecule);
        }

        if ((std::string)(XC(a->getNodeName())) == "bond") {
          readBond(a, pMolecule);
        }

        if ((std::string)(XC(a->getNodeName())) == "angle") {
          readAngle(a, pMolecule);
        }

        if ((std::string)(XC(a->getNodeName())) == "torsion") {
          readTorsion(a, pMolecule);
        }
      }
    }
}

// ============================================================
// Function : readSubmolecule
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::readSubmolecule(DOMNode* e, molecule* pMolecule)
{
    submolecule* pSubmolecule = pMolecule->addSubMolecule();

    if (e->hasAttributes()) {
      DOMNamedNodeMap *pAttributes = e->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());

        if (std::string(name) == std::string("name")) {
          char *cSMolName = XMLString::transcode(pAttributeNode->getNodeValue());
          pSubmolecule->setName(std::string(cSMolName));
          delete cSMolName;
        }

      }
    }

    DOMNode *a;
    for (a = e->getFirstChild(); a != 0; a = a->getNextSibling()) {
      if (a->getNodeType() == 1) {
        if ((std::string)(XC(a->getNodeName())) == "atom") {
          readAtom(a, pSubmolecule);
        }
      }
    }
}

// ============================================================
// Function : readAtom
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::readAtom(DOMNode* e, submolecule* pSubmolecule)
{
    atom* pAtom = pSubmolecule->addAtom();

    if (e->hasAttributes()) {
      DOMNamedNodeMap *pAttributes = e->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());

        if (std::string(name) == std::string("name")) {
          char *cName = XMLString::transcode(pAttributeNode->getNodeValue());
          pAtom->setName(std::string(cName));
          delete cName;
        }

        if (std::string(name) == std::string("index")) {
          char *cIndex = XMLString::transcode(pAttributeNode->getNodeValue());
          pAtom->setIndex(atoi(cIndex));
          delete cIndex;
        }

        if (std::string(name) == std::string("element")) {
          char *cEle = XMLString::transcode(pAttributeNode->getNodeValue());
          element* ele = pCollection->pElements->getElement(std::string(cEle));
          if (ele) {
            pAtom->setElement(ele);
          }
          delete cEle;
        }

        if (std::string(name) == std::string("valence")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          pAtom->setValence(atoi(cI));
          delete cI;
        }

        if (std::string(name) == std::string("hybridization")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          pAtom->setHybridization(atoi(cI));
          delete cI;
        }

        if (std::string(name) == std::string("type")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          pAtom->setType(atoi(cI));
          delete cI;
        }

        if (std::string(name) == std::string("formalCharge")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          pAtom->setFormalCharge(atoi(cI));
          delete cI;
        }
      }
    }

    DOMNode *a;
    double xCoord = 0.0;
    double yCoord = 0.0;
    double zCoord = 0.0;
    for (a = e->getFirstChild(); a != 0; a = a->getNextSibling()) {
      if (a->getNodeType() == 1) {
        if ((std::string)(XC(a->getNodeName())) == "coord") {
          if (a->hasAttributes()) {
            DOMNamedNodeMap *pAttributes = a->getAttributes();
            int nSize = pAttributes->getLength();

            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              char *name = XMLString::transcode(pAttributeNode->getName());

              if (std::string(name) == std::string("x")) {
                xCoord = strtod(XC(pAttributeNode->getNodeValue()), 0);
              }

              if (std::string(name) == std::string("y")) {
                yCoord = strtod(XC(pAttributeNode->getNodeValue()), 0);
              }

              if (std::string(name) == std::string("z")) {
                zCoord = strtod(XC(pAttributeNode->getNodeValue()), 0);
              }
              pAtom->setCoords(xCoord, yCoord, zCoord);
            }
          }
        }
      }
    }
}
// ============================================================
// Function : readBond
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::readBond(DOMNode* e, molecule* pMolecule)
{
    int at1 = 0;
    int at2 = 0;
    int kind = 0;
    int rotatable = 0;
    int stereo = 0;
    int type = 0;
    int topology= 0;

    if (e->hasAttributes()) {
      DOMNamedNodeMap *pAttributes = e->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());

        if (std::string(name) == std::string("at1")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at1 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("at2")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at2 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("kind")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          kind = atoi(cI);
          delete cI;
        }

        if (std::string(name) == std::string("rotatable")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          rotatable = atoi(cI);
          delete cI;
        }

        if (std::string(name) == std::string("stereo")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          stereo = atoi(cI);
          delete cI;
        }

        if (std::string(name) == std::string("type")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          type = atoi(cI);
          delete cI;
        }

        if (std::string(name) == std::string("topology")) {
          char *cI = XMLString::transcode(pAttributeNode->getNodeValue());
          topology = atoi(cI);
          delete cI;
        }
      }
    }

    atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
    atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);

    if (pAtom1 and pAtom2) {
      pMolecule->addBond(pAtom1, pAtom2, type, stereo, topology, 0.0);
    }
}

// ============================================================
// Function : readAngle
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::readAngle(DOMNode* e, molecule* pMolecule)
{
    int at1 = 0;
    int at2 = 0;
    int at3 = 0;

    if (e->hasAttributes()) {
      DOMNamedNodeMap *pAttributes = e->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());

        if (std::string(name) == std::string("at1")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at1 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("at2")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at2 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("at3")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at3 = atoi(cS);
          delete cS;
        }
      }
    }

    atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
    atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);
    atom* pAtom3 = pMolecule->getAtom(at3, 1, 0, 0);

    if (pAtom1 and pAtom2 and pAtom3) {
      pMolecule->addAngle(pAtom1, pAtom2, pAtom3, 0.0);
    }
}

// ============================================================
// Function : readTorsion
// ------------------------------------------------------------
// Parses MTK++ State xml files using xerces-c
// ------------------------------------------------------------
void mtkppParser::readTorsion(DOMNode* e, molecule* pMolecule)
{
    int at1 = 0;
    int at2 = 0;
    int at3 = 0;
    int at4 = 0;

    if (e->hasAttributes()) {
      DOMNamedNodeMap *pAttributes = e->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());

        if (std::string(name) == std::string("at1")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at1 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("at2")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at2 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("at3")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at3 = atoi(cS);
          delete cS;
        }

        if (std::string(name) == std::string("at4")) {
          char *cS = XMLString::transcode(pAttributeNode->getNodeValue());
          at4 = atoi(cS);
          delete cS;
        }
      }
    }

    atom* pAtom1 = pMolecule->getAtom(at1, 1, 0, 0);
    atom* pAtom2 = pMolecule->getAtom(at2, 1, 0, 0);
    atom* pAtom3 = pMolecule->getAtom(at3, 1, 0, 0);
    atom* pAtom4 = pMolecule->getAtom(at4, 1, 0, 0);

    if (pAtom1 and pAtom2 and pAtom3 and pAtom4) {
      pMolecule->addTorsion(pAtom1, pAtom2, pAtom3, pAtom4, 0.0);
    }
}

#endif // USE_XERCES

#ifdef USE_QT
    // ----------------------------- //
    // -    Qt WRITE FUNCTIONS     - //
    // ----------------------------- //

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using Qt
// ------------------------------------------------------------
void mtkppParser::Write(const std::string &fileName, collection* pCollection)
{
    std::cout << " IMPLEMENT mtkppParser::Write using Qt" << std::endl;
}

// ============================================================
// Function : writeMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeMolecule(QDomDocument doc, molecule* pMolecule)
{
    QDomElement rootElem = doc.documentElement();
    QDomElement moleculeElem = doc.createElement("molecule");
    rootElem.appendChild(moleculeElem);
    moleculeElem.setAttribute("identity", QString::fromStdString(pMolecule->getName()));

    // submolecule's
    std::vector<submolecule*> smolList = pMolecule->getSubMoleculeList();
    for (unsigned int j = 0; j < smolList.size(); j++) {
      this->writeSubMolecule(doc, moleculeElem, smolList[j]);
    }

    // Bonds
    std::map<int, Bond*> bondMap = pMolecule->getBondMap();
    this->writeBonds(doc, moleculeElem, bondMap);

    // Angles
    std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();
    this->writeAngles(doc, moleculeElem, angleMap);

    // Torsions
    std::map<ULONG_KIND, Torsion*> torsionMap = pMolecule->getTorsionMap();
    this->writeTorsions(doc, moleculeElem, torsionMap);

    // Impropers
    std::map<int, Improper*> improperMap = pMolecule->getImproperMap();
    for (ImproperMapIterator i = improperMap.begin(); i != improperMap.end(); i++) {
      pImproper = i->second;
    }
}

// ============================================================
// Function : writeMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeMolecule(QDomDocument doc, QDomElement groupElem, molecule* pMolecule)
{
    QDomElement moleculeElem = doc.createElement("molecule");
    groupElem.appendChild(moleculeElem);
    moleculeElem.setAttribute("identity", QString::fromStdString(pMolecule->getName()));

    // submolecule's
    std::vector<submolecule*> smolList = pMolecule->getSubMoleculeList();
    for (unsigned int j = 0; j < smolList.size(); j++) {
      this->writeSubMolecule(doc, moleculeElem, smolList[j]);
    }

    // Bonds
    std::map<int, Bond*> bondMap = pMolecule->getBondMap();
    this->writeBonds(doc, moleculeElem, bondMap);

    // Angles
    std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();
    this->writeAngles(doc, moleculeElem, angleMap);

    // Torsions
    std::map<ULONG_KIND, Torsion*> torsionMap = pMolecule->getTorsionMap();
    this->writeTorsions(doc, moleculeElem, torsionMap);

    // Impropers
    std::map<int, Improper*> improperMap = pMolecule->getImproperMap();
    for (ImproperMapIterator i = improperMap.begin(); i != improperMap.end(); i++) {
      pImproper = i->second;
    }
}

// ============================================================
// Function : writeSubMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeSubMolecule(QDomDocument doc, QDomElement moleculeElem, submolecule* pSubMolecule)
{
    QDomElement subMolElem = doc.createElement("submolecule");
    moleculeElem.appendChild(subMolElem);
    subMolElem.setAttribute("name", QString::fromStdString(pSubMolecule->getName()));
    //subMolElem.setAttribute("index", QString::number(pSubMolecule->getSubMolId(), 10));
    //subMolElem.setAttribute("fileID", QString::fromStdString(pSubMolecule->getiCode()));

    // atom's
    std::vector<atom*> atomList = pSubMolecule->getAtomList();
    for (unsigned int k = 0; k < atomList.size(); k++) {
      this->writeAtom(doc, subMolElem, pSubMolecule, atomList[k]);
    }

    // stdFrag
    if (pSubMolecule->hasStdFrag()) {
      QDomElement stdFragElem = doc.createElement("std");
      subMolElem.appendChild(stdFragElem);

      stdFrag* pStdFrag = pSubMolecule->getStdFrag();
      stdFragElem.setAttribute("code", QString::fromStdString(pStdFrag->getSymbol()));
      stdFragElem.setAttribute("type", QString::fromStdString(pStdFrag->getType()));
    }
}

// ============================================================
// Function : writeAtom
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeAtom(QDomDocument doc, QDomElement smolElem, submolecule* pSubMolecule, atom* pAtom)
{
    QDomElement atomElem = doc.createElement("atom");
    smolElem.appendChild(atomElem);

    atomElem.setAttribute("element", QString::fromStdString(pAtom->getElementSymbol()));
    atomElem.setAttribute("name", QString::fromStdString(pAtom->getName()));

    atomElem.setAttribute("index", QString::number(pAtom->getIndex(), 10));
    //atomElem.setAttribute("colIndex", QString::number(pAtom->getColIndex(), 10));
    //atomElem.setAttribute("fileID", QString::number(pAtom->getFileID(), 10));
    atomElem.setAttribute("valence", QString::number(pAtom->getValence(), 10));
    atomElem.setAttribute("hybridization", QString::number(pAtom->getHybridization(), 10));
    atomElem.setAttribute("type", QString::number(pAtom->getType(), 10));
    atomElem.setAttribute("formalCharge", QString::number(pAtom->getFormalCharge(), 10));

    // Coordinates
    QDomElement coordElem = doc.createElement("coord");
    atomElem.appendChild(coordElem);
    coordElem.setAttribute("x", QString::number(pAtom->getX()));
    coordElem.setAttribute("y", QString::number(pAtom->getY()));
    coordElem.setAttribute("z", QString::number(pAtom->getZ()));

    // Properties
    std::map<std::string, double> atomDMap = pAtom->getPropertyMap();
    std::map<std::string, int> atomIMap = pAtom->getIntPropertyMap();

    for (PropertyMapIterator p = atomDMap.begin(); p != atomDMap.end(); p++ ) {
      QDomElement propElem = doc.createElement("property");
      atomElem.appendChild(propElem);
      propElem.setAttribute("name", QString::fromStdString(p->first));
      propElem.setAttribute("value", QString::number(p->second));
      propElem.setAttribute("type", "d");
    }

    for (intPropertyMapIterator p = atomIMap.begin(); p != atomIMap.end(); p++ ) {
      QDomElement propElem = doc.createElement("property");
      atomElem.appendChild(propElem);
      propElem.setAttribute("name", QString::fromStdString(p->first));
      propElem.setAttribute("value", QString::number(p->second, 10));
      propElem.setAttribute("type", "i");
    }
}

// ============================================================
// Function : writeBonds
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeBonds(QDomDocument doc, QDomElement moleculeElem, std::map<int, Bond*> bondMap)
{
    for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
      pBond = b->second;
      QDomElement bondElem = doc.createElement("bond");
      moleculeElem.appendChild(bondElem);

      //bondElem.setAttribute("index", QString::number(b->first, 10));
      bondElem.setAttribute("at1", QString::number(pBond->atom1->getIndex(), 10));
      bondElem.setAttribute("at2", QString::number(pBond->atom2->getIndex(), 10));
      bondElem.setAttribute("type", QString::number(pBond->type, 10));
      bondElem.setAttribute("stereo", QString::number(pBond->stereo, 10));
      bondElem.setAttribute("topology", QString::number(pBond->topology, 10));
      bondElem.setAttribute("kind", QString::number(pBond->kind, 10));
      bondElem.setAttribute("rotatable", QString::number(pBond->rotatable, 10));
    }
}

// ============================================================
// Function : writeAngles
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeAngles(QDomDocument doc, QDomElement moleculeElem, std::map<ULONG_KIND, Angle*> angleMap)
{
    for (AngleMapIterator a = angleMap.begin(); a != angleMap.end(); a++) {
      pAngle = a->second;
      QDomElement angleElem = doc.createElement("angle");
      moleculeElem.appendChild(angleElem);
      //angleElem.setAttribute("index", QString::number(a->first, 10));

      angleElem.setAttribute("index", string2QString(uLongKind2String(a->first)));

      angleElem.setAttribute("at1", QString::number(pAngle->atom1->getIndex(), 10));
      angleElem.setAttribute("at2", QString::number(pAngle->atom2->getIndex(), 10));
      angleElem.setAttribute("at3", QString::number(pAngle->atom3->getIndex(), 10));
    }
}

// ============================================================
// Function : writeTorsions
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeTorsions(QDomDocument doc, QDomElement moleculeElem, std::map<ULONG_KIND, Torsion*> torsionMap)
{
    for (TorsionMapIterator t = torsionMap.begin(); t != torsionMap.end(); t++) {
      pTorsion = t->second;
      QDomElement torElem = doc.createElement("torsion");
      moleculeElem.appendChild(torElem);
      //torElem.setAttribute("index", QString::number(t->first, 10)); //uLongKind2String
      torElem.setAttribute("at1", QString::number(pTorsion->atom1->getIndex(), 10));
      torElem.setAttribute("at2", QString::number(pTorsion->atom2->getIndex(), 10));
      torElem.setAttribute("at3", QString::number(pTorsion->atom3->getIndex(), 10));
      torElem.setAttribute("at4", QString::number(pTorsion->atom4->getIndex(), 10));
    }
}

#endif // USE_QT

#ifdef USE_XERCES
    // ----------------------------------- //
    // -    XERCES-C WRITE FUNCTIONS     - //
    // ----------------------------------- //

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using xercesc
// ------------------------------------------------------------
void mtkppParser::Write(const std::string &fileName, collection* pCollection)
{
    char* xmlFile = (char*)(fileName.c_str());

    try {
      XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
      char* message = XMLString::transcode(toCatch.getMessage());
      std::cout << "Error during initialization! :\n"
                << message << "\n";
      XMLString::release(&message);
      setError(1);
      return;
    }

    int errorCode = 0;

    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);

    // xml settings
    bool gDoNamespaces          = false;
    bool gDoSchema              = false;
    bool gSchemaFullChecking    = false;
    bool gDoCreate              = false;

#if (XERCES_VERSION_MAJOR == 3)
#else
    bool gSplitCdataSections    = true;
    bool gDiscardDefaultContent = true;
    bool gFormatPrettyPrint     = true;
    bool gWriteBOM              = false;
#endif

    XercesDOMParser::ValSchemes gValScheme  = XercesDOMParser::Val_Auto;

    //  Create our parser, then attach an error handler to the parser.
    //  The parser will call back to methods of the ErrorHandler if it
    //  discovers errors during the course of parsing the XML document.
    XercesDOMParser *parser = new XercesDOMParser;
    parser->setValidationScheme(gValScheme);
    parser->setDoNamespaces(gDoNamespaces);
    parser->setDoSchema(gDoSchema);
    parser->setValidationSchemaFullChecking(gSchemaFullChecking);
    parser->setCreateEntityReferenceNodes(gDoCreate);

    //DOMTreeErrorReporter *errReporter = new DOMTreeErrorReporter();
    //parser->setErrorHandler(errReporter);

    {
      if (impl != 0) {
        try {
#if (XERCES_VERSION_MAJOR == 3)
          DOMLSSerializer *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer(XMLPlatformUtils::fgMemoryManager);
#else
          DOMWriter *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
#endif

          XERCES_CPP_NAMESPACE::DOMDocument* doc = impl->createDocument(
            0,           // root element namespace URI.
            X("MTKpp"), // root element name
            0);          // document type object (DTD).

          // Elements

          // Libraries

          // Parameters

          // Molecules
          std::vector<molecule*> molList = pCollection->getMoleculeList();
          for (unsigned int i = 0; i < molList.size(); i++) {
            this->writeMolecule(doc, molList[i]);
          }

          // Metal Centers
          std::vector<metalCenter*> metalCenterList = pCollection->getMetalCenters();
          for (unsigned int i = 0; i < metalCenterList.size(); i++) {
            this->writeMetalCenter(doc, metalCenterList[i]);
          }

#if (XERCES_VERSION_MAJOR == 3)
#else
          // set feature if the serializer supports the feature/mode
          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections))
            theSerializer->setFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent))
            theSerializer->setFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint))
            theSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTBOM, gWriteBOM))
            theSerializer->setFeature(XMLUni::fgDOMWRTBOM, gWriteBOM);
#endif

          XMLFormatTarget *myFormTarget;
          myFormTarget = new LocalFileFormatTarget(xmlFile);

#if (XERCES_VERSION_MAJOR == 3)
          DOMLSOutput         *lsOutput = ((DOMImplementationLS*)impl)->createLSOutput(XMLPlatformUtils::fgMemoryManager);
          lsOutput->setByteStream(myFormTarget);

          // do the serialization to xml file.
          theSerializer->write(doc, lsOutput);
#else
          // do the serialization through DOMWriter::writeNode();
          theSerializer->writeNode(myFormTarget, *doc);
#endif

          delete theSerializer;
          delete myFormTarget;
          //delete myErrorHandler;

          doc->release();
        }
        catch (const OutOfMemoryException&) {
          XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
          errorCode = 5;
        }
        catch (const DOMException& e) {
          XERCES_STD_QUALIFIER cerr << "DOMException code is:  " << e.code << XERCES_STD_QUALIFIER endl;
          errorCode = 2;
        }
        catch (...) {
          XERCES_STD_QUALIFIER cerr << "An error occurred creating the document" << XERCES_STD_QUALIFIER endl;
          errorCode = 3;
        }
      }
      else {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorCode = 4;
      }
    }

    XMLPlatformUtils::Terminate();
    setError(errorCode);
    return;
}

// ============================================================
// Function : writeMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeMolecule(XERCES_CPP_NAMESPACE::DOMDocument* doc, molecule* pMolecule)
{
    //std::cout << " mtkppParser::writeMolecule (Xerces-c) " << std::endl;
    DOMElement* rootElem = doc->getDocumentElement();
    DOMElement* moleculeElem = doc->createElement(X("molecule"));
    rootElem->appendChild(moleculeElem);
    moleculeElem->setAttribute(X("identity"), X(pMolecule->getName().c_str()));

    // submolecule's
    std::vector<submolecule*> smolList = pMolecule->getSubMoleculeList();
    for (unsigned int j = 0; j < smolList.size(); j++) {
      this->writeSubMolecule(doc, moleculeElem, smolList[j]);
    }

    // Bonds
    std::map<int, Bond*> bondMap = pMolecule->getBondMap();
    this->writeBonds(doc, moleculeElem, bondMap);

    // Angles
    std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();
    this->writeAngles(doc, moleculeElem, angleMap);

    // Torsions
    std::map<ULONG_KIND, Torsion*> torsionMap = pMolecule->getTorsionMap();
    this->writeTorsions(doc, moleculeElem, torsionMap);

    // Impropers
    std::map<int, Improper*> improperMap = pMolecule->getImproperMap();
    for (ImproperMapIterator i = improperMap.begin(); i != improperMap.end(); i++) {
      pImproper = i->second;
    }
}

// ============================================================
// Function : writeSubMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeSubMolecule(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, submolecule* pSubMolecule)
{
    DOMElement* smolElem = doc->createElement(X("submolecule"));
    moleculeElem->appendChild(smolElem);

    smolElem->setAttribute(X("name"), X(pSubMolecule->getName().c_str()));
    smolElem->setAttribute(X("index"), X(int2String(pSubMolecule->getSubMolId()).c_str()));
    smolElem->setAttribute(X("fileID"), X(pSubMolecule->getiCode().c_str()));

    // atom's
    std::vector<atom*> atomList = pSubMolecule->getAtomList();
    for (unsigned int k = 0; k < atomList.size(); k++) {
      this->writeAtom(doc, smolElem, pSubMolecule, atomList[k]);
    }
}

// ============================================================
// Function : writeAtom
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeAtom(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement*  smolElem, submolecule* pSubMolecule, atom* pAtom)
{
    DOMElement* atomElem = doc->createElement(X("atom"));
    smolElem->appendChild(atomElem);

    atomElem->setAttribute(X("element"), X(pAtom->getElementSymbol().c_str()));
    atomElem->setAttribute(X("name"), X(pAtom->getName().c_str()));
    atomElem->setAttribute(X("index"), X(int2String(pAtom->getIndex()).c_str()));
    atomElem->setAttribute(X("colIndex"), X(int2String(pAtom->getColIndex()).c_str()));
    atomElem->setAttribute(X("fileID"), X(int2String(pAtom->getFileID()).c_str()));
    atomElem->setAttribute(X("valence"), X(int2String(pAtom->getValence()).c_str()));
    atomElem->setAttribute(X("hybridization"), X(int2String(pAtom->getHybridization()).c_str()));
    atomElem->setAttribute(X("type"), X(int2String(pAtom->getType()).c_str()));
    atomElem->setAttribute(X("formalCharge"), X(int2String(pAtom->getFormalCharge()).c_str()));

    std::map<std::string, double> atomDMap = pAtom->getPropertyMap();
    std::map<std::string, int> atomIMap = pAtom->getIntPropertyMap();

    for (PropertyMapIterator p = atomDMap.begin(); p != atomDMap.end(); p++ ) {
      std::cout << "dProp: '" << p->first << "', Value: "
                << p->second << std::endl;

      DOMElement* propElem = doc->createElement(X("property"));
      atomElem->appendChild(propElem);
      propElem->setAttribute(X("name"), X(p->first.c_str()));
      propElem->setAttribute(X("value"), X(double2String(p->second).c_str()));
      propElem->setAttribute(X("type"), X("d"));
    }

    for (intPropertyMapIterator p = atomIMap.begin(); p != atomIMap.end(); p++ ) {
      std::cout << "iProp: '" << p->first << "', Value: "
                << p->second << std::endl;

      DOMElement* propElem = doc->createElement(X("property"));
      atomElem->appendChild(propElem);
      propElem->setAttribute(X("name"), X(p->first.c_str()));
      propElem->setAttribute(X("value"), X(int2String(p->second).c_str()));
      propElem->setAttribute(X("type"), X("i"));
    }
}

// ============================================================
// Function : writeMetalCenter
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeMetalCenter(XERCES_CPP_NAMESPACE::DOMDocument* doc, metalCenter* pMetalCenter)
{
    std::cout << " mtkppParser::writeMetalCenter " << std::endl;
    DOMElement* rootElem = doc->getDocumentElement();
    DOMElement* mcElem = doc->createElement(X("metalCenter"));
    rootElem->appendChild(mcElem);

    mcElem->setAttribute(X("identity"), X(pMetalCenter->getName().c_str()));
}

// ============================================================
// Function : writeBonds
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeBonds(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, std::map<int, Bond*> bondMap)
{
    for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
      pBond = b->second;
      DOMElement* bondElem = doc->createElement(X("bond"));
      moleculeElem->appendChild(bondElem);
      bondElem->setAttribute(X("index"), X(int2String(b->first).c_str()));
      bondElem->setAttribute(X("at1"), X(int2String(pBond->atom1->getIndex()).c_str()));
      bondElem->setAttribute(X("at2"), X(int2String(pBond->atom2->getIndex()).c_str()));
      bondElem->setAttribute(X("type"), X(int2String(pBond->type).c_str()));
      bondElem->setAttribute(X("stereo"), X(int2String(pBond->stereo).c_str()));
      bondElem->setAttribute(X("topology"), X(int2String(pBond->topology).c_str()));
      bondElem->setAttribute(X("kind"), X(int2String(pBond->kind).c_str()));
      bondElem->setAttribute(X("rotatable"), X(int2String(pBond->rotatable).c_str()));
/*
    bondParam* pBondParam;
*/
    }
}

// ============================================================
// Function : writeAngles
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeAngles(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, std::map<ULONG_KIND, Angle*> angleMap)
{
    for (AngleMapIterator a = angleMap.begin(); a != angleMap.end(); a++) {
      pAngle = a->second;
      DOMElement* angleElem = doc->createElement(X("angle"));
      moleculeElem->appendChild(angleElem);
      angleElem->setAttribute(X("index"), X(uLongKind2String(a->first).c_str()));
      angleElem->setAttribute(X("at1"), X(int2String(pAngle->atom1->getIndex()).c_str()));
      angleElem->setAttribute(X("at2"), X(int2String(pAngle->atom2->getIndex()).c_str()));
      angleElem->setAttribute(X("at3"), X(int2String(pAngle->atom3->getIndex()).c_str()));
/*
    angleParam* pAngleParam;
*/
    }
}

// ============================================================
// Function : writeTorsions
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeTorsions(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* moleculeElem, std::map<ULONG_KIND, Torsion*> torsionMap)
{
    for (TorsionMapIterator t = torsionMap.begin(); t != torsionMap.end(); t++) {
      pTorsion = t->second;
      DOMElement* torElem = doc->createElement(X("torsion"));
      moleculeElem->appendChild(torElem);
      torElem->setAttribute(X("index"), X(uLongKind2String(t->first).c_str()));
      torElem->setAttribute(X("at1"), X(int2String(pTorsion->atom1->getIndex()).c_str()));
      torElem->setAttribute(X("at2"), X(int2String(pTorsion->atom2->getIndex()).c_str()));
      torElem->setAttribute(X("at3"), X(int2String(pTorsion->atom3->getIndex()).c_str()));
      torElem->setAttribute(X("at4"), X(int2String(pTorsion->atom4->getIndex()).c_str()));
    }
/*
    torsionParam* pTorsionParam;
*/
}
#endif // USE_XERCES


#ifdef USE_TINYXML
    // ----------------------------------- //
    // -    XERCES-C WRITE FUNCTIONS     - //
    // ----------------------------------- //

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using tinyxml
// ------------------------------------------------------------
void mtkppParser::Write(const std::string &fileName, collection* pCollection)
{
    //std::string errMessage = " Writing " + fileName;
    //errorLogger.throwError("mtkppParser", errMessage, INFO);

    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

    TiXmlElement* root = new TiXmlElement( "MTKpp" );
    doc.LinkEndChild(root);

    // Elements

    // Libraries

    // Parameters

    // Molecules
    std::vector<molecule*> molList = pCollection->getMoleculeList();
    for (unsigned int i = 0; i < molList.size(); i++) {
      this->writeMolecule(root, molList[i]);
    }

    // Metal Centers
    std::vector<metalCenter*> metalCenterList = pCollection->getMetalCenters();
    for (unsigned int i = 0; i < metalCenterList.size(); i++) {
      this->writeMetalCenter(root, metalCenterList[i]);
    }

    doc.SaveFile(fileName);

    //dump_to_stdout(&doc);

    return;
}

// ============================================================
// Function : writeMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeMolecule(TiXmlElement* root, molecule* pMolecule)
{
    TiXmlElement* molElem = new TiXmlElement("molecule");
    root->LinkEndChild(molElem);
    molElem->SetAttribute("identity", pMolecule->getName());

    // submolecule's
    std::vector<submolecule*> smolList = pMolecule->getSubMoleculeList();
    for (unsigned int j = 0; j < smolList.size(); j++) {
      this->writeSubMolecule(molElem, smolList[j]);
    }

    // Bonds
    std::map<int, Bond*> bondMap = pMolecule->getBondMap();
    this->writeBonds(molElem, bondMap);

    // Angles
    std::map<ULONG_KIND, Angle*> angleMap = pMolecule->getAngleMap();
    this->writeAngles(molElem, angleMap);

    // Torsions
    std::map<ULONG_KIND, Torsion*> torsionMap = pMolecule->getTorsionMap();
    this->writeTorsions(molElem, torsionMap);

    // Impropers
    std::map<int, Improper*> improperMap = pMolecule->getImproperMap();
    for (ImproperMapIterator i = improperMap.begin(); i != improperMap.end(); i++) {
      pImproper = i->second;
    }
}

// ============================================================
// Function : writeSubMolecule
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeSubMolecule(TiXmlElement* root, submolecule* pSubMolecule)
{
    TiXmlElement* smolElem = new TiXmlElement("submolecule");
    root->LinkEndChild(smolElem);

    smolElem->SetAttribute("name", pSubMolecule->getName());
    smolElem->SetAttribute("index", int2String(pSubMolecule->getSubMolId()));
    smolElem->SetAttribute("fileID", pSubMolecule->getiCode());

    // atom's
    std::vector<atom*> atomList = pSubMolecule->getAtomList();
    for (unsigned int k = 0; k < atomList.size(); k++) {
      this->writeAtom(smolElem, pSubMolecule, atomList[k]);
    }
}

// ============================================================
// Function : writeAtom
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeAtom(TiXmlElement* root, submolecule* smol, atom* pAtom)
{
    TiXmlElement* atomElem = new TiXmlElement("atom");
    root->LinkEndChild(atomElem);

    atomElem->SetAttribute("element", pAtom->getElementSymbol());
    atomElem->SetAttribute("name", pAtom->getName());
    atomElem->SetAttribute("index", int2String(pAtom->getIndex()));
    atomElem->SetAttribute("colIndex", int2String(pAtom->getColIndex()));
    atomElem->SetAttribute("fileID", int2String(pAtom->getFileID()));
    atomElem->SetAttribute("valence", int2String(pAtom->getValence()));
    atomElem->SetAttribute("hybridization", int2String(pAtom->getHybridization()));
    atomElem->SetAttribute("type", int2String(pAtom->getType()));
    atomElem->SetAttribute("formalCharge", int2String(pAtom->getFormalCharge()));

    // Coordinates
    TiXmlElement* coordElem = new TiXmlElement("coord");
    atomElem->LinkEndChild(coordElem);

    coordElem->SetAttribute("x", double2String(pAtom->getX()));
    coordElem->SetAttribute("y", double2String(pAtom->getY()));
    coordElem->SetAttribute("z", double2String(pAtom->getZ()));

    // Properties
    std::map<std::string, double> atomDMap = pAtom->getPropertyMap();
    std::map<std::string, int> atomIMap = pAtom->getIntPropertyMap();

    for (PropertyMapIterator p = atomDMap.begin(); p != atomDMap.end(); p++ ) {

      TiXmlElement* propElem = new TiXmlElement("property");
      atomElem->LinkEndChild(propElem);

      propElem->SetAttribute("name", p->first);
      propElem->SetAttribute("value", double2String(p->second));
      propElem->SetAttribute("type", "d");
    }

    for (intPropertyMapIterator p = atomIMap.begin(); p != atomIMap.end(); p++ ) {

      TiXmlElement* propElem = new TiXmlElement("property");
      atomElem->LinkEndChild(propElem);

      propElem->SetAttribute("name", p->first);
      propElem->SetAttribute("value", int2String(p->second));
      propElem->SetAttribute("type", "i");
    }
}

// ============================================================
// Function : writeMetalCenter
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeMetalCenter(TiXmlElement* root, metalCenter* pMetalCenter)
{
    TiXmlElement* mcElem = new TiXmlElement("metalCenter");
    root->LinkEndChild(mcElem);

    mcElem->SetAttribute("identity", pMetalCenter->getName());
}

// ============================================================
// Function : writeBonds
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeBonds(TiXmlElement* root, std::map<int, Bond*> bondMap)
{
    for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
      pBond = b->second;

      TiXmlElement* bondElem = new TiXmlElement("bond");
      root->LinkEndChild(bondElem);

      bondElem->SetAttribute("index", int2String(b->first));
      bondElem->SetAttribute("at1", int2String(pBond->atom1->getIndex()));
      bondElem->SetAttribute("at2", int2String(pBond->atom2->getIndex()));
      bondElem->SetAttribute("type", int2String(pBond->type));
      bondElem->SetAttribute("stereo", int2String(pBond->stereo));
      bondElem->SetAttribute("topology", int2String(pBond->topology));
      bondElem->SetAttribute("kind", int2String(pBond->kind));
      bondElem->SetAttribute("rotatable", int2String(pBond->rotatable));

      //bondParam* pBondParam;
    }
}

// ============================================================
// Function : writeAngles
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeAngles(TiXmlElement* root, std::map<ULONG_KIND, Angle*> angleMap)
{
    for (AngleMapIterator a = angleMap.begin(); a != angleMap.end(); a++) {
      pAngle = a->second;

      TiXmlElement* angleElem = new TiXmlElement("angle");
      root->LinkEndChild(angleElem);

      angleElem->SetAttribute("index", uLongKind2String(a->first));
      angleElem->SetAttribute("at1", int2String(pAngle->atom1->getIndex()));
      angleElem->SetAttribute("at2", int2String(pAngle->atom2->getIndex()));
      angleElem->SetAttribute("at3", int2String(pAngle->atom3->getIndex()));

     //angleParam* pAngleParam;
    }
}

// ============================================================
// Function : writeTorsions
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void mtkppParser::writeTorsions(TiXmlElement* root, std::map<ULONG_KIND, Torsion*> torsionMap)
{
    for (TorsionMapIterator t = torsionMap.begin(); t != torsionMap.end(); t++) {
      pTorsion = t->second;

      TiXmlElement* torElem = new TiXmlElement("torsion");
      root->LinkEndChild(torElem);

      torElem->SetAttribute("index", uLongKind2String(t->first));
      torElem->SetAttribute("at1", int2String(pTorsion->atom1->getIndex()));
      torElem->SetAttribute("at2", int2String(pTorsion->atom2->getIndex()));
      torElem->SetAttribute("at3", int2String(pTorsion->atom3->getIndex()));
      torElem->SetAttribute("at4", int2String(pTorsion->atom4->getIndex()));
    }
    //torsionParam* pTorsionParam;
}

#endif // USE_TINYXML

} // MTKpp namespace
