/*!
   \file paramParser.cpp
   \brief Parses parameter xml files using XERCES-C and Trolltech's Qt
   \author Martin Peters

   $Date: 2010/04/29 19:06:19 $
   $Revision: 1.20 $

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

#include "paramParser.h"
#include "Molecule/parameters.h"
#include "StringManip.h"
#include "Utils/constants.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

#include "config.h"

#ifdef USE_XERCES
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#endif // USE_XERCES

namespace MTKpp
{

// ============================================================
// Function : paramParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
paramParser::paramParser(parameters *c):pParameters(c)
{
    setError(0);
    if (!pParameters) {
      setError(1);
      std::string errorMessage = "paramParser::Error, can't find parameters object";
      setErrorMessage(errorMessage);
    }

#ifdef USE_TINYXML
    errorLogger.throwError("paramParser", " Reading XML with TINYXML ", MESSAGE);
#endif // USE_TINYXML

#ifdef USE_XERCES
    errorLogger.throwError("paramParser", " Reading XML with XERCES-C ", MESSAGE);
#endif // USE_XERCES

#ifdef USE_QT
    errorLogger.throwError("paramParser", " Reading XML with QtXml ", MESSAGE);
#endif // USE_QT
}

// ============================================================
// Function : paramParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
paramParser::~paramParser() {}

    // ---------------------------- //
    // -    QT READ FUNCTIONS     - //
    // ---------------------------- //

#ifdef USE_TINYXML
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses parameter xml files using tinyxml
// ------------------------------------------------------------
int paramParser::Read(std::string fileName)
{
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("paramParser", errMessage, INFO);

    TiXmlDocument doc(fileName);
    bool loadOkay = doc.LoadFile();

    if (!loadOkay) {
      //printf ("Could not load parameter file. Error='%s'. Exiting.\n", doc.ErrorDesc() );
      errMessage = " Reading " + fileName;
      errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss << "paramParser"<< errMessage;
      throw MTKException(ss.str());
    }
    else {
      TiXmlNode* node = 0;
      //TiXmlNode* node2 = 0;
      TiXmlElement* element;

      node = doc.RootElement();
      element = node->ToElement();

      if (element->Attribute("name")) {
        this->groupName = element->Attribute("name");
      }

      TiXmlNode* typesNode = node->FirstChild("types");
      if (typesNode) {
        for (element = typesNode->FirstChildElement(); element; element = element->NextSiblingElement()) {
          pAtomType = pParameters->addAtomType();
          pAtomType->groupName = this->groupName;

          if (element->Attribute("mass")) {
            pAtomType->mass = string2Double(element->Attribute("mass"));
          }

          if (element->Attribute("name")) {
            pAtomType->name = element->Attribute("name");
          }

          if (element->Attribute("vdwRadius")) {
            pAtomType->rvalue = string2Double(element->Attribute("vdwRadius"));
          }

          if (element->Attribute("potentialWellDepth")) {
            pAtomType->evalue = string2Double(element->Attribute("potentialWellDepth"));
          }

          if (element->Attribute("atomPolarizability")) {
            pAtomType->atomPolarizability = string2Double(element->Attribute("atomPolarizability"));
          }

          if (element->Attribute("element")) {
            pAtomType->element = element->Attribute("element");
          }

          if (element->Attribute("hybridization")) {
            pAtomType->hybridization = element->Attribute("hybridization");
          }

          if (element->Attribute("description")) {
            pAtomType->description = element->Attribute("description");
          }

          if (element->Attribute("groupName")) {
            pAtomType->groupName = element->Attribute("groupName");
          }

          //std::cout << "|"<< pAtomType->name << "| |" << pAtomType->element << std::endl;
 
          if (pAtomType->element != "") {
            pParameters->setAtomNumber(pAtomType);
          }
        }
      }

      TiXmlNode* bondNode = node->FirstChild("bondLengths");
      if (bondNode) {
        for (element = bondNode->FirstChildElement(); element; element = element->NextSiblingElement()) {
          pBondParam = pParameters->addBondParam();
          pBondParam->groupName = this->groupName;

          if (element->Attribute("t1")) {
            pBondParam->atomType1 = element->Attribute("t1");
          }

          if (element->Attribute("t2")) {
            pBondParam->atomType2 = element->Attribute("t2");
          }

          if (element->Attribute("keq")) {
            pBondParam->keq = string2Double(element->Attribute("keq"));
          }

          if (element->Attribute("req")) {
            pBondParam->req = string2Double(element->Attribute("req"));
          }

          if (element->Attribute("optimize")) {
            std::string optim = element->Attribute("optimize");
            if (optim == "t") {
              pBondParam->optimize = true;
            }
          }
        }
      }

      TiXmlNode* angleNode = node->FirstChild("bondAngles");
      if (angleNode) {
        for (element = angleNode->FirstChildElement(); element; element = element->NextSiblingElement()) {
          pAngleParam = pParameters->addAngleParam();
          pAngleParam->groupName = this->groupName;

          if (element->Attribute("t1")) {
            pAngleParam->atomType1 = element->Attribute("t1");
          }

          if (element->Attribute("t2")) {
            pAngleParam->atomType2 = element->Attribute("t2");
          }

          if (element->Attribute("t3")) {
            pAngleParam->atomType3 = element->Attribute("t3");
          }

          if (element->Attribute("keq")) {
            pAngleParam->keq = string2Double(element->Attribute("keq"));
          }

          if (element->Attribute("req")) {
            pAngleParam->req = string2Double(element->Attribute("req")) * DEG2RAD;
          }

          if (element->Attribute("optimize")) {
            std::string optim = element->Attribute("optimize");
            if (optim == "t") {
              pAngleParam->optimize = true;
            }
          }
        }
      }

      TiXmlNode* torsionNode = node->FirstChild("bondTorsions");
      if (torsionNode) {
        for (element = torsionNode->FirstChildElement(); element; element = element->NextSiblingElement()) {
          pTorsionParam = pParameters->addTorsionParam();
          pTorsionParam->groupName = this->groupName;

          if (element->Attribute("t1")) {
            pTorsionParam->atomType1 = element->Attribute("t1");
          }

          if (element->Attribute("t2")) {
            pTorsionParam->atomType2 = element->Attribute("t2");
          }

          if (element->Attribute("t3")) {
            pTorsionParam->atomType3 = element->Attribute("t3");
          }

          if (element->Attribute("t4")) {
            pTorsionParam->atomType4 = element->Attribute("t4");
          }

          if (element->Attribute("Nt")) {
            pTorsionParam->Nt = string2Double(element->Attribute("Nt"));
          }

          if (element->Attribute("Vn")) {
            pTorsionParam->Vn = string2Double(element->Attribute("Vn"));
          }

          if (element->Attribute("gamma")) {
            pTorsionParam->gamma = string2Double(element->Attribute("gamma")) * DEG2RAD;
          }

          if (element->Attribute("npth")) {
            pTorsionParam->npth = string2Int(element->Attribute("npth"));
          }
        }
      }

      TiXmlNode* impNode = node->FirstChild("bondImpropers");
      if (impNode) {
        for (element = impNode->FirstChildElement(); element; element = element->NextSiblingElement()) {
          pImproperParam = pParameters->addImproperParam();
          pImproperParam->groupName = this->groupName;

          if (element->Attribute("t1")) {
            pImproperParam->atomType1 = element->Attribute("t1");
          }

          if (element->Attribute("t2")) {
            pImproperParam->atomType2 = element->Attribute("t2");
          }

          if (element->Attribute("t3")) {
            pImproperParam->atomType3 = element->Attribute("t3");
          }

          if (element->Attribute("t4")) {
            pImproperParam->atomType4 = element->Attribute("t4");
          }

          if (element->Attribute("Nt")) {
            pImproperParam->Nt = string2Double(element->Attribute("Nt"));
          }

          if (element->Attribute("Vn")) {
            pImproperParam->Vn = string2Double(element->Attribute("Vn"));
          }

          if (element->Attribute("gamma")) {
            pImproperParam->gamma = string2Double(element->Attribute("gamma")) * DEG2RAD;
          }
        }
      }

      // <entry p10="0.0" p12="0.0" t1="HW" t2="OW"/>
      //TiXmlNode* hbondsNode = node->FirstChild("hbonds");

      TiXmlNode* eqAtomsNode = node->FirstChild("equivalentAtoms");
      if (eqAtomsNode) {
        for (node = eqAtomsNode->FirstChild(); node; node = node->NextSibling()) {
          pEquivalentAtomsParam = pParameters->addEquivalentAtomsParam();

          TiXmlElement* nodeElement = node->ToElement();
          if (nodeElement) {
            if (nodeElement->Attribute("org")) {
              pEquivalentAtomsParam->original = nodeElement->Attribute("org");
              //std::cout << "eq atom = " << pEquivalentAtomsParam->original << std::endl;
            }
          }

          for (element = node->FirstChildElement(); element; element = element->NextSiblingElement()) {
            pImproperParam = pParameters->addImproperParam();
            pImproperParam->groupName = this->groupName;

            if (element->Attribute("equil")) {
              pEquivalentAtomsParam->itsEquivalentList.push_back(element->Attribute("equil"));
              //std::cout <<  "     " << element->Attribute("equil") << std::endl;
            }
          }
        }
      }
    }
    this->updateEquivalentAtoms();

    return 0;
/*
      for (node = doc.FirstChild(); node; node = node->NextSibling()) {
        for (node2 = node->FirstChild(); node2; node2 = node2->NextSibling()) {
        std::cout << (*node) << std::endl;
          for (element = node2->FirstChildElement(); element; element = element->NextSiblingElement()) {
            //std::cout << (*element) << std::endl;
            if (element->Attribute("mass")) {
              std::string mass = element->Attribute("mass");
              std::cout << mass << std::endl;
    } } } } }
*/
}
#endif

#ifdef USE_QT
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses parameter xml files using Qt
// ------------------------------------------------------------
int paramParser::Read(std::string fileName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("paramParser", errMessage, INFO);

    QString domErrorMessage = "";
    int domErrorLine = 0;
    int domErrorColumn = 0;

    // Read elements file using DOM into memory
    QDomDocument doc("mydocument");
    QFile file(qFileName);
    if (!file.open(QIODevice::ReadOnly)) {
      errMessage = " Reading " + fileName;
      errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
      return 1;
    }
    if (!doc.setContent(&file, &domErrorMessage, &domErrorLine, &domErrorColumn)) {
      errMessage = " Reading " + fileName + " " + domErrorMessage.toStdString()
                 + " " + int2String(domErrorLine) + " " + int2String(domErrorColumn);
      errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
      file.close();
      return 1;
    }
    file.close();

    // Parse the contents of the elements file
    QString att;
    QDomElement docElem = doc.documentElement();
    if (docElem.tagName() == "parms") {

      if (docElem.hasAttribute("name")) {
        att = docElem.attribute("name");
        this->groupName = att.toStdString();
      }

      QDomNode node = docElem.firstChild();
      while (!node.isNull()) {
        if (node.nodeName() == "types") {
          typeFiller(node);
        }
        if (node.nodeName() == "bondLengths") {
          bondFiller(node);
        }
        if (node.nodeName() == "bondAngles") {
          angleFiller(node);
        }
        if (node.nodeName() == "bondTorsions") {
          torsionFiller(node);
        }
        if (node.nodeName() == "bondImpropers") {
          improperFiller(node);
        }
        if (node.nodeName() == "hbonds") {
          hBondFiller(node);
        }
        if (node.nodeName() == "equivalentAtoms") {
          equivalAtomFiller(node);
        }
        node = node.nextSibling();
      }
    }
    this->updateEquivalentAtoms();

    return 0;
}

// ============================================================
// Function : typeFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::typeFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {
        QDomElement entryElement = entryNode.toElement();
        pAtomType = pParameters->addAtomType();
        pAtomType->groupName = this->groupName;

        if (!entryElement.isNull()) {
          QString att;
          bool bConvertOk = true;

          if (entryElement.hasAttribute("mass")) {
            att = entryElement.attribute("mass");
            double mass = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pAtomType->mass = mass;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("name")) {
            att = entryElement.attribute("name");
            pAtomType->name = att.toStdString();
          }

          if (entryElement.hasAttribute("vdwRadius")) {
            att = entryElement.attribute("vdwRadius");
            double rvalue = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pAtomType->rvalue = rvalue;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("potentialWellDepth")) {
            att = entryElement.attribute("potentialWellDepth");
            double evalue = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pAtomType->evalue = evalue;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("element")) {
            att = entryElement.attribute("element");
            pAtomType->element = att.toStdString();
          }

          if (entryElement.hasAttribute("hybridization")) {
            att = entryElement.attribute("hybridization");
            pAtomType->hybridization = att.toStdString();
          }

          if (entryElement.hasAttribute("description")) {
            att = entryElement.attribute("description");
            pAtomType->description = att.toStdString();
          }

          if (entryElement.hasAttribute("groupName")) {
            att = entryElement.attribute("groupName");
            pAtomType->groupName = att.toStdString();
          }

          if (pAtomType->element != "") {
            pParameters->setAtomNumber(pAtomType);
          }
        }
        entryNode = entryNode.nextSibling();
      }
    }
}

// ============================================================
// Function : bondFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::bondFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {

        QDomElement entryElement = entryNode.toElement();
        pBondParam = pParameters->addBondParam();
        pBondParam->groupName = this->groupName;

        if (!entryElement.isNull()) {
          QString att;
          bool bConvertOk = true;

          if (entryElement.hasAttribute("t1")) {
            att = entryElement.attribute("t1");
            pBondParam->atomType1 = att.toStdString();
          }

          if (entryElement.hasAttribute("t2")) {
            att = entryElement.attribute("t2");
            pBondParam->atomType2 = att.toStdString();
          }

          if (entryElement.hasAttribute("keq")) {
            att = entryElement.attribute("keq");
            double keq = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pBondParam->keq = keq;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("req")) {
            att = entryElement.attribute("req");
            double req = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pBondParam->req = req;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("optimize")) {
            att = entryElement.attribute("optimize");
            if (att == "t") {
              pBondParam->optimize = true;
            }
          }

        }
        entryNode = entryNode.nextSibling();
      }
    }
}

// ============================================================
// Function : angleFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::angleFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {
        QDomElement entryElement = entryNode.toElement();
        pAngleParam = pParameters->addAngleParam();
        pAngleParam->groupName = this->groupName;

        if (!entryElement.isNull()) {
          QString att;
          bool bConvertOk = true;

          if (entryElement.hasAttribute("t1")) {
            att = entryElement.attribute("t1");
            pAngleParam->atomType1 = att.toStdString();
          }

          if (entryElement.hasAttribute("t2")) {
            att = entryElement.attribute("t2");
            pAngleParam->atomType2 = att.toStdString();
          }

          if (entryElement.hasAttribute("t3")) {
            att = entryElement.attribute("t3");
            pAngleParam->atomType3 = att.toStdString();
          }

          if (entryElement.hasAttribute("keq")) {
            att = entryElement.attribute("keq");
            double keq = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pAngleParam->keq = keq;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("req")) {
            att = entryElement.attribute("req");
            double req = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pAngleParam->req = req * DEG2RAD;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("optimize")) {
            att = entryElement.attribute("optimize");
            if (att == "t") {
              pAngleParam->optimize = true;
            }
          }

        }
        entryNode = entryNode.nextSibling();
      }
    }
}

// ============================================================
// Function : torsionFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::torsionFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {
        QDomElement entryElement = entryNode.toElement();
        pTorsionParam = pParameters->addTorsionParam();
        pTorsionParam->groupName = this->groupName;

        if (!entryElement.isNull()) {
          QString att;
          bool bConvertOk = true;

          if (entryElement.hasAttribute("t1")) {
            att = entryElement.attribute("t1");
            pTorsionParam->atomType1 = att.toStdString();
          }

          if (entryElement.hasAttribute("t2")) {
            att = entryElement.attribute("t2");
            pTorsionParam->atomType2 = att.toStdString();
          }

          if (entryElement.hasAttribute("t3")) {
            att = entryElement.attribute("t3");
            pTorsionParam->atomType3 = att.toStdString();
          }

          if (entryElement.hasAttribute("t4")) {
            att = entryElement.attribute("t4");
            pTorsionParam->atomType4 = att.toStdString();
          }

          if (entryElement.hasAttribute("Nt")) {
            att = entryElement.attribute("Nt");
            double Nt = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pTorsionParam->Nt = Nt;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("Vn")) {
            att = entryElement.attribute("Vn");
            double Vn = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pTorsionParam->Vn = Vn;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("gamma")) {
            att = entryElement.attribute("gamma");
            double gamma = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pTorsionParam->gamma = gamma * DEG2RAD;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("npth")) {
            att = entryElement.attribute("npth");
            int npth = att.toInt(&bConvertOk, 10);
            if (bConvertOk) {
              pTorsionParam->npth = npth;
            }
            else {
              errMessage = " Error converting string to integer ";
              errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }
        }
        entryNode = entryNode.nextSibling();
      }
    }
}

// ============================================================
// Function : improperFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::improperFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {
        QDomElement entryElement = entryNode.toElement();
        pImproperParam = pParameters->addImproperParam();
        pImproperParam->groupName = this->groupName;

        if (!entryElement.isNull()) {
          QString att;
          bool bConvertOk = true;

          if (entryElement.hasAttribute("t1")) {
            att = entryElement.attribute("t1");
            pImproperParam->atomType1 = att.toStdString();
          }

          if (entryElement.hasAttribute("t2")) {
            att = entryElement.attribute("t2");
            pImproperParam->atomType2 = att.toStdString();
          }

          if (entryElement.hasAttribute("t3")) {
            att = entryElement.attribute("t3");
            pImproperParam->atomType3 = att.toStdString();
          }

          if (entryElement.hasAttribute("t4")) {
            att = entryElement.attribute("t4");
            pImproperParam->atomType4 = att.toStdString();
          }

          if (entryElement.hasAttribute("Nt")) {
            att = entryElement.attribute("Nt");
            double Nt = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pImproperParam->Nt = Nt;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("Vn")) {
            att = entryElement.attribute("Vn");
            double Vn = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pImproperParam->Vn = Vn;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("gamma")) {
            att = entryElement.attribute("gamma");
            double gamma = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pImproperParam->gamma = gamma * DEG2RAD;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }
        }
        entryNode = entryNode.nextSibling();
      }
    }
}

// ============================================================
// Function : hBondFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::hBondFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {
        QDomElement entryElement = entryNode.toElement();
        pHBondParam = pParameters->addHBondParam();

        if (!entryElement.isNull()) {
          QString att;
          bool bConvertOk = true;

          if (entryElement.hasAttribute("t1")) {
            att = entryElement.attribute("t1");
            pHBondParam->atomType1 = att.toStdString();
          }

          if (entryElement.hasAttribute("t2")) {
            att = entryElement.attribute("t2");
            pHBondParam->atomType2 = att.toStdString();
          }

          if (entryElement.hasAttribute("p10")) {
            att = entryElement.attribute("p10");
            double p10 = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pHBondParam->p10 = p10;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }

          if (entryElement.hasAttribute("p12")) {
            att = entryElement.attribute("p12");
            double p12 = att.toDouble(&bConvertOk);
            if (bConvertOk) {
              pHBondParam->p12 = p12;
            }
            else {
              errMessage = " Error converting string to double ";
              errorLogger.throwError("paramParser", errMessage, MTK_ERROR);
              //exit(1);
              std::stringstream ss;
              ss << "paramParser"<< errMessage;
              throw MTKException(ss.str());
            }
          }
        }
        entryNode = entryNode.nextSibling();
      }
    }
}

// ============================================================
// Function : equivalAtomFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::equivalAtomFiller(QDomNode a)
{
    std::string errMessage;
    QDomElement typeElement = a.toElement();
    QDomNode entryNode = typeElement.firstChild();

    while (!entryNode.isNull()) {
      if (entryNode.nodeName() == "entry" and entryNode.hasAttributes()) {
        QDomElement entryElement = entryNode.toElement();
        pEquivalentAtomsParam = pParameters->addEquivalentAtomsParam();

        if (!entryElement.isNull()) {
          QString att;

          if (entryElement.hasAttribute("org")) {
            att = entryElement.attribute("org");
            pEquivalentAtomsParam->original = att.toStdString();
          }
          QDomNode equilNode = entryElement.firstChild();
          while (!equilNode.isNull()) {
            if (equilNode.nodeName() == "equil" and equilNode.hasAttributes()) {
              QDomElement equilElement = equilNode.toElement();
              if (equilElement.hasAttribute("equil")) {
                att = equilElement.attribute("equil");
                pEquivalentAtomsParam->itsEquivalentList.push_back(att.toStdString());
              }
            }
            equilNode = equilNode.nextSibling();
          }
        }
        entryNode = entryNode.nextSibling();
      }
    }
}
#endif // USE_QT

    // ---------------------------------- //
    // -    XERCES-C READ FUNCTIONS     - //
    // ---------------------------------- //

#ifdef USE_XERCES
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses parameter xml files using XERCES-C
// ------------------------------------------------------------
int paramParser::Read(std::string fileName)
{
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("paramParser", errMessage, INFO);

    setError(0);
    if (!fileExists(fileName)) {
      setError(1);
      std::string errorMessage = "  paramParser::Read Error, Can't Find " + fileName;
      setErrorMessage(errorMessage);
      return 1;
    }

    try {
      XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
      char* message = XMLString::transcode(toCatch.getMessage());
      std::cout << "Error during initialization! :\n"
                << message << "\n";
      XMLString::release(&message);
      return 1;
    }

    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);

#if (XERCES_VERSION_MAJOR == 3)
    DOMLSParser* parser=((DOMImplementationLS*)impl)->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS,0, XMLPlatformUtils::fgMemoryManager, NULL);
#else
    DOMBuilder* parser=((DOMImplementationLS*)impl)->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS,0);

    // optionally you can set some features on this builder
    if (parser->canSetFeature(XMLUni::fgDOMValidation, true))
      parser->setFeature(XMLUni::fgDOMValidation, true);
    if (parser->canSetFeature(XMLUni::fgDOMNamespaces, true))
      parser->setFeature(XMLUni::fgDOMNamespaces, true);
    if (parser->canSetFeature(XMLUni::fgDOMDatatypeNormalization, true))
      parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
#endif

    char* xmlFile = (char*)(fileName.c_str());
    XERCES_CPP_NAMESPACE::DOMDocument *doc = 0;

    // ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
    // parser->setErrorHandler(errHandler);
    // MyDOMErrorHandler* errHandler = new MyDOMErrorHandler();
    // parser->setErrorHandler(errHandler);

    try {
      doc = parser->parseURI(xmlFile);
      this->paramFiller((DOMNode*)doc->getDocumentElement());
      this->updateEquivalentAtoms();
    }
    catch (const XMLException& toCatch) {
      char* message = XMLString::transcode(toCatch.getMessage());
      std::cout << "Exception message is: \n"
                << message << "\n";
      XMLString::release(&message);
      return -1;
    }
    catch (const DOMException& toCatch) {
      char* message = XMLString::transcode(toCatch.msg);
      std::cout << "Exception message is: \n"
                << message << "\n";
      XMLString::release(&message);
      return -1;
    }
    catch (...) {
      std::cout << "Unexpected Exception \n" ;
      return -1;
    }
    //delete xmlFile;

    parser->release();
    //delete errHandler;
    XMLPlatformUtils::Terminate();
    return 0;
}

// ============================================================
// Function : paramFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::paramFiller(DOMNode *rootnode) {
    char *name = XMLString::transcode(rootnode->getNodeName());
    DOMNode *d;
    this->groupName = "";
    if (std::string(name) == std::string("parms")) {
      if (rootnode->hasAttributes()) {
        // get all the attributes of the node
        DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
        int nSize = pAttributes->getLength();
        for (int i = 0; i < nSize; ++i) {
          DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
          // get attribute name
          char *_name = XMLString::transcode(pAttributeNode->getName());
          if (std::string(_name) == std::string("name")) {
            char *paramName = XMLString::transcode(pAttributeNode->getNodeValue());
            this->groupName = std::string(paramName);
            delete paramName;
          }
          delete _name;
        }
      }
      if (this->groupName == "") {
        return;
      }
      for (d = rootnode->getFirstChild(); d != 0; d = d->getNextSibling()) {
         if (d->getNodeType() == 1) {
           if ((std::string)(XC(d->getNodeName())) == "types") {
             typeFiller(d);
           }
           if ((std::string)(XC(d->getNodeName())) == "bondLengths") {
             bondFiller(d);
           }
           if ((std::string)(XC(d->getNodeName())) == "bondAngles") {
             angleFiller(d);
           }
           if ((std::string)(XC(d->getNodeName())) == "bondTorsions") {
             torsionFiller(d);
           }
           if ((std::string)(XC(d->getNodeName())) == "bondImpropers") {
             improperFiller(d);
           }
           if ((std::string)(XC(d->getNodeName())) == "hbonds") {
             hBondFiller(d);
           }
           if ((std::string)(XC(d->getNodeName())) == "equivalentAtoms") {
             equivalAtomFiller(d);
           }
         }
       }
    }
    delete name;
}

// ============================================================
// Function : typeFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::typeFiller(DOMNode *rootnode) {
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pAtomType = pParameters->addAtomType();

          pAtomType->groupName = this->groupName;
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("mass")) {
                double dMass = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pAtomType->mass = dMass;
              }
              if (std::string(name) == std::string("name")) {
                char *cAtomType = XMLString::transcode(pAttributeNode->getNodeValue());
                pAtomType->name = std::string(cAtomType);
                delete cAtomType;
              }
              if (std::string(name) == std::string("vdwRadius")) {
                double dvdw = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pAtomType->rvalue = dvdw;
              }
              if (std::string(name) == std::string("potentialWellDepth")) {
                double dp = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pAtomType->evalue = dp;
              }
              if (std::string(name) == std::string("element")) {
                char *cElement = XMLString::transcode(pAttributeNode->getNodeValue());
                pAtomType->element = std::string(cElement);
                delete cElement;
              }
              if (std::string(name) == std::string("hybridization")) {
                char *cH = XMLString::transcode(pAttributeNode->getNodeValue());
                pAtomType->hybridization = std::string(cH);
                delete cH;
              }
              if (std::string(name) == std::string("description")) {
                char *cD = XMLString::transcode(pAttributeNode->getNodeValue());
                pAtomType->description = std::string(cD);
                delete cD;
              }
              if (std::string(name) == std::string("groupName")) {
                char *cGN = XMLString::transcode(pAttributeNode->getNodeValue());
                pAtomType->groupName = std::string(cGN);
                delete cGN;
              }
              delete name;
            }
          }
          if (pAtomType->element != "") {
            pParameters->setAtomNumber(pAtomType);
          }
        }
      }
    }
}

// ============================================================
// Function : bondFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::bondFiller(DOMNode *rootnode) {
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pBondParam = pParameters->addBondParam();

          pBondParam->groupName = this->groupName;
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("t1")) {
                char *cAtomType1 = XMLString::transcode(pAttributeNode->getNodeValue());
                pBondParam->atomType1 = std::string(cAtomType1);
                delete cAtomType1;
              }
              if (std::string(name) == std::string("t2")) {
                char *cAtomType2 = XMLString::transcode(pAttributeNode->getNodeValue());
                pBondParam->atomType2 = std::string(cAtomType2);
                delete cAtomType2;
              }
              if (std::string(name) == std::string("keq")) {
                double dkeq = strtod(XC(pAttributeNode->getNodeValue()),0);
                pBondParam->keq = dkeq;
              }
              if (std::string(name) == std::string("req")) {
                double dreq = strtod(XC(pAttributeNode->getNodeValue()),0);
                pBondParam->req = dreq;
              }
              if (std::string(name) == std::string("optimize")) {
                char *coptimize = XMLString::transcode(pAttributeNode->getNodeValue());
                if (std::string(coptimize) == "t") {
                  pBondParam->optimize = true;
                }
                delete coptimize;
              }
              delete name;
            }
          }
        }
      }
    }
}

// ============================================================
// Function : angleFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::angleFiller(DOMNode *rootnode) {
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pAngleParam = pParameters->addAngleParam();

          pAngleParam->groupName = this->groupName;
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("t1")) {
                char *cAtomType1 = XMLString::transcode(pAttributeNode->getNodeValue());
                pAngleParam->atomType1 = std::string(cAtomType1);
                delete cAtomType1;
              }
              if (std::string(name) == std::string("t2")) {
                char *cAtomType2 = XMLString::transcode(pAttributeNode->getNodeValue());
                pAngleParam->atomType2 = std::string(cAtomType2);
                delete cAtomType2;
              }
              if (std::string(name) == std::string("t3")) {
                char *cAtomType3 = XMLString::transcode(pAttributeNode->getNodeValue());
                pAngleParam->atomType3 = std::string(cAtomType3);
                delete cAtomType3;
              }
              if (std::string(name) == std::string("keq")) {
                double dkeq = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pAngleParam->keq = dkeq;
              }
              if (std::string(name) == std::string("req")) {
                double dreq = strtod(XC(pAttributeNode->getNodeValue()), 0);
                //pAngleParam->req = dreq / (180.0 / PI);
                pAngleParam->req = dreq * DEG2RAD;
              }
              if (std::string(name) == std::string("optimize")) {
                char *coptimize = XMLString::transcode(pAttributeNode->getNodeValue());
                if (std::string(coptimize) == "t") {
                  pAngleParam->optimize = true;
                }
                delete coptimize;
              }
              delete name;
            }
          }
          //std::cout << "|" << pAngleParam->atomType1
          //          << "|" << pAngleParam->atomType2
          //          << "|" << pAngleParam->atomType3
          //          << "|" << std::endl;
        }
      }
    }
}

// ============================================================
// Function : torsionFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::torsionFiller(DOMNode *rootnode){
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1){
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pTorsionParam = pParameters->addTorsionParam();
          pTorsionParam->groupName = this->groupName;
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("t1")) {
                char *cAtomType1 = XMLString::transcode(pAttributeNode->getNodeValue());
                pTorsionParam->atomType1 = std::string(cAtomType1);
                delete cAtomType1;
              }
              if (std::string(name) == std::string("t2")) {
                char *cAtomType2 = XMLString::transcode(pAttributeNode->getNodeValue());
                pTorsionParam->atomType2 = std::string(cAtomType2);
                delete cAtomType2;
              }
              if (std::string(name) == std::string("t3")) {
                char *cAtomType3 = XMLString::transcode(pAttributeNode->getNodeValue());
                pTorsionParam->atomType3 = std::string(cAtomType3);
                delete cAtomType3;
              }
              if (std::string(name) == std::string("t4")) {
                char *cAtomType4 = XMLString::transcode(pAttributeNode->getNodeValue());
                pTorsionParam->atomType4 = std::string(cAtomType4);
                delete cAtomType4;
              }
              if (std::string(name) == std::string("Nt")) {
                double dnt = strtod(XC(pAttributeNode->getNodeValue()),0);
                pTorsionParam->Nt = dnt;
              }
              if (std::string(name) == std::string("Vn")) {
                double dvn = strtod(XC(pAttributeNode->getNodeValue()),0);
                pTorsionParam->Vn = dvn;
              }
              if (std::string(name) == std::string("gamma")) {
                double dgamma = strtod(XC(pAttributeNode->getNodeValue()),0);
                pTorsionParam->gamma = dgamma * DEG2RAD; //  deg2rad
              }
              if (std::string(name) == std::string("npth")) {
                char *cNpth = XMLString::transcode(pAttributeNode->getNodeValue());
                pTorsionParam->npth = atoi(cNpth);
                delete cNpth;
              }
              delete name;
            }
          }
        }
      }
    }
}

// ============================================================
// Function : improperFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::improperFiller(DOMNode *rootnode) {
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pImproperParam = pParameters->addImproperParam();
          pImproperParam->groupName = this->groupName;
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("t1")) {
                char *cAtomType1 = XMLString::transcode(pAttributeNode->getNodeValue());
                pImproperParam->atomType1 = std::string(cAtomType1);
                delete cAtomType1;
              }
              if (std::string(name) == std::string("t2")) {
                char *cAtomType2 = XMLString::transcode(pAttributeNode->getNodeValue());
                pImproperParam->atomType2 = std::string(cAtomType2);
                delete cAtomType2;
              }
              if (std::string(name) == std::string("t3")) {
                char *cAtomType3 = XMLString::transcode(pAttributeNode->getNodeValue());
                pImproperParam->atomType3 = std::string(cAtomType3);
                delete cAtomType3;
              }
              if (std::string(name) == std::string("t4")) {
                char *cAtomType4 = XMLString::transcode(pAttributeNode->getNodeValue());
                pImproperParam->atomType4 = std::string(cAtomType4);
                delete cAtomType4;
              }
              if (std::string(name) == std::string("Nt")) {
                double dnt = strtod(XC(pAttributeNode->getNodeValue()),0);
                pImproperParam->Nt = dnt;
              }
              if (std::string(name) == std::string("Vn")) {
                double dvn = strtod(XC(pAttributeNode->getNodeValue()),0);
                pImproperParam->Vn = dvn;
              }
              if (std::string(name) == std::string("gamma")) {
                double dgamma = strtod(XC(pAttributeNode->getNodeValue()),0);
                pImproperParam->gamma = dgamma * DEG2RAD; // deg2rad
              }
              delete name;
            }
          }
        }
      }
    }
}

// ============================================================
// Function : hBondFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::hBondFiller(DOMNode *rootnode) {
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pHBondParam = pParameters->addHBondParam();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("t1")) {
                char *cAtomType1 = XMLString::transcode(pAttributeNode->getNodeValue());
                pHBondParam->atomType1 = std::string(cAtomType1);
                delete cAtomType1;
              }
              if (std::string(name) == std::string("t2")) {
                char *cAtomType2 = XMLString::transcode(pAttributeNode->getNodeValue());
                pHBondParam->atomType2 = std::string(cAtomType2);
                delete cAtomType2;
              }
              if (std::string(name) == std::string("p10")) {
                double dp10 = strtod(XC(pAttributeNode->getNodeValue()),0);
                pHBondParam->p10 = dp10;
              }
              if (std::string(name) == std::string("p12")) {
                double dp12 = strtod(XC(pAttributeNode->getNodeValue()),0);
                pHBondParam->p12 = dp12;
              }
              delete name;
            }
          }
        }
      }
    }
}

// ============================================================
// Function : equivalAtomFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void paramParser::equivalAtomFiller(DOMNode *rootnode) {
    DOMNode *c;
    DOMNode *e;
    for (c=rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "entry") {
          pEquivalentAtomsParam = pParameters->addEquivalentAtomsParam();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("org")){
                char *cOrg = XMLString::transcode(pAttributeNode->getNodeValue());
                pEquivalentAtomsParam->original = std::string(cOrg);
                delete cOrg;
              }
              delete name;
            }
          }
          for (e = c->getFirstChild(); e != 0; e = e->getNextSibling()) {
            if (e->getNodeType() == 1) {
              if ((std::string)(XC(e->getNodeName())) == "equil") {
                if (e->hasAttributes()){
                  DOMNamedNodeMap *pAttributes2 = e->getAttributes();
                  int nSize2 = pAttributes2->getLength();
                  for (int i = 0; i < nSize2; ++i) {
                    DOMAttr *pAttributeNode2 = (DOMAttr*) pAttributes2->item(i);
                    // get attribute name
                    char *name2 = XMLString::transcode(pAttributeNode2->getName());
                    if (std::string(name2) == std::string("equil")) {
                      char *cEquil = XMLString::transcode(pAttributeNode2->getNodeValue());
                      pEquivalentAtomsParam->itsEquivalentList.push_back(std::string(cEquil));
                      delete cEquil;
                    }
                    delete name2;
                  }
                }
              }
            }
          }
        }
      }
    }
}
#endif // USE_XERCES

// ============================================================
// Function : updateEquivalentAtoms
// ------------------------------------------------------------
// Assign equivalent atom types vdw radii and potential well depths
// ------------------------------------------------------------
void paramParser::updateEquivalentAtoms()
{
    std::vector<equivalentAtomsParam*> equivalentAtomsList;
    equivalentAtomsList = pParameters->getEquivalentAtomList();

    atomType* pAt1;
    atomType* pAt2;
    std::string a;

    for (unsigned int i = 0; i < equivalentAtomsList.size(); i++) {
      pEquivalentAtomsParam = equivalentAtomsList[i];
      pAt1 = pParameters->getAtomType(pEquivalentAtomsParam->original);

      for (unsigned int j = 0; j < pEquivalentAtomsParam->itsEquivalentList.size(); j++) {
        a = pEquivalentAtomsParam->itsEquivalentList[j];
        pAt2 = pParameters->getAtomType(a);
        if (pAt2) {
          pAt2->mass = pAt1->mass;
          pAt2->rvalue = pAt1->rvalue;
          pAt2->evalue = pAt1->evalue;
        }
      }
    }
}

#ifdef USE_QT
    // ----------------------------- //
    // -    Qt WRITE FUNCTIONS     - //
    // ----------------------------- //

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write library xml files using Qt
// ------------------------------------------------------------
int paramParser::Write(std::string fileName, std::string groupName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("paramParser", errMessage, INFO);

    // Write parameter file using DOM into memory
    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QDomElement rootElem = doc.createElement("parms");
    rootElem.setAttribute("name", QString::fromStdString(groupName));
    doc.appendChild(rootElem);

    // ATOM TYPES
    QDomElement typesElem = doc.createElement("types");
    rootElem.appendChild(typesElem);

    typedef std::vector<atomType*>::iterator atomTypeIterator;
    std::vector<atomType*> itsTypeList = pParameters->getAtomTypes();
    for (atomTypeIterator a = itsTypeList.begin(); a != itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->groupName == groupName) {
        QDomElement typeElem = doc.createElement("entry");
        typesElem.appendChild(typeElem);
        typeElem.setAttribute("name"          , QString::fromStdString(pAtomType->name));
        typeElem.setAttribute("mass"          , QString::fromStdString(double2String(pAtomType->mass).c_str()));
        typeElem.setAttribute("vdwRadius"     , QString::fromStdString(double2String(pAtomType->rvalue).c_str()));
        typeElem.setAttribute("potentialWellDepth" , QString::fromStdString(double2String(pAtomType->evalue).c_str()));
        typeElem.setAttribute("element"       , QString::fromStdString(pAtomType->element));
        typeElem.setAttribute("hybridization" , QString::fromStdString(pAtomType->hybridization));
        typeElem.setAttribute("groupName"     , QString::fromStdString(pAtomType->groupName));
        typeElem.setAttribute("description"   , QString::fromStdString(pAtomType->description));
      }
    }

    // BONDS
    QDomElement bondsElem = doc.createElement("bondLengths");
    rootElem.appendChild(bondsElem);

    typedef std::vector<bondParam*>::iterator bondParamIterator;
    std::vector<bondParam*> bBList = pParameters->getBondParams();
    for (bondParamIterator b = bBList.begin(); b != bBList.end(); b++) {
      pBondParam = *b;
      if (pBondParam->groupName == groupName) {
        QDomElement bondElem = doc.createElement("entry");
        bondsElem.appendChild(bondElem);
        bondElem.setAttribute("t1", QString::fromStdString(pBondParam->atomType1));
        bondElem.setAttribute("t2", QString::fromStdString(pBondParam->atomType2));

        //bondElem.setAttribute("keq", QString::fromStdString(double2String(pBondParam->keq).c_str()));
        //bondElem.setAttribute("req", QString::fromStdString(double2String(pBondParam->req).c_str()));

        bondElem.setAttribute("keq", QString::number(pBondParam->keq, 'f', 2));
        bondElem.setAttribute("req", QString::number(pBondParam->req, 'f', 3)); // doesn't work for OW-HW, HW-HW

        bondElem.setAttribute("groupName", QString::fromStdString(pBondParam->groupName));

        if (pBondParam->optimize) {
          bondElem.setAttribute("optimize", "t");
        }
      }
    }

    // ANGLES
    QDomElement anglesElem = doc.createElement("bondAngles");
    rootElem.appendChild(anglesElem);

    typedef std::vector<angleParam*>::iterator angleParamIterator;
    std::vector<angleParam*> bAList = pParameters->getAngleParams();
    for (angleParamIterator a = bAList.begin(); a != bAList.end(); a++) {
      pAngleParam = *a;
      if (pAngleParam->groupName == groupName) {
        QDomElement angleElem = doc.createElement("entry");
        anglesElem.appendChild(angleElem);
        angleElem.setAttribute("t1", QString::fromStdString(pAngleParam->atomType1));
        angleElem.setAttribute("t2", QString::fromStdString(pAngleParam->atomType2));
        angleElem.setAttribute("t3", QString::fromStdString(pAngleParam->atomType3));

        //angleElem.setAttribute("keq", QString::fromStdString(double2String(pAngleParam->keq).c_str()));
        //angleElem.setAttribute("req", QString::fromStdString(double2String(pAngleParam->req * RAD2DEG).c_str()));

        angleElem.setAttribute("keq", QString::number(pAngleParam->keq, 'f', 2));
        angleElem.setAttribute("req", QString::number(pAngleParam->req * RAD2DEG, 'f', 2));

        angleElem.setAttribute("groupName", QString::fromStdString(pAngleParam->groupName));
        if (pAngleParam->optimize) {
          angleElem.setAttribute("optimize", "t");
        }
      }
    }

    // TORSIONS
    QDomElement torsionsElem = doc.createElement("bondTorsions");
    rootElem.appendChild(torsionsElem);

    typedef std::vector<torsionParam*>::iterator torsionParamIterator;
    std::vector<torsionParam*> bTList = pParameters->getTorsionParams();
    for (torsionParamIterator a = bTList.begin(); a != bTList.end(); a++) {
      pTorsionParam = *a;
      if (pTorsionParam->groupName == groupName) {
        QDomElement torElem = doc.createElement("entry");
        torsionsElem.appendChild(torElem);
        torElem.setAttribute("t1", QString::fromStdString(pTorsionParam->atomType1));
        torElem.setAttribute("t2", QString::fromStdString(pTorsionParam->atomType2));
        torElem.setAttribute("t3", QString::fromStdString(pTorsionParam->atomType3));
        torElem.setAttribute("t4", QString::fromStdString(pTorsionParam->atomType4));
        torElem.setAttribute("Nt", QString::fromStdString(double2String(pTorsionParam->Nt).c_str()));
        torElem.setAttribute("Vn", QString::fromStdString(double2String(pTorsionParam->Vn).c_str()));
        torElem.setAttribute("gamma", QString::fromStdString(double2String(pTorsionParam->gamma * RAD2DEG).c_str()));
        torElem.setAttribute("npth", QString::fromStdString(int2String(pTorsionParam->npth).c_str()));
        torElem.setAttribute("groupName", QString::fromStdString(pTorsionParam->groupName));
        //if (torElem->optimize) {
        //  torElem.setAttribute("optimize", "t");
        //}
      }
    }

    // IMPROPERS
    QDomElement impsElem = doc.createElement("bondImpropers");
    rootElem.appendChild(impsElem);

    typedef std::vector<improperParam*>::iterator improperParamIterator;
    std::vector<improperParam*> bIList = pParameters->getImproperParams();
    for (improperParamIterator a = bIList.begin(); a != bIList.end(); a++) {
      pImproperParam = *a;
      if (pImproperParam->groupName == groupName) {
        QDomElement impElem = doc.createElement("entry");
        impsElem.appendChild(impElem);
        impElem.setAttribute("t1", QString::fromStdString(pImproperParam->atomType1));
        impElem.setAttribute("t2", QString::fromStdString(pImproperParam->atomType2));
        impElem.setAttribute("t3", QString::fromStdString(pImproperParam->atomType3));
        impElem.setAttribute("t4", QString::fromStdString(pImproperParam->atomType4));
        impElem.setAttribute("Nt", QString::fromStdString(double2String(pImproperParam->Nt).c_str()));
        impElem.setAttribute("Vn", QString::fromStdString(double2String(pImproperParam->Vn).c_str()));
        impElem.setAttribute("gamma", QString::fromStdString(double2String(pImproperParam->gamma * RAD2DEG).c_str()));
        impElem.setAttribute("groupName", QString::fromStdString(pImproperParam->groupName));
        //if (impElem->optimize) {
        //  impElem.setAttribute("optimize", "t");
        //}
      }
    }

    // Write dom xml to file
    QString domXml = doc.toString();

    QFile file(qFileName);
    if (!file.open(QIODevice::WriteOnly)) {
      errMessage = " Writing " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      return 1;
    }

    QTextStream out(&file);
    out << domXml;
    file.close();

    return 0;
}
#endif // USE_QT

    // ----------------------------------- //
    // -     TINYXML WRITE FUNCTIONS     - //
    // ----------------------------------- //

#ifdef USE_TINYXML
// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using tinyxml
// ------------------------------------------------------------
int paramParser::Write(std::string fileName, std::string groupName)
{
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("paramParser", errMessage, INFO);

    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

    // Write parameter file using DOM into memory
    //QString mtkppInfo = PACKAGE_TARNAME;
    //QDomDocument doc(mtkppInfo);

    TiXmlElement* root = new TiXmlElement( "parms" );
    doc.LinkEndChild( root );
    root->SetAttribute("name", groupName);

    // ATOM TYPES
    TiXmlElement* types = new TiXmlElement("types");
    root->LinkEndChild(types);

    typedef std::vector<atomType*>::iterator atomTypeIterator;
    std::vector<atomType*> itsTypeList = pParameters->getAtomTypes();
    for (atomTypeIterator a = itsTypeList.begin(); a != itsTypeList.end(); a++) {
      pAtomType = *a;
      if (pAtomType->groupName == groupName) {

        TiXmlElement* entry = new TiXmlElement("entry");
        types->LinkEndChild(entry);

        entry->SetAttribute("name"          , pAtomType->name);
        entry->SetAttribute("mass"          , double2String(pAtomType->mass).c_str());
        entry->SetAttribute("vdwRadius"     , double2String(pAtomType->rvalue).c_str());
        entry->SetAttribute("potentialWellDepth" , double2String(pAtomType->evalue).c_str());
        entry->SetAttribute("atomPolarizability" , double2String(pAtomType->atomPolarizability).c_str());
        entry->SetAttribute("element"       , pAtomType->element);
        entry->SetAttribute("hybridization" , pAtomType->hybridization);
        entry->SetAttribute("groupName"     , pAtomType->groupName);
        entry->SetAttribute("description"   , pAtomType->description);
      }
    }

    // BONDS
    TiXmlElement* bonds = new TiXmlElement("bondLengths");
    root->LinkEndChild(bonds);

    typedef std::vector<bondParam*>::iterator bondParamIterator;
    std::vector<bondParam*> bBList = pParameters->getBondParams();
    for (bondParamIterator b = bBList.begin(); b != bBList.end(); b++) {
      pBondParam = *b;
      if (pBondParam->groupName == groupName) {

        TiXmlElement* entry = new TiXmlElement("entry");
        bonds->LinkEndChild(entry);

        entry->SetAttribute("t1", pBondParam->atomType1);
        entry->SetAttribute("t2", pBondParam->atomType2);

        entry->SetAttribute("keq", double2String(pBondParam->keq));
        entry->SetAttribute("req", double2String(pBondParam->req));

        //bondElem.setAttribute("keq", QString::number(pBondParam->keq, 'f', 2));
        //bondElem.setAttribute("req", QString::number(pBondParam->req, 'f', 3)); // doesn't work for OW-HW, HW-HW

        entry->SetAttribute("groupName", pBondParam->groupName);

        if (pBondParam->optimize) {
          entry->SetAttribute("optimize", "t");
        }
      }
    }

    // ANGLES
    TiXmlElement* angles = new TiXmlElement("bondAngles");
    root->LinkEndChild(angles);

    typedef std::vector<angleParam*>::iterator angleParamIterator;
    std::vector<angleParam*> bAList = pParameters->getAngleParams();
    for (angleParamIterator a = bAList.begin(); a != bAList.end(); a++) {
      pAngleParam = *a;
      if (pAngleParam->groupName == groupName) {

        TiXmlElement* entry = new TiXmlElement("entry");
        angles->LinkEndChild(entry);

        entry->SetAttribute("t1", pAngleParam->atomType1);
        entry->SetAttribute("t2", pAngleParam->atomType2);
        entry->SetAttribute("t3", pAngleParam->atomType3);

        entry->SetAttribute("keq", double2String(pAngleParam->keq));
        entry->SetAttribute("req", double2String(pAngleParam->req * RAD2DEG));

        //entry->SetAttribute("keq", QString::number(pAngleParam->keq, 'f', 2));
        //entry->SetAttribute("req", QString::number(pAngleParam->req * RAD2DEG, 'f', 2));

        entry->SetAttribute("groupName", pAngleParam->groupName);
        if (pAngleParam->optimize) {
          entry->SetAttribute("optimize", "t");
        }
      }
    }

    // TORSIONS
    TiXmlElement* torsions = new TiXmlElement("bondTorsions");
    root->LinkEndChild(torsions);

    typedef std::vector<torsionParam*>::iterator torsionParamIterator;
    std::vector<torsionParam*> bTList = pParameters->getTorsionParams();
    for (torsionParamIterator a = bTList.begin(); a != bTList.end(); a++) {
      pTorsionParam = *a;
      if (pTorsionParam->groupName == groupName) {

        TiXmlElement* entry = new TiXmlElement("entry");
        torsions->LinkEndChild(entry);

        entry->SetAttribute("t1", pTorsionParam->atomType1);
        entry->SetAttribute("t2", pTorsionParam->atomType2);
        entry->SetAttribute("t3", pTorsionParam->atomType3);
        entry->SetAttribute("t4", pTorsionParam->atomType4);
        entry->SetAttribute("Nt", double2String(pTorsionParam->Nt));
        entry->SetAttribute("Vn", double2String(pTorsionParam->Vn));
        entry->SetAttribute("gamma", double2String(pTorsionParam->gamma * RAD2DEG));
        entry->SetAttribute("npth", int2String(pTorsionParam->npth));
        entry->SetAttribute("groupName", pTorsionParam->groupName);
        //if (torElem->optimize) {
        //  entry->SetAttribute("optimize", "t");
        //}
      }
    }

    // IMPROPERS
    TiXmlElement* impropers = new TiXmlElement("bondImpropers");
    root->LinkEndChild(impropers);

    typedef std::vector<improperParam*>::iterator improperParamIterator;
    std::vector<improperParam*> bIList = pParameters->getImproperParams();
    for (improperParamIterator a = bIList.begin(); a != bIList.end(); a++) {
      pImproperParam = *a;
      if (pImproperParam->groupName == groupName) {

        TiXmlElement* entry = new TiXmlElement("entry");
        impropers->LinkEndChild(entry);

        entry->SetAttribute("t1", pImproperParam->atomType1);
        entry->SetAttribute("t2", pImproperParam->atomType2);
        entry->SetAttribute("t3", pImproperParam->atomType3);
        entry->SetAttribute("t4", pImproperParam->atomType4);
        entry->SetAttribute("Nt", double2String(pImproperParam->Nt));
        entry->SetAttribute("Vn", double2String(pImproperParam->Vn));
        entry->SetAttribute("gamma", double2String(pImproperParam->gamma * RAD2DEG));
        entry->SetAttribute("groupName", pImproperParam->groupName);
        //if (impElem->optimize) {
        //  impElem.setAttribute("optimize", "t");
        //}
      }
    }

    doc.SaveFile(fileName);

    //dump_to_stdout(&doc);

    return 0;
}

#endif // USE_TINYXML

    // ----------------------------------- //
    // -    XERCES-C WRITE FUNCTIONS     - //
    // ----------------------------------- //

#ifdef USE_XERCES
// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using xercesc
// ------------------------------------------------------------
int paramParser::Write(std::string fileName, std::string groupName)
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
      return 1;
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
          DOMLSSerializer         *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer(XMLPlatformUtils::fgMemoryManager);
#else
          DOMWriter         *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
#endif

          //DOMErrorHandler *myErrorHandler = new DOMPrintErrorHandler();
          //theSerializer->setErrorHandler(myErrorHandler);

          XERCES_CPP_NAMESPACE::DOMDocument* doc = impl->createDocument(
            0,          // root element namespace URI.
            X("parms"), // root element name
            0);         // document type object (DTD).

          XERCES_CPP_NAMESPACE::DOMElement* rootElem = doc->getDocumentElement();
          rootElem->setAttribute(X("name"), X(groupName.c_str()));

          // ATOM TYPES
          DOMElement* typesElem = doc->createElement(X("types"));
          rootElem->appendChild(typesElem);

          typedef std::vector<atomType*>::iterator atomTypeIterator;
          std::vector<atomType*> itsTypeList = pParameters->getAtomTypes();
          for (atomTypeIterator a = itsTypeList.begin(); a != itsTypeList.end(); a++) {
            pAtomType = *a;
            if (pAtomType->groupName == groupName) {
              DOMElement* typeElem = doc->createElement(X("entry"));
              typesElem->appendChild(typeElem);
              typeElem->setAttribute(X("name"), X(pAtomType->name.c_str()));
              typeElem->setAttribute(X("mass"), X(double2String(pAtomType->mass).c_str()));
              typeElem->setAttribute(X("vdwRadius"), X(double2String(pAtomType->rvalue).c_str()));
              typeElem->setAttribute(X("potentialWellDepth"), X(double2String(pAtomType->evalue).c_str()));
              typeElem->setAttribute(X("atomPolarizability"), X(double2String(pAtomType->atomPolarizability).c_str()));
              typeElem->setAttribute(X("element"), X(pAtomType->element.c_str()));
              typeElem->setAttribute(X("hybridization"), X(pAtomType->hybridization.c_str()));
              typeElem->setAttribute(X("groupName"), X(pAtomType->groupName.c_str()));
              typeElem->setAttribute(X("description"), X(pAtomType->description.c_str()));
            }
          }

          // BONDS
          DOMElement* bondsElem = doc->createElement(X("bondLengths"));
          rootElem->appendChild(bondsElem);

          typedef std::vector<bondParam*>::iterator bondParamIterator;
          std::vector<bondParam*> bBList = pParameters->getBondParams();
          for (bondParamIterator b = bBList.begin(); b != bBList.end(); b++) {
            pBondParam = *b;
            if (pBondParam->groupName == groupName) {
              DOMElement* bondElem = doc->createElement(X("entry"));
              bondsElem->appendChild(bondElem);
              bondElem->setAttribute(X("t1"), X(pBondParam->atomType1.c_str()));
              bondElem->setAttribute(X("t2"), X(pBondParam->atomType2.c_str()));

              //bondElem->setAttribute(X("keq"), X(double2String(pBondParam->keq).c_str()));
              //bondElem->setAttribute(X("req"), X(double2String(pBondParam->req).c_str()));

              bondElem->setAttribute(X("keq"), X(double2String(pBondParam->keq, 2).c_str()));
              bondElem->setAttribute(X("req"), X(double2String(pBondParam->req, 3).c_str())); // // doesn't work for OW-HW, HW-HW

              bondElem->setAttribute(X("groupName"), X(pBondParam->groupName.c_str()));

              if (pBondParam->optimize) {
                bondElem->setAttribute(X("optimize"), X("t"));
              }
            }
          }

          // ANGLES
          DOMElement* anglesElem = doc->createElement(X("bondAngles"));
          rootElem->appendChild(anglesElem);

          typedef std::vector<angleParam*>::iterator angleParamIterator;
          std::vector<angleParam*> bAList = pParameters->getAngleParams();
          for (angleParamIterator a = bAList.begin(); a != bAList.end(); a++) {
            pAngleParam = *a;
            if (pAngleParam->groupName == groupName) {
              DOMElement* angleElem = doc->createElement(X("entry"));
              anglesElem->appendChild(angleElem);
              angleElem->setAttribute(X("t1"), X(pAngleParam->atomType1.c_str()));
              angleElem->setAttribute(X("t2"), X(pAngleParam->atomType2.c_str()));
              angleElem->setAttribute(X("t3"), X(pAngleParam->atomType3.c_str()));

              //angleElem->setAttribute(X("keq"), X(double2String(pAngleParam->keq).c_str()));
              //angleElem->setAttribute(X("req"), X(double2String(pAngleParam->req * RAD2DEG).c_str()));

              angleElem->setAttribute(X("keq"), X(double2String(pAngleParam->keq, 2).c_str()));
              angleElem->setAttribute(X("req"), X(double2String(pAngleParam->req * RAD2DEG, 2).c_str()));

              angleElem->setAttribute(X("groupName"), X(pAngleParam->groupName.c_str()));
              if (pAngleParam->optimize) {
                angleElem->setAttribute(X("optimize"), X("t"));
              }
            }
          }

          // TORSIONS
          DOMElement* torsionsElem = doc->createElement(X("bondTorsions"));
          rootElem->appendChild(torsionsElem);

          typedef std::vector<torsionParam*>::iterator torsionParamIterator;
          std::vector<torsionParam*> bTList = pParameters->getTorsionParams();
          for (torsionParamIterator a = bTList.begin(); a != bTList.end(); a++) {
            pTorsionParam = *a;
            if (pTorsionParam->groupName == groupName) {
              DOMElement* torElem = doc->createElement(X("entry"));
              torsionsElem->appendChild(torElem);
              torElem->setAttribute(X("t1"), X(pTorsionParam->atomType1.c_str()));
              torElem->setAttribute(X("t2"), X(pTorsionParam->atomType2.c_str()));
              torElem->setAttribute(X("t3"), X(pTorsionParam->atomType3.c_str()));
              torElem->setAttribute(X("t4"), X(pTorsionParam->atomType4.c_str()));
              torElem->setAttribute(X("Nt"), X(double2String(pTorsionParam->Nt).c_str()));
              torElem->setAttribute(X("Vn"), X(double2String(pTorsionParam->Vn).c_str()));
              torElem->setAttribute(X("gamma"), X(double2String(pTorsionParam->gamma * RAD2DEG).c_str()));
              torElem->setAttribute(X("npth"), X(int2String(pTorsionParam->npth).c_str()));
              torElem->setAttribute(X("groupName"), X(pTorsionParam->groupName.c_str()));
              //if (torElem->optimize) {
              //  torElem->setAttribute(X("optimize"), X("t"));
              //}
            }
          }

          // IMPROPERS
          DOMElement* impsElem = doc->createElement(X("bondImpropers"));
          rootElem->appendChild(impsElem);

          typedef std::vector<improperParam*>::iterator improperParamIterator;
          std::vector<improperParam*> bIList = pParameters->getImproperParams();
          for (improperParamIterator a = bIList.begin(); a != bIList.end(); a++) {
            pImproperParam = *a;
            if (pImproperParam->groupName == groupName) {
              DOMElement* impElem = doc->createElement(X("entry"));
              impsElem->appendChild(impElem);
              impElem->setAttribute(X("t1"), X(pImproperParam->atomType1.c_str()));
              impElem->setAttribute(X("t2"), X(pImproperParam->atomType2.c_str()));
              impElem->setAttribute(X("t3"), X(pImproperParam->atomType3.c_str()));
              impElem->setAttribute(X("t4"), X(pImproperParam->atomType4.c_str()));
              impElem->setAttribute(X("Nt"), X(double2String(pImproperParam->Nt).c_str()));
              impElem->setAttribute(X("Vn"), X(double2String(pImproperParam->Vn).c_str()));
              impElem->setAttribute(X("gamma"), X(double2String(pImproperParam->gamma * RAD2DEG).c_str()));
              impElem->setAttribute(X("groupName"), X(pImproperParam->groupName.c_str()));
              //if (impElem->optimize) {
              //  impElem->setAttribute(X("optimize"), X("t"));
              //}
            }
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
          std::cout << " paramParser::Write(" << fileName << ", " << groupName << std::endl;
          errorCode = 3;
        }
      }
      else {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorCode = 4;
      }
    }

    XMLPlatformUtils::Terminate();
    return errorCode;
}
#endif // USE_XERCES

} // MTKpp namespace

