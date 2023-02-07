/*!
   \file elementParser.cpp
   \brief Parses elements xml file using XERCES-C and Trolltech's Qt
   \author Martin Peters

   $Date: 2010/03/29 20:39:34 $
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

#include "elementParser.h"
#include "Molecule/element.h"
#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

#include "StringManip.h"

namespace MTKpp
{

// ============================================================
// Function : elementParser()
// ------------------------------------------------------------
// Constructor for the elementParser with tinyxml xml library
// ============================================================
#ifdef USE_TINYXML
elementParser::elementParser(elements *e):pElements(e)
{
    errorLogger.throwError("elementParser", " Reading XML with TINYXML ", MESSAGE);
}
#endif // USE_TINYXML

// ============================================================
// Function : elementParser()
// ------------------------------------------------------------
// Constructor for the elementParser with xerces-c xml library
// ============================================================
#ifdef USE_XERCES
elementParser::elementParser(elements *e):pElements(e)
{
    errorLogger.throwError("elementParser", " Reading XML with XERCES-C ", MESSAGE);
}
#endif // USE_XERCES

// ============================================================
// Function : elementParser()
// ------------------------------------------------------------
// Constructor for the elementParser with Qt xml library
// ============================================================
#ifdef USE_QT
elementParser::elementParser(elements *e):pElements(e)
{
    errorLogger.throwError("elementParser", " Reading XML with QtXml ", MESSAGE);
}
#endif // USE_QT

// ============================================================
// Function :elementParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
elementParser::~elementParser() {}

// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses element xml file using tinyxml
// ------------------------------------------------------------
#ifdef USE_TINYXML
int elementParser::Read(std::string fileName)
{
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("elementParser", errMessage, INFO);

    TiXmlDocument doc(fileName);
    bool loadOkay = doc.LoadFile();

    if ( !loadOkay ) {
      //printf ("Could not load elements file. Error='%s'. Exiting.\n", doc.ErrorDesc() );
      //exit( 1 );
      errMessage += " error ";
      std::stringstream ss;
      ss << "elementParser" << errMessage;
      throw MTKException(ss.str());
    }
    else {
      //std::string errMessage = " Read file ";
      //errorLogger.throwError("elementParser", errMessage, INFO);

      //std::cout << doc << std::endl;

      // Parse the contents of the elements file
/*
      TiXmlElement* root = doc.FirstChildElement("elementlist");
      if ( root ) {
        std::cout << " root " << std::endl;

        TiXmlElement* element = root->FirstChildElement("element");
        if ( element ) {
          std::cout << " element " << std::endl;
        }
      }
*/
      TiXmlNode* node = 0;
      TiXmlElement* element;

      node = doc.RootElement();

      //int count = 0;
      for (node = doc.FirstChild(); node; node = node->NextSibling()) {
        //count++;

        //int count2 = 0;
        for (element = node->FirstChildElement(); element; element = element->NextSiblingElement()) {
          //std::cout << (*element) << std::endl;

          pElement = pElements->addElement();

          if (element->Attribute("symbol")) {
            std::string symbol = element->Attribute("symbol");
            pElements->setElementName(symbol);
          }

          if (element->Attribute("number")) {
            std::string numberStr = element->Attribute("number");
            int number = string2Int(numberStr);
            pElement->number = number;
          }

          if (element->Attribute("name")) {
            std::string name = element->Attribute("name");
            pElement->name = name;
          }

          if (element->Attribute("mass")) {
            std::string massStr = element->Attribute("mass");
            double mass = string2Double(massStr);
            pElement->mass = mass;
          }

          if (element->Attribute("group")) {
            std::string groupStr = element->Attribute("group");
            int group = string2Int(groupStr);
            pElement->group = group;
          }

          if (element->Attribute("period")) {
            std::string periodStr = element->Attribute("period");
            int period = string2Int(periodStr);
            pElement->period = period;
          }

          if (element->Attribute("vdWRadius")) {
            std::string vdWRadiusStr = element->Attribute("vdWRadius");
            double vdWRadius = string2Double(vdWRadiusStr);
            pElement->vdWRadius = vdWRadius;
          }

          if (element->Attribute("red")) {
            std::string redStr = element->Attribute("red");
            double red = string2Double(redStr);
            pElement->red = red;
          }

          if (element->Attribute("green")) {
            std::string greenStr = element->Attribute("green");
            double green = string2Double(greenStr);
            pElement->green = green;
          }

          if (element->Attribute("blue")) {
            std::string blueStr = element->Attribute("blue");
            double blue = string2Double(blueStr);
            pElement->blue = blue;
          }

          if (element->Attribute("valence")) {
            std::string valenceStr = element->Attribute("valence");
            int valence = string2Int(valenceStr);
            pElement->valence = valence;
          }

          if (element->Attribute("filledShell")) {
            std::string filledShellStr = element->Attribute("filledShell");
            int filledShell = string2Int(filledShellStr);
            pElement->filledShell = filledShell;
          }

          if (element->Attribute("covalentRadius")) {
            std::string covalentRadiusStr = element->Attribute("covalentRadius");
            double covalentRadius = string2Double(covalentRadiusStr);
            pElement->covalentRadius = covalentRadius;
            //std::cout << pElement->symbol << " " << pElement->covalentRadius << std::endl;
          }

          if (element->Attribute("paulingEN")) {
            std::string paulingENStr = element->Attribute("paulingEN");
            double paulingEN = string2Double(paulingENStr);
            pElement->paulingEN = paulingEN;
          }

          if (element->Attribute("seHams")) {
            std::string seHamsStr = element->Attribute("seHams");
            std::vector<std::string> vsSEHams;
            splitString(seHamsStr, " ", vsSEHams, 0);
            for (unsigned int x = 0; x < vsSEHams.size(); x++) {
              pElement->seHams.push_back(vsSEHams[x]);
            }
          }
          //count2++;
        }
        //std::cout << " count2 " << count2 << std::endl;
      }
      //std::cout << "Top level nodes, using First / Next. " << count << std::endl;
    }
    return 0;
}
#endif // USE_TINYXML

// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses element xml file using Trolltech's Qt
// ------------------------------------------------------------
#ifdef USE_QT
int elementParser::Read(std::string fileName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("elementParser", errMessage, INFO);

    // Read elements file using DOM into memory
    QDomDocument doc("mydocument");
    QFile file(qFileName);
    if (!file.open(QIODevice::ReadOnly)) {
      errMessage = " Reading " + fileName;
      errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
      return 1;
    }
    if (!doc.setContent(&file)) {
      errMessage = " Reading " + fileName;
      errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
      file.close();
      return 1;
    }
    file.close();

    // Parse the contents of the elements file
    QString att;
    bool bConvertOk = true;
    QDomElement docElem = doc.documentElement();
    if (docElem.tagName() == "elementlist") {

      QDomNode elemNode = docElem.firstChild();
      while (!elemNode.isNull()) {
        if (elemNode.nodeName() == "element" and elemNode.hasAttributes()) {
          QDomElement elem = elemNode.toElement();
          if (!elem.isNull()) {
            pElement = pElements->addElement();
            //std::cout << elem.tagName().toStdString() << std::endl;

            if (elem.hasAttribute("symbol")) {
              att = elem.attribute("symbol");
              //std::cout << att.toStdString() << std::endl;
              pElements->setElementName(att.toStdString());
            }

            if (elem.hasAttribute("number")) {
              att = elem.attribute("number");
              int number = att.toInt(&bConvertOk, 10);
              if (bConvertOk) {
                pElement->number = number;
              }
              else {
                errMessage = " Error converting string to integer ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser" << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("name")) {
              att = elem.attribute("name");
              pElement->name = att.toStdString();
            }

            if (elem.hasAttribute("mass")) {
              att = elem.attribute("mass");
              double mass = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->mass = mass;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("group")) {
              att = elem.attribute("group");
              int group = att.toInt(&bConvertOk, 10);
              if (bConvertOk) {
                pElement->group = group;
              }
              else {
                errMessage = " Error converting string to integer ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("period")) {
              att = elem.attribute("period");
              int period = att.toInt(&bConvertOk, 10);
              if (bConvertOk) {
                pElement->period = period;
              }
              else {
                errMessage = " Error converting string to integer ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("vdWRadius")) {
              att = elem.attribute("vdWRadius");
              double vdWRadius = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->vdWRadius = vdWRadius;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("red")) {
              att = elem.attribute("red");
              double red = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->red = red;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("green")) {
              att = elem.attribute("green");
              double green = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->green = green;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("blue")) {
              att = elem.attribute("blue");
              double blue = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->blue = blue;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("valence")) {
              att = elem.attribute("valence");
              int valence = att.toInt(&bConvertOk, 10);
              if (bConvertOk) {
                pElement->valence = valence;
              }
              else {
                errMessage = " Error converting string to integer ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("filledShell")) {
              att = elem.attribute("filledShell");
              int filledShell = att.toInt(&bConvertOk, 10);
              if (bConvertOk) {
                pElement->filledShell = filledShell;
              }
              else {
                errMessage = " Error converting string to integer ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("covalentRadius")) {
              att = elem.attribute("covalentRadius");
              double covalentRadius = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->covalentRadius = covalentRadius;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("paulingEN")) {
              att = elem.attribute("paulingEN");
              double paulingEN = att.toDouble(&bConvertOk);
              if (bConvertOk) {
                pElement->paulingEN = paulingEN;
              }
              else {
                errMessage = " Error converting string to double ";
                errorLogger.throwError("elementParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss << "elementParser " << errMessage;
                throw MTKException(ss.str());
              }
            }

            if (elem.hasAttribute("seHams")) {
              att = elem.attribute("seHams");
              std::string sSEHams = att.toStdString();
              std::vector<std::string> vsSEHams;
              splitString(sSEHams, " ", vsSEHams, 0);
              for (unsigned int x = 0; x < vsSEHams.size(); x++) {
                pElement->seHams.push_back(vsSEHams[x]);
              }
            }
          }
          elemNode = elemNode.nextSibling();
        }
      }
    }
    return 0;
}
#endif // USE_QT

// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parsers element xml file using xerces-c
// ============================================================
#ifdef USE_XERCES
int elementParser::Read(std::string fileName)
{
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("elementParser", errMessage, INFO);

    setError(0);
    if (!fileExists(fileName)) {
      setError(1);
      std::string errorMessage = "  elementParser::Read Error, Can't Find " + fileName;
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

    //ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
    //parser->setErrorHandler(errHandler);
    //MyDOMErrorHandler* errHandler = new MyDOMErrorHandler();
    //parser->setErrorHandler(errHandler);

//    try {
       doc = parser->parseURI(xmlFile);
       if (doc == NULL) {
         throw MTKException("elementParser: parse effort turned up nothing for "+fileName+".");
       }

       elementsFiller((DOMNode*)doc->getDocumentElement());
/*
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
*/
    //delete xmlFile;
    parser->release();
    //delete errHandler;
    XMLPlatformUtils::Terminate();

    return 0;
}

// ============================================================
// Function : elementsfiller
// ------------------------------------------------------------
// Parses basis set xml files using xercesc
// Given the rootnode ie the project node
// ------------------------------------------------------------
void elementParser::elementsFiller(DOMNode *rootnode) {
    char *name = XMLString::transcode(rootnode->getNodeName());

    if (std::string(name) == std::string("elementlist")) {
      elementFiller(rootnode);
    }
    delete name;
}

// ============================================================
// Function : elementFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void elementParser::elementFiller(DOMNode *rootnode) {
    DOMNode *child;
    for (child = rootnode->getFirstChild(); child != 0; child = child->getNextSibling()) {
      if (child->getNodeType() == 1) {
        if ((std::string)(XC(child->getNodeName())) == "element") {
          pElement = pElements->addElement();
          if (child->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = child->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("number")) {
                char *cNumber = XMLString::transcode(pAttributeNode->getNodeValue());
                pElement->number = atoi(cNumber);
                delete cNumber;
              }
              if (std::string(name) == std::string("name")) {
                char *cName = XMLString::transcode(pAttributeNode->getNodeValue());
                pElement->name = std::string(cName);
                delete cName;
              }
              if (std::string(name) == std::string("symbol")) {
                char *cSymbol = XMLString::transcode(pAttributeNode->getNodeValue());
                pElements->setElementName(std::string(cSymbol));
                delete cSymbol;
              }
              if (std::string(name) == std::string("mass")) {
                double dMass = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->mass = dMass;
              }
              if (std::string(name) == std::string("group")) {
                char *cGroup = XMLString::transcode(pAttributeNode->getNodeValue());
                pElement->group = atoi(cGroup);
                delete cGroup;
              }
              if (std::string(name) == std::string("period")) {
                char *cPeriod = XMLString::transcode(pAttributeNode->getNodeValue());
                pElement->period = atoi(cPeriod);
                delete cPeriod;
              }
              if (std::string(name) == std::string("red")) {
                double dRed = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->red = dRed;
              }
              if (std::string(name) == std::string("green")) {
                double dGreen = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->green = dGreen;
              }
              if (std::string(name) == std::string("blue")) {
                double dBlue = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->blue = dBlue;
              }
              if (std::string(name) == std::string("valence")) {
                char *cValence = XMLString::transcode(pAttributeNode->getNodeValue());
                pElement->valence = atoi(cValence);
                delete cValence;
              }
              if (std::string(name) == std::string("filledShell")) {
                char *cFS = XMLString::transcode(pAttributeNode->getNodeValue());
                pElement->filledShell = atoi(cFS);
                delete cFS;
              }
              if (std::string(name) == std::string("covalentRadius")) {
                double dCR = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->covalentRadius = dCR;
              }
              if (std::string(name) == std::string("vdWRadius")) {
                double dVR = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->vdWRadius = dVR;
              }
              if (std::string(name) == std::string("paulingEN")) {
                double dPEN = strtod(XC(pAttributeNode->getNodeValue()), 0);
                pElement->paulingEN = dPEN;
              }
              if (std::string(name) == std::string("seHams")) {
                char *cSEHams= XMLString::transcode(pAttributeNode->getNodeValue());
                std::string sSEHams = std::string(cSEHams);
                std::vector<std::string> vsSEHams;
                splitString(sSEHams, " ", vsSEHams, 0);
                for (unsigned int x = 0; x < vsSEHams.size(); x++) {
                  pElement->seHams.push_back(vsSEHams[x]);
                }
                delete cSEHams;
              }
              delete name;
            }
          }
        }
        else {
        }
      }
      else {
      }
    }
}
#endif // USE_XERCES

} // MTKpp namespace
