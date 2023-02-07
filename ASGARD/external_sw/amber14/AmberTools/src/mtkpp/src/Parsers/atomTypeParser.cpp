/*!
   \file atomTypeParser.cpp
   \brief atom type parser
   \author Martin Peters

   Parses atom type xml files using xercesc

   $Date: 2010/03/29 20:39:34 $
   $Revision: 1.7 $

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

#include "atomTypeParser.h"
#include "Molecule/atomType.h"

namespace MTKpp
{

// ============================================================
// Function : atomTypeParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
atomTypeParser::atomTypeParser(atomTypes *c):pAtomTypes(c) {}

// ============================================================
// Function : atomTypeParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
atomTypeParser::~atomTypeParser() {}

// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses basis set xml files using xercesc
// ------------------------------------------------------------
int atomTypeParser::Read(std::string fileName)
{
#ifdef USE_XERCES
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
    //DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);
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
    DOMDocument *doc = 0;

//    ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
//    parser->setErrorHandler(errHandler);
    //MyDOMErrorHandler* errHandler = new MyDOMErrorHandler();
    //parser->setErrorHandler(errHandler);

    try {
       doc = parser->parseURI(xmlFile);
       atomTypesFiller((DOMNode*)doc->getDocumentElement());
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

    parser->release();
    //delete errHandler;
#endif // USE_XERCES
    return 0;
}

#ifdef USE_XERCES
// ============================================================
// Function : atomTypesFiller
// ------------------------------------------------------------
// Parses basis set xml files using xercesc
// Given the rootnode ie the project node
// ------------------------------------------------------------
void atomTypeParser::atomTypesFiller(DOMNode *rootnode){
    char *name = XMLString::transcode(rootnode->getNodeName());

    if (std::string(name) == std::string("atomtypes")){
       if(rootnode->hasAttributes()) {
          // get all the attributes of the node
          DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
          int nSize = pAttributes->getLength();
          for(int i=0; i<nSize; ++i) {
             DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
             // get attribute name
             char *name = XMLString::transcode(pAttributeNode->getName());
             if (std::string(name) == std::string("name")){
                char *atomTypeName = XMLString::transcode(pAttributeNode->getNodeValue());
                pAtomTypes->setName(std::string(atomTypeName));
             }
          }
       }
       atomTypeFiller(rootnode);
    }
}

// ============================================================
// Function : atomTypefiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void atomTypeParser::atomTypeFiller(DOMNode *rootnode){
    DOMNode *c;
    for (c=rootnode->getFirstChild(); c!=0; c=c->getNextSibling()){
       if(c->getNodeType()==1){
          if((std::string)(XC(c->getNodeName()))=="atomtype"){
             pAtomType = pAtomTypes->addAtomType();
             if (c->hasAttributes()){
                // get all the attributes of the node
                DOMNamedNodeMap *pAttributes = c->getAttributes();
                int nSize = pAttributes->getLength();
                for(int i=0;i<nSize;++i) {
                   DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
                   // get attribute name
                   char *name = XMLString::transcode(pAttributeNode->getName());
                   if (std::string(name) == std::string("name")){
                      char *atomTypeName = XMLString::transcode(pAttributeNode->getNodeValue());
                      pAtomType->name = std::string(atomTypeName);
                      pAtomTypes->setAtomTypeName(std::string(atomTypeName));
                   }
                   if (std::string(name) == std::string("element")){
                      char *atomTypeElement = XMLString::transcode(pAttributeNode->getNodeValue());
                      pAtomType->element = std::string(atomTypeElement);
                   }
                   if (std::string(name) == std::string("hybridization")){
                      char *atomTypeHyb = XMLString::transcode(pAttributeNode->getNodeValue());
                      pAtomType->hybridization = std::string(atomTypeHyb);
                   }
                   if (std::string(name) == std::string("description")){
                      char *atomTypeDes = XMLString::transcode(pAttributeNode->getNodeValue());
                      pAtomType->description = std::string(atomTypeDes);
                   }
                }
             }
             parameterFiller(c);
          }
          else{
          }
       }
       else{
       }
    }
}

// ============================================================
// Function : parameterFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void atomTypeParser::parameterFiller(DOMNode *rootnode){
    DOMNode *child;
    double value = 0.0;
    char *param = 0;
    for (child=rootnode->getFirstChild(); child!=0; child=child->getNextSibling()){
       if(child->getNodeType()==1){
          if((std::string)(XC(child->getNodeName()))=="parameter"){
             if (child->hasAttributes()){
                // get all the attributes of the node
                DOMNamedNodeMap *pAttributes = child->getAttributes();
                int nSize = pAttributes->getLength();
                for(int i=0;i<nSize;++i) {
                   DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
                   // get attribute name
                   char *name = XMLString::transcode(pAttributeNode->getName());
                   if (std::string(name) == std::string("name")){
                      param = XMLString::transcode(pAttributeNode->getNodeValue());
                   }
                   if (std::string(name) == std::string("value")){
                      value = strtod(XC(pAttributeNode->getNodeValue()), 0);
                   }
                }
                if (param == std::string("rvalue")){
                   pAtomType->rvalue = value;
                }
                else if (param == std::string("evalue")){
                   pAtomType->evalue = value;
                }
             }
          }
          else{
          }
       }
       else{
       }
    }
}
#endif // USE_XERCES

} // MTKpp namespace

