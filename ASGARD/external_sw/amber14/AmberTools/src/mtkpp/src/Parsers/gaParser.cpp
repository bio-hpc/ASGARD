/*!
   \file gaParser.cpp
   \brief Parses GA xml files using xercesc
   \author Martin Peters

   $Date: 2010/07/22 11:57:13 $
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
#include <sstream>

#include "gaParser.h"
#include "GA/gaWorld.h"
#include "GA/gaRegion.h"
#include "GA/gaPopulation.h"
#include "GA/gaIndividual.h"
#include "GA/gaChromosome.h"
#include "GA/gaGene.h"

#include "StringManip.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : gaParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
gaParser::gaParser(gaWorld *c):myWorld(c)
{
    bWorldFileRead = false;
}

// ============================================================
// Function : ~gaParser()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
gaParser::~gaParser() {}

// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses gaWorld and input xml files using xercesc
// ------------------------------------------------------------
int gaParser::Read(std::string worldFileName, std::string inputFileName)
{
    int success = this->ReadWorld(worldFileName);
    if (success < 0) return success;
    success = this->ReadInput(inputFileName);
    return success;
}

#ifdef USE_QT
// ============================================================
// Function : ReadWorld
// ------------------------------------------------------------
// Parses gaWorld xml files using Qt
// ------------------------------------------------------------
int gaParser::ReadWorld(std::string fileName)
{
    return 0;
}

// ============================================================
// Function : ReadInput
// ------------------------------------------------------------
// Parses gaInput xml files using Qt
// ------------------------------------------------------------
int gaParser::ReadInput(std::string fileName)
{
    return 0;
}
#endif // USE_QT

#ifdef USE_TINYXML
// ============================================================
// Function : ReadWorld
// ------------------------------------------------------------
// Parses gaWorld xml files using Qt
// ------------------------------------------------------------
int gaParser::ReadWorld(std::string fileName)
{
    return 0;
}

// ============================================================
// Function : ReadInput
// ------------------------------------------------------------
// Parses gaInput xml files using Qt
// ------------------------------------------------------------
int gaParser::ReadInput(std::string fileName)
{
    return 0;
}
#endif // USE_TINYXML

#ifdef USE_XERCES
// ============================================================
// ===                                                      ===
// ===              g a W o r l d   R e a d e r             ===
// ===                                                      ===
// ============================================================

// ============================================================
// Function : ReadWorld
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// ------------------------------------------------------------
int gaParser::ReadWorld(std::string fileName)
{
#ifdef DEBUG
    std::cout << "gaParser::ReadWorld" << std::endl;
#endif
    bWorldFileRead = true;

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

    try {
       doc = parser->parseURI(xmlFile);
       gaWorldFiller((DOMNode*)doc->getDocumentElement());
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
    XMLPlatformUtils::Terminate();

    return 0;
}

// ============================================================
// Function : gaWorldFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaWorldFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaWorldFiller" << std::endl;
#endif

    char *name = XMLString::transcode(rootnode->getNodeName());
    if (std::string(name) == std::string("world")) {
      if (rootnode->hasAttributes()) {
        // get all the attributes of the node
        DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
        int nSize = pAttributes->getLength();
        for (int i = 0; i<nSize; ++i) {
          DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
          // get attribute name
          name = XMLString::transcode(pAttributeNode->getName());
          if (std::string(name) == std::string("name")) {
            char *worldName = XMLString::transcode(pAttributeNode->getNodeValue());
            myWorld->setName(std::string(worldName));
            delete worldName;
          }
        }
      }
      gaRegionFiller(rootnode);
    }
    delete name;
}

// ============================================================
// Function : gaRegionFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaRegionFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaRegionFiller" << std::endl;
#endif

    DOMNode *c;

    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "region") {
          myRegion = myWorld->addRegion();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *regionName = XMLString::transcode(pAttributeNode->getNodeValue());
                myRegion->setName(std::string(regionName));
                delete regionName;
              }
              if (std::string(name) == std::string("id")) {
                char *id = XMLString::transcode(pAttributeNode->getNodeValue());
                myRegion->setId(atoi(id));
                delete id;
              }
              delete name;
            }
          }
          gaPopulationFiller(c);
        }
        else{
        }
      }
      else{
      }
    }
}

// ============================================================
// Function : gaPopulationFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaPopulationFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaPopulationFiller" << std::endl;
#endif

    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "population") {
          myPopulation = myRegion->addPopulation();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *popName = XMLString::transcode(pAttributeNode->getNodeValue());
                myPopulation->setName(std::string(popName));
                delete popName;
              }
              if (std::string(name) == std::string("id")) {
                char *id = XMLString::transcode(pAttributeNode->getNodeValue());
                myPopulation->setId(atoi(id));
                delete id;
              }
              delete name;
            }
          }
          gaIndividualFiller(c);
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : gaIndividualFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaIndividualFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaIndividualFiller" << std::endl;
#endif

    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName()))=="individual") {
          myIndividual = myPopulation->addIndividual();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *indName = XMLString::transcode(pAttributeNode->getNodeValue());
                myIndividual->setName(std::string(indName));
                delete indName;
              }
              //if (string(name) == string("id")){
              //  char *id = XMLString::transcode(pAttributeNode->getNodeValue());
              //  myIndividual->setId(atoi(id));
              //}
              delete name;
            }
          }
          gaChromosomeFiller(c);
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : gaChromosomeFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaChromosomeFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaChromosomeFiller" << std::endl;
#endif

    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType()==1){
        if ((std::string)(XC(c->getNodeName())) == "chromosome") {
          myChromosome = myIndividual->addChromosome();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *chrName = XMLString::transcode(pAttributeNode->getNodeValue());
                myChromosome->setName(std::string(chrName));
                delete chrName;
              }
              //if (string(name) == string("id")){
              //  char *id = XMLString::transcode(pAttributeNode->getNodeValue());
              //  myChromosome->setId(atoi(id));
              //}
              delete name;
            }
          }
          gaGeneFiller(c);
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : gaGeneFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaGeneFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaGeneFiller" << std::endl;
#endif

    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "gene") {
          myGene = myChromosome->addGene();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *genName = XMLString::transcode(pAttributeNode->getNodeValue());
                myGene->setName(std::string(genName));
                delete genName;
              }
              //if (string(name) == string("id")){
              //  char *id = XMLString::transcode(pAttributeNode->getNodeValue());
              //  myGene->setId(atoi(id));
              //}
              delete name;
            }
          }
          char *gene = XMLString::transcode(c->getTextContent());
          std::string sGene = std::string(gene);
          sGene = removeCharacter(sGene,'[');
          sGene = removeCharacter(sGene,']');
          std::vector<std::string> bits;
          splitString(sGene, ",", bits, 0);
          for (unsigned int b = 0; b < bits.size(); b++) {
            myGene->addBit( strtod(bits[b].c_str(), 0) );
          }
          delete gene;
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// ===                                                      ===
// ===              g a I n p u t   R e a d e r             ===
// ===                                                      ===
// ============================================================

// ============================================================
// Function : ReadInput
// ------------------------------------------------------------
// Parses gaInput xml files using xercesc
// ------------------------------------------------------------
int gaParser::ReadInput(std::string fileName)
{
#ifdef DEBUG
    std::cout << "\n\n\n gaParser::ReadInput" << std::endl;
#endif
    if (!bWorldFileRead) {
      std::cout << " Run time error.  Developer must read in world xml file first " << std::endl;
      //std::cout << " Exiting " << std::endl;
      throw MTKException(" Run time error.  Developer must read in world xml file first ");
      //exit(0);
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

    try {
       doc = parser->parseURI(xmlFile);
       gaWorldInputFiller((DOMNode*)doc->getDocumentElement());
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
    XMLPlatformUtils::Terminate();

    return 0;
}

// ============================================================
// Function : gaWorldInputFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaWorldInputFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaWorldInputFiller" << std::endl;
#endif

    char *name = XMLString::transcode(rootnode->getNodeName());

    std::string worldName = "";
    int chrPerInd = 1;
    std::string genePerChr;
    std::vector<std::string> vGenePerChr;
    std::vector<int> vdGenePerChr;

    std::string geneSizes;
    std::vector<std::string> vGeneSizes;
    std::vector<int> vdGeneSizes;

    std::string outputFileName = "";
    std::string restartFileName = "";
    std::string convergFileName = "";
    int levelOfOutput = 1;

    if (std::string(name) == std::string("world")) {
      if (rootnode->hasAttributes()) {
        // get all the attributes of the node
        DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
        int nSize = pAttributes->getLength();
        for (int i = 0; i < nSize; ++i) {
          DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
          // get attribute name
          name = XMLString::transcode(pAttributeNode->getName());
          if (std::string(name) == std::string("name")) {
            char *cWorldName = XMLString::transcode(pAttributeNode->getNodeValue());
            worldName = std::string(cWorldName);
            delete cWorldName;
          }
          if (std::string(name) == std::string("chrPerInd")) {
            char *cChrPerInd = XMLString::transcode(pAttributeNode->getNodeValue());
            chrPerInd = atoi(cChrPerInd);
            delete cChrPerInd;
          }
          if (std::string(name) == std::string("genePerChr")) {
            char *cGenePerChr = XMLString::transcode(pAttributeNode->getNodeValue());
            genePerChr = std::string(cGenePerChr);
            splitString(genePerChr, " ", vGenePerChr, 0);
            for (unsigned int j = 0; j < vGenePerChr.size(); j++) {
              vdGenePerChr.push_back(atoi(vGenePerChr[j].c_str()));
            }
            delete cGenePerChr;
          }
          if (std::string(name) == std::string("geneSizes")) {
            char *cGeneSizes = XMLString::transcode(pAttributeNode->getNodeValue());
            geneSizes = std::string(cGeneSizes);
            splitString(geneSizes, " ", vGeneSizes, 0);
            for (unsigned int j = 0; j < vGeneSizes.size(); j++) {
              vdGeneSizes.push_back(atoi(vGeneSizes[j].c_str()));
            }
            delete cGeneSizes;
          }
          if (std::string(name) == std::string("outputFileName")) {
            char *cOutputFileName = XMLString::transcode(pAttributeNode->getNodeValue());
            outputFileName = std::string(cOutputFileName);
            delete cOutputFileName;
          }
          if (std::string(name) == std::string("restartFileName")) {
            char *cRestartFileName = XMLString::transcode(pAttributeNode->getNodeValue());
            restartFileName = std::string(cRestartFileName);
            delete cRestartFileName;
          }
          if (std::string(name) == std::string("convergFileName")) {
            char *cConvergFileName = XMLString::transcode(pAttributeNode->getNodeValue());
            convergFileName = std::string(cConvergFileName);
            delete cConvergFileName;
          }
          if (std::string(name) == std::string("levelOfOutput")) {
            char *cLevelOfOutput = XMLString::transcode(pAttributeNode->getNodeValue());
            levelOfOutput = atoi(cLevelOfOutput);
            delete cLevelOfOutput;
          }
        }
      }
      if (worldName == myWorld->getName()) {
#ifdef DEBUG
    std::cout << " chrPerInd " << chrPerInd << std::endl;
    std::cout << " outputFileName " << outputFileName << std::endl;
    std::cout << " restartFileName " << restartFileName << std::endl;
    std::cout << " convergFileName " << convergFileName << std::endl;
    std::cout << " levelOfOutput " << levelOfOutput << std::endl;

#endif
        myWorld->setup(chrPerInd, vdGenePerChr, vdGeneSizes, outputFileName,
                 restartFileName, convergFileName, levelOfOutput);
        gaRegionInputFiller(rootnode);
      }
      else {
        std::cout << " gaWorld names do not match ... exiting " << std::endl;
        //exit(0);
        throw MTKException(" gaWorld names do not match ... ");
      }
    }
    delete name;
}

// ============================================================
// Function : gaRegionInputFiller
// ------------------------------------------------------------
// Parses gaWorld xml files using xercesc
// Given the rootnode
// ------------------------------------------------------------
void gaParser::gaRegionInputFiller(DOMNode *rootnode)
{
#ifdef DEBUG
    std::cout << "gaParser::gaRegionInputFiller" << std::endl;
#endif

    DOMNode *c;
    std::string regionName;
    int maxInds = 0;
    int seed = 0;
    int maxGens = 0;
    int elitism = 0;
    double selectionPressure = 0.0;
    int curPop = 0;
    int popKeep = 1;
    int nChild = 1;
    double chreDiff = 0.01;
    double pKeep = 0.5;
    double pCrossover = 0.2;
    double pMutate = 0.18;
    double pAverage = 0.12;
    std::string selection = "semi-random";
    std::string crossover = "single-gene";
    std::string mutate = "single-gene";
    std::string average = "single-gene";
    std::string funcDir = "maximize";
    std::string maxParameters = "";
    std::string minParameters = "";
    std::string stepSize = "";
    std::vector<std::string> vMaxParameters;
    std::vector<double> vdMaxParameters;
    std::vector<std::string> vMinParameters;
    std::vector<double> vdMinParameters;
    std::vector<std::string> vStepSize;
    std::vector<double> vdStepSize;

    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "region") {
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")){
                char *regName = XMLString::transcode(pAttributeNode->getNodeValue());
                regionName = std::string(regName);
                delete regName;
              }
              //if (string(name) == string("id")){
              //  char *cId = XMLString::transcode(pAttributeNode->getNodeValue());
              //  id = atoi(cId);
              //}
              if (std::string(name) == std::string("maxInds")) {
                char *cMaxInds = XMLString::transcode(pAttributeNode->getNodeValue());
                maxInds = atoi(cMaxInds);
                delete cMaxInds;
              }
              if (std::string(name) == std::string("seed")) {
                char *cSeed = XMLString::transcode(pAttributeNode->getNodeValue());
                seed = atoi(cSeed);
                delete cSeed;
              }
              if (std::string(name) == std::string("maxGens")) {
                char *cMaxGens = XMLString::transcode(pAttributeNode->getNodeValue());
                maxGens = atoi(cMaxGens);
                delete cMaxGens;
              }
              if (std::string(name) == std::string("elitism")) {
                char *cElitism = XMLString::transcode(pAttributeNode->getNodeValue());
                elitism = atoi(cElitism);
                delete cElitism;
              }
              if (std::string(name) == std::string("selectionPressure")) {
                char *cSelectionPressure = XMLString::transcode(pAttributeNode->getNodeValue());
                selectionPressure = strtod(cSelectionPressure, 0);
                delete cSelectionPressure;
              }
              if (std::string(name) == std::string("curPop")) {
                char *cCurPop = XMLString::transcode(pAttributeNode->getNodeValue());
                curPop = atoi(cCurPop);
                delete cCurPop;
              }
              if (std::string(name) == std::string("popKeep")) {
                char *cPopKeep = XMLString::transcode(pAttributeNode->getNodeValue());
                popKeep = atoi(cPopKeep);
                delete cPopKeep;
              }
              if (std::string(name) == std::string("nChild")) {
                char *cnChild = XMLString::transcode(pAttributeNode->getNodeValue());
                nChild = atoi(cnChild);
                delete cnChild;
              }
              if (std::string(name) == std::string("chreDiff")) {
                char *cChreDiff = XMLString::transcode(pAttributeNode->getNodeValue());
                chreDiff = strtod(cChreDiff, 0);
                delete cChreDiff;
              }
              if (std::string(name) == std::string("pKeep")) {
                char *cpKeep = XMLString::transcode(pAttributeNode->getNodeValue());
                pKeep = strtod(cpKeep, 0);
                delete cpKeep;
              }
              if (std::string(name) == std::string("pCrossover")) {
                char *cpCrossover = XMLString::transcode(pAttributeNode->getNodeValue());
                pCrossover = strtod(cpCrossover, 0);
                delete cpCrossover;
              }
              if (std::string(name) == std::string("pMutate")) {
                char *cpMutate = XMLString::transcode(pAttributeNode->getNodeValue());
                pMutate = strtod(cpMutate, 0);
                delete cpMutate;
              }
              if (std::string(name) == std::string("pAverage")) {
                char *cpAverage = XMLString::transcode(pAttributeNode->getNodeValue());
                pAverage = strtod(cpAverage, 0);
                delete cpAverage;
              }
              if (std::string(name) == std::string("selection")) {
                char *cSelection = XMLString::transcode(pAttributeNode->getNodeValue());
                selection = std::string(cSelection);
                delete cSelection;
              }
              if (std::string(name) == std::string("crossover")) {
                char *cCrossover = XMLString::transcode(pAttributeNode->getNodeValue());
                crossover = std::string(cCrossover);
                delete cCrossover;
              }
              if (std::string(name) == std::string("mutate")) {
                char *cMutate = XMLString::transcode(pAttributeNode->getNodeValue());
                mutate = std::string(cMutate);
                delete cMutate;
              }
              if (std::string(name) == std::string("average")) {
                char *cAverage = XMLString::transcode(pAttributeNode->getNodeValue());
                average = std::string(cAverage);
                delete cAverage;
              }
              if (std::string(name) == std::string("funcDir")) {
                char *cFuncDir = XMLString::transcode(pAttributeNode->getNodeValue());
                funcDir = std::string(cFuncDir);
                delete cFuncDir;
              }
              if (std::string(name) == std::string("maxParameters")) {
                char *cMaxParameters = XMLString::transcode(pAttributeNode->getNodeValue());
                maxParameters = std::string(cMaxParameters);
                splitString(maxParameters, " ", vMaxParameters, 0);
                for (unsigned int j = 0; j < vMaxParameters.size(); j++) {
                  vdMaxParameters.push_back(strtod(vMaxParameters[j].c_str(), 0));
                }
                delete cMaxParameters;
              }
              if (std::string(name) == std::string("minParameters")) {
                char *cMinParameters = XMLString::transcode(pAttributeNode->getNodeValue());
                minParameters = std::string(cMinParameters);
                splitString(minParameters, " ", vMinParameters, 0);
                for (unsigned int j = 0; j < vMinParameters.size(); j++) {
                  vdMinParameters.push_back(strtod(vMinParameters[j].c_str(), 0));
                }
                delete cMinParameters;
              }
              if (std::string(name) == std::string("stepSize")) {
                char *cStepSize = XMLString::transcode(pAttributeNode->getNodeValue());
                stepSize = std::string(cStepSize);
                splitString(stepSize, " ", vStepSize, 0);
                for (unsigned int j = 0; j < vStepSize.size(); j++) {
                  vdStepSize.push_back(strtod(vStepSize[j].c_str(), 0));
                }
                delete cStepSize;
              }
            }
          }

          myRegion = myWorld->getRegion(regionName);
          if (myRegion != 0) {
            for (unsigned int i = 0; i < vdMinParameters.size(); i++) {
              if (vdMinParameters[i] > vdMaxParameters[i]) { 
                //std::cout << " gaParser:: error user defined minParameter: " << i+1
                //          << "  larger than maxParameter ... exiting " << std::endl;
                //exit(0);
                std::stringstream ss;
                ss << " gaParser:: error user defined minParameter: " << i+1
                   << "  larger than maxParameter ... exiting " << std::endl;
                std::cout << ss.str();
                throw MTKException(ss.str());
              }
            }
#ifdef DEBUG
            std::cout << " maxInds " << maxInds << std::endl;
            std::cout << " seed " << seed << std::endl;
            std::cout << " maxGens " << maxGens << std::endl;
            std::cout << " elitism " << elitism << std::endl;
            std::cout << " selectionPressure " << selectionPressure << std::endl;
            std::cout << " curPop " << curPop << std::endl;
            std::cout << " popKeep " << popKeep << std::endl;
            std::cout << " chreDiff " << chreDiff << std::endl;
            std::cout << " pKeep " << pKeep << std::endl;
            std::cout << " pCrossover " << pCrossover << std::endl;
            std::cout << " pMutate " << pMutate << std::endl;
            std::cout << " pAverage " << pAverage << std::endl;
            std::cout << " selection " << selection << std::endl;
            std::cout << " crossover " << crossover << std::endl;
            std::cout << " mutate " << mutate << std::endl;
            std::cout << " average " << average << std::endl;
            std::cout << " funcDir " << funcDir << std::endl;
#endif

            myRegion->setup(maxInds, seed, maxGens, elitism, selectionPressure, curPop, popKeep,
                            nChild, chreDiff, pKeep, pCrossover, pMutate, pAverage,
                            selection, crossover, mutate, average, funcDir,
                            vdMaxParameters, vdMinParameters, vdStepSize);
          }
          else {
            std::cout << " region is null, check names " << std::endl;
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
