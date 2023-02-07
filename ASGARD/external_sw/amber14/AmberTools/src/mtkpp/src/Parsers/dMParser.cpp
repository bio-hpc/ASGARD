/*!
   \file dMParser.cpp
   \brief matrix xml parser
   \author Martin Peters

   $Date: 2010/04/29 19:06:19 $
   $Revision: 1.9 $

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
#include "dMParser.h"
#include "Statistics/sheet.h"
#include "Statistics/table.h"
#include "StringManip.h"

#include "Log/errorHandler.h"
#include "config.h"

#ifdef USE_XERCES
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
XERCES_CPP_NAMESPACE_USE
#endif // USE_XERCES

namespace MTKpp
{

// ============================================================
// Class : dMParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
dMParser::dMParser() {}

// ============================================================
// Function : ~dMParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
dMParser::~dMParser() {}

// ============================================================
// Function : import
// ------------------------------------------------------------
// 
// ============================================================
int dMParser::import(sheet* s, std::string fileName)
{
    mySheet = s;
    if (!mySheet) return 1;

    std::ifstream ifile;
    ifile.open(fileName.c_str());

    if (!ifile) {
      std::cout << "\nUNABLE TO OPEN FILE"
           << "\nFILENAME = " << fileName
           << "\nEXITING...\n" << std::endl;
      return 1;
    }

    int end   = fileName.length();
    int slash = fileName.find_last_of("/");
    std::string file_name = fileName.substr(slash+1,(end-slash-5));

    std::string fileline;

    // 1st line
    getline(ifile,fileline);
    std::vector<std::string> splitLine1;
    splitString(fileline, " ", splitLine1, 0);
    if (splitLine1.size() == 3) {
      myDoubleTable = mySheet->addTable();
      if (!myDoubleTable) return 1;

      myDoubleTable->setName(splitLine1[0]);
      myDoubleTable->setSizes(atoi(splitLine1[1].c_str()), atoi(splitLine1[2].c_str()));
    }

    // 2nd line
    getline(ifile,fileline);
    std::vector<std::string> splitLine2;
    splitString(fileline, " ", splitLine2, 0);
    for (unsigned int i = 0; i < splitLine2.size(); i++) {
      myDoubleTable->setColumnLabel(i, splitLine2[i]);
    }

    int rowIndex = 0;
    while (ifile) {
      getline(ifile,fileline);
      if (fileline == "") break;
      std::vector<std::string> splitLine3;
      splitString(fileline, " ", splitLine3, 0);
      myDoubleTable->setRowLabel(rowIndex, splitLine3[0]);

      for (unsigned int i = 1; i < splitLine3.size(); i++) {
        double cellValue = strtod(splitLine3[i].c_str(), 0);
        myDoubleTable->setCellValue(rowIndex, i-1, cellValue);
      }
      rowIndex++;
    }
    return 0;
}

    // ---------------------------- //
    // -    QT READ FUNCTIONS     - //
    // ---------------------------- //

#ifdef USE_QT
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses dM xml files using Qt
// ------------------------------------------------------------
int dMParser::read(sheet* s, std::string fileName)
{
    mySheet = s;

    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("dMParser", errMessage, INFO);

    // Read elements file using DOM into memory
    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QFile file(qFileName);
    if (!file.open(QIODevice::ReadOnly)) {
      errMessage = " Reading " + fileName;
      errorLogger.throwError("dMParser", errMessage, MTK_ERROR);
    }

    if (!doc.setContent(&file)) {
      errMessage = " Reading " + fileName;
      errorLogger.throwError("dMParser", errMessage, MTK_ERROR);
      file.close();
    }
    file.close();

    // Parse the contents of the elements file
    QString att;
    QDomElement docElem = doc.documentElement();
    if (docElem.tagName() == "Sheet") {

      if (docElem.hasAttribute("name")) {
        att = docElem.attribute("name");
        mySheet->setName(att.toStdString());
      }

      QDomNode tableNode = docElem.firstChild();
      while (!tableNode.isNull()) {
        if (tableNode.nodeName() == "Table" and tableNode.hasAttributes()) {
          tableFiller(tableNode);
        }
        tableNode = tableNode.nextSibling();
      }
    }
    return 0;
}

// ============================================================
// Function : tableFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::tableFiller(QDomNode g)
{
    QDomElement tableElement = g.toElement();

    bool bName = false;
    bool bnCols = false;
    bool bnRows = false;
    bool bType = false;

    std::string tableName = "";
    int nCols = 0;
    int nRows = 0;
    std::string tableType = "d";

    if (!tableElement.isNull()) {
      QString att;

      if (tableElement.hasAttribute("name")) {
        att = tableElement.attribute("name");
        bName = true;
        tableName = att.toStdString();
      }

      if (tableElement.hasAttribute("num_cols")) {
        att = tableElement.attribute("num_cols");
        bnCols = true;
        nCols = att.toInt();
      }

      if (tableElement.hasAttribute("num_rows")) {
        att = tableElement.attribute("num_rows");
        bnRows = true;
        nRows = att.toInt();
      }

      if (tableElement.hasAttribute("type")) {
        att = tableElement.attribute("type");
        bType = true;
        tableType = att.toStdString();
      }

      if (tableType == "d") {
        myDoubleTable = mySheet->addTable();
        myDoubleTable->setName(tableName);
        myDoubleTable->setType(tableType);
        myDoubleTable->setNumColumns(nCols);
        myDoubleTable->setNumRows(nRows);
        myDoubleTable->setup();
        currentType = "d";
      }
      else if (tableType == "i") {
        myIntTable = mySheet->addIntTable();
        myIntTable->setName(tableName);
        myIntTable->setType(tableType); // should this be doubleTable ??
        myIntTable->setNumColumns(nCols);
        myIntTable->setNumRows(nRows);
        myIntTable->setup();
        currentType = "i";
      }

      QDomNode tNode = g.firstChild();
      while (!tNode.isNull()) {
        if (tNode.nodeName() == "Row" and tNode.hasAttributes()) {
          rowColFiller(tNode);
        }
        if (tNode.nodeName() == "Col" and tNode.hasAttributes()) {
          rowColFiller(tNode);
        }
        tNode = tNode.nextSibling();
      }
    }
}

// ============================================================
// Function : rowColFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::rowColFiller(QDomNode g)
{
    QDomElement rowColElement = g.toElement();

    std::string rowName = "";
    int rowId = -1;
    std::string colName = "";
    int colId = -1;

    if (g.nodeName() == "Row") { // <Row id="0" name="CH2Cl2" >
      if (!rowColElement.isNull()) {
        QString att;

        if (rowColElement.hasAttribute("name")) {
          att = rowColElement.attribute("name");
          rowName = att.toStdString();
        }

        if (rowColElement.hasAttribute("id")) {
          att = rowColElement.attribute("id");
          rowId = att.toInt();
        }

        if ((rowId > -1) and (rowName != "")) {
          if (currentType == "d") {
            myDoubleTable->setRowLabel(rowId, rowName);
          }
          else if (currentType == "i") {
            myIntTable->setRowLabel(rowId, rowName);
          }
        }

        if (currentType == "d") {
          cellFiller(g, rowId);
        }
        else if (currentType == "i") {
          cellIntFiller(g, rowId);
        }
        rowName = "";
        rowId = -1;
      }
    }
    else if (g.nodeName() == "Col") { // <Col id="0" name="LD25" />
      if (!rowColElement.isNull()) {
        QString att;

        if (rowColElement.hasAttribute("name")) {
          att = rowColElement.attribute("name");
          colName = att.toStdString();
        }

        if (rowColElement.hasAttribute("id")) {
          att = rowColElement.attribute("id");
          colId = att.toInt();
        }

        if ((colId > -1) and (colName != "")) {
          if (currentType == "d") {
            myDoubleTable->setColumnLabel(colId, colName);
          }
          else if (currentType == "i") {
            myIntTable->setColumnLabel(colId, colName);
          }
        }
        colName = "";
        colId = -1;
      }
    }
}

// ============================================================
// Function : cellFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::cellFiller(QDomNode g, int rowId)
{
    QDomElement rowColElement = g.toElement();

    int colId = -1;
    double cellValue = 0.0;

    QDomNode tNode = g.firstChild();
    while (!tNode.isNull()) {
      if (tNode.nodeName() == "Cell" and tNode.hasAttributes()) { // <Cell col="0" value="0.96" />
        QDomElement cellElement = tNode.toElement();

        if (!cellElement.isNull()) {
          QString att;

          if (cellElement.hasAttribute("col")) {
            att = cellElement.attribute("col");
            colId = att.toInt();
          }

          if (cellElement.hasAttribute("value")) {
            att = cellElement.attribute("value");
            cellValue = att.toDouble();
          }
          if ((colId > -1)) {
            myDoubleTable->setCellValue(rowId, colId, cellValue);
          }
          cellValue = 0.0;
          colId = -1;
          
        }
      }
      tNode = tNode.nextSibling();
    }
}

// ============================================================
// Function : cellIntFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::cellIntFiller(QDomNode g, int rowId)
{
    QDomElement rowColElement = g.toElement();

    int colId = -1;
    int cellValue = 0;

    QDomNode tNode = g.firstChild();
    while (!tNode.isNull()) {
      if (tNode.nodeName() == "Cell" and tNode.hasAttributes()) { // <Cell col="0" value="0.96" />
        QDomElement cellElement = tNode.toElement();

        if (!cellElement.isNull()) {
          QString att;

          if (cellElement.hasAttribute("col")) {
            att = cellElement.attribute("col");
            colId = att.toInt();
          }

          if (cellElement.hasAttribute("value")) {
            att = cellElement.attribute("value");
            cellValue = att.toInt();
          }
          if ((colId > -1)) {
            myIntTable->setCellValue(rowId, colId, cellValue);
          }
          cellValue = 0.0;
          colId = -1;
        }
      }
      tNode = tNode.nextSibling();
    }
}

#endif // USE_QT

    // ---------------------------------- //
    // -     TINYXML READ FUNCTIONS     - //
    // ---------------------------------- //

#ifdef USE_TINYXML
// ============================================================
// Function : read
// ------------------------------------------------------------
// 
// ============================================================
int dMParser::read(sheet* s, std::string fileName)
{
    mySheet = s;

    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("dMParser", errMessage, INFO);

    TiXmlDocument doc(fileName);
    bool loadOkay = doc.LoadFile();

    if (!loadOkay) {
      //printf ("Could not load parameter file. Error='%s'. Exiting.\n", doc.ErrorDesc() );
      errMessage = " Reading " + fileName;
      errorLogger.throwError("dMParser", errMessage, MTK_ERROR);
      //exit( 1 );
      std::stringstream ss;
      ss << "dMParser error " << errMessage;
      throw MTKException(ss.str());
    }
    else {
      TiXmlNode* node = 0;
      TiXmlElement* element;

      node = doc.RootElement();
      element = node->ToElement();

      if (element->Attribute("name")) {
        mySheet->setName(element->Attribute("name"));
        //std::cout << " sheet " << element->Attribute("name") << std::endl;
      }

      for (TiXmlNode* tableNode = node->FirstChild("Table"); tableNode; tableNode = tableNode->NextSibling("Table")) {
        TiXmlElement* tableElement = tableNode->ToElement();

        bool bName = false;
        bool bnCols = false;
        bool bnRows = false;
        bool bType = false;

        std::string tableName = "";
        int nCols = 0;
        int nRows = 0;
        std::string tableType = "d";
 
        if (tableElement->Attribute("name")) {
          //std::cout << "    table " << tableElement->Attribute("name") << std::endl;
          tableName = tableElement->Attribute("name");
          bName = true;
        }

        if (tableElement->Attribute("num_cols")) {
          nCols = string2Int(tableElement->Attribute("num_cols"));
          bnCols = true;
        }

        if (tableElement->Attribute("num_rows")) {
          nRows = string2Int(tableElement->Attribute("num_rows"));
          bnRows = true;
        }

        if (tableElement->Attribute("type")) {
          tableType = tableElement->Attribute("type");
          bType = true;
        }

        //std::cout << tableName << " " << tableType << " " << nCols << " " << nRows << std::endl;

        if (tableType == "d") {
          myDoubleTable = mySheet->addTable();
          myDoubleTable->setName(tableName);
          myDoubleTable->setType(tableType);
          myDoubleTable->setNumColumns(nCols);
          myDoubleTable->setNumRows(nRows);
          myDoubleTable->setup();
          currentType = "d";
        }
        else if (tableType == "i") {
          myIntTable = mySheet->addIntTable();
          myIntTable->setName(tableName);
          myIntTable->setType(tableType);
          myIntTable->setNumColumns(nCols);
          myIntTable->setNumRows(nRows);
          myIntTable->setup();
          currentType = "i";
        }

        //std::cout << " Doing Cols " << std::endl;
        for (TiXmlNode* colNode = tableNode->FirstChild("Col"); colNode; colNode = colNode->NextSibling("Col")) {
          TiXmlElement* colElement = colNode->ToElement();

          std::string colName = "";
          int colId = -1;
          if (colElement->Attribute("name")) {
            //std::cout << "        col name " << colElement->Attribute("name") << std::endl;
            colName = colElement->Attribute("name");
          }

          if (colElement->Attribute("id")) {
            //std::cout << "        col id " << colElement->Attribute("id") << std::endl;
            colId = string2Int(colElement->Attribute("id"));
          }

          if ((colId > -1) and (colName != "")) {
            if (currentType == "d") {
              myDoubleTable->setColumnLabel(colId, colName);
            }
            else if (currentType == "i") {
              myIntTable->setColumnLabel(colId, colName);
            }
          }
          colName = "";
          colId = -1;
        }

        //std::cout << " Doing Rows " << std::endl;
        for (TiXmlNode* rowNode = tableNode->FirstChild("Row"); rowNode; rowNode = rowNode->NextSibling("Row")) {
          TiXmlElement* rowElement = rowNode->ToElement();

          std::string rowName = "";
          int rowId = -1;

          if (rowElement->Attribute("name")) {
            //std::cout << "        row  name " << rowElement->Attribute("name") << std::endl;
            rowName = rowElement->Attribute("name");
          }

          if (rowElement->Attribute("id")) {
            //std::cout << "        row id " << rowElement->Attribute("id") << std::endl;
            rowId = string2Int(rowElement->Attribute("id"));
          }

          if ((rowId > -1) and (rowName != "")) {
            if (currentType == "d") {
              myDoubleTable->setRowLabel(rowId, rowName);
            }
            else if (currentType == "i") {
              myIntTable->setRowLabel(rowId, rowName);
            }
          }

          //std::cout << " Doing Cells " << std::endl;
          for (TiXmlNode* cellNode = rowNode->FirstChild("Cell"); cellNode; cellNode = cellNode->NextSibling("Cell")) {
            TiXmlElement* cellElement = cellNode->ToElement();

            int colId = -1;
            double cellValue_d = 0.0;
            int cellValue_i = 0;

            if (cellElement->Attribute("col")) {
              //std::cout << "        cell col " << cellElement->Attribute("col") << std::endl;
              colId = string2Int(cellElement->Attribute("col"));
            }

            if (cellElement->Attribute("value")) {

              if (tableType == "d") {
                cellValue_d = string2Double(cellElement->Attribute("value"));
              }
              else if (tableType == "i") {
                cellValue_i = string2Int(cellElement->Attribute("value"));
              }
              //std::cout << "        cell value " << cellElement->Attribute("value") << std::endl;
            }

            if ((colId > -1)) {
              if (tableType == "d") {
                myDoubleTable->setCellValue(rowId, colId, cellValue_d);
              }
              else if (tableType == "i") {
                myIntTable->setCellValue(rowId, colId, cellValue_i);
              }
            }
            cellValue_d = 0.0;
            cellValue_i = 0;
            colId = -1;
          }
          rowName = "";
          rowId = -1;
        }
      }
    }
    return 0;
}
#endif // USE_TINYXML


    // ---------------------------------- //
    // -    XERCES-C READ FUNCTIONS     - //
    // ---------------------------------- //

#ifdef USE_XERCES
// ============================================================
// Function : read
// ------------------------------------------------------------
// 
// ============================================================
int dMParser::read(sheet* s, std::string fileName)
{
    mySheet = s;

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
       sheetFiller((DOMNode*)doc->getDocumentElement());
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
// Function : sheetFiller
// ------------------------------------------------------------
// Parses xml files using xercesc, given the rootnode
// ------------------------------------------------------------
void dMParser::sheetFiller(DOMNode *rootnode)
{
    char *name = XMLString::transcode(rootnode->getNodeName());

    if (std::string(name) == std::string("Sheet")) {
      if (rootnode->hasAttributes()) {
        // get all the attributes of the node
        DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
        int nSize = pAttributes->getLength();
        for (int i = 0; i < nSize; ++i) {
          DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
          // get attribute name
          char *name2 = XMLString::transcode(pAttributeNode->getName());
          if (std::string(name2) == std::string("name")) {
            char *sheetName = XMLString::transcode(pAttributeNode->getNodeValue());
            mySheet->setName(std::string(sheetName));
          }
          delete name2;
        }
      }
      tableFiller(rootnode);
    }
    delete name;
}

// ============================================================
// Function : tableFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::tableFiller(DOMNode *rootnode)
{
    DOMNode *c;

    bool bName = false;
    bool bnCols = false;
    bool bnRows = false;
    bool bType = false;

    std::string tableName = "";
    int nCols = 0;
    int nRows = 0;
    std::string tableType = "d";

    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "Table") {
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *tName = XMLString::transcode(pAttributeNode->getNodeValue());
                tableName = std::string(tName);
                bName = true;
                delete tName;
              }
              if (std::string(name) == std::string("type")) {
                char *tType = XMLString::transcode(pAttributeNode->getNodeValue());
                tableType = std::string(tType);
                bType = true;
                delete tType;
              }
              if (std::string(name) == std::string("num_cols")) {
                char *nc = XMLString::transcode(pAttributeNode->getNodeValue());
                nCols = atoi(nc);
                bnCols = true;
                delete nc;
              }
              if (std::string(name) == std::string("num_rows")) {
                char *nr = XMLString::transcode(pAttributeNode->getNodeValue());
                nRows = atoi(nr);
                bnRows = true;
                delete nr;
              }
            }
          }

          if (tableType == "d") {
            myDoubleTable = mySheet->addTable();
            myDoubleTable->setName(tableName);
            myDoubleTable->setType(tableType);
            myDoubleTable->setNumColumns(nCols);
            myDoubleTable->setNumRows(nRows);
            myDoubleTable->setup();
            currentType = "d";
          }
          else if (tableType == "i") {
            myIntTable = mySheet->addIntTable();
            myIntTable->setName(tableName);
            myDoubleTable->setType(tableType); // should this be doubleTable ??
            myIntTable->setNumColumns(nCols);
            myIntTable->setNumRows(nRows);
            myIntTable->setup();
            currentType = "i";
          }
          rowColFiller(c);
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : rowColFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::rowColFiller(DOMNode *rootnode)
{
    DOMNode *c;
    std::string rowName = "";
    int rowId = -1;
    std::string colName = "";
    int colId = -1;

    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "Row") { // <Row id="0" name="CH2Cl2" >
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *rName = XMLString::transcode(pAttributeNode->getNodeValue());
                rowName = std::string(rName);
                delete rName;
              }
              if (std::string(name) == std::string("id")) {
                char *rId = XMLString::transcode(pAttributeNode->getNodeValue());
                rowId = atoi(rId);
                delete rId;
              }
            }
            if ((rowId > -1) and (rowName != "")) {
              if (currentType == "d") {
                myDoubleTable->setRowLabel(rowId, rowName);
              }
              else if (currentType == "i") {
                myIntTable->setRowLabel(rowId, rowName);
              }
            }
          }
          if (currentType == "d") {
            cellFiller(c, rowId);
          }
          else if (currentType == "i") {
            cellIntFiller(c, rowId);
          }
          rowName = "";
          rowId = -1;
        }
        else if ((std::string)(XC(c->getNodeName())) == "Col") { // <Col id="0" name="LD25" />
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("name")) {
                char *cName = XMLString::transcode(pAttributeNode->getNodeValue());
                colName = std::string(cName);
                delete cName;
              }
              if (std::string(name) == std::string("id")) {
                char *cId = XMLString::transcode(pAttributeNode->getNodeValue());
                colId = atoi(cId);
                delete cId;
              }
            }
            if ((colId > -1) and (colName != "")) {
              if (currentType == "d") {
                myDoubleTable->setColumnLabel(colId, colName);
              }
              else if (currentType == "i") {
                myIntTable->setColumnLabel(colId, colName);
              }
            }
            colName = "";
            colId = -1;
          }
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : cellFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::cellFiller(DOMNode *rootnode, int rowId)
{
    DOMNode *c;
    int colId = -1;
    double cellValue = 0.0;

    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "Cell") { //    <Cell col="0" value="0.96" />
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("value")) {
                char *cValue = XMLString::transcode(pAttributeNode->getNodeValue());
                cellValue = strtod(cValue, 0);
              }
              if (std::string(name) == std::string("col")) {
                char *cId = XMLString::transcode(pAttributeNode->getNodeValue());
                colId = atoi(cId);
              }
              delete name;
            }
          }
          if ((colId > -1)) {
            myDoubleTable->setCellValue(rowId, colId, cellValue);
          }
          cellValue = 0.0;
          colId = -1;
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : cellIntFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::cellIntFiller(DOMNode *rootnode, int rowId)
{
    DOMNode *c;
    int colId = -1;
    int cellValue = 0;

    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "Cell") { //    <Cell col="0" value="0.96" />
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("value")) {
                char *cValue = XMLString::transcode(pAttributeNode->getNodeValue());
                cellValue = atoi(cValue);
              }
              if (std::string(name) == std::string("col")) {
                char *cId = XMLString::transcode(pAttributeNode->getNodeValue());
                colId = atoi(cId);
              }
              delete name;
            }
          }
          if ((colId > -1)) {
            myIntTable->setCellValue(rowId, colId, cellValue);
          }
          cellValue = 0;
          colId = -1;
        }
        else {
        }
      }
      else {
      }
    }
}
#endif // USE_XERCES

    // ----------------------------- //
    // -    Qt WRITE FUNCTIONS     - //
    // ----------------------------- //
#ifdef USE_QT
// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write library xml files using Qt
// ------------------------------------------------------------
int dMParser::write(sheet* s, std::string fileName, bool bComments)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("dMParser", errMessage, INFO);

    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QDomElement rootElem = doc.createElement("Sheet");
    rootElem.setAttribute("name", QString::fromStdString(s->getName()));
    doc.appendChild(rootElem);

    std::vector<table<double>*> doubleTables = s->getTables();
    std::vector<table<int>*> intTables = s->getIntTables();

    for (unsigned int i = 0; i < doubleTables.size(); i++) {
      this->writeTable(doc, doubleTables[i]);
    }

    for (unsigned int i = 0; i < intTables.size(); i++) {
      this->writeIntTable(doc, intTables[i]);
    }

    // Write dom xml to file
    QString domXml = doc.toString();

    QFile file(qFileName);
    if (!file.open(QIODevice::WriteOnly)) {
      errMessage = " Writing " + fileName;
      errorLogger.throwError("dMParser", errMessage, MTK_ERROR);
      return 1;
    }

    QTextStream out(&file);
    out << domXml;
    file.close();
    return 0;
}

// ============================================================
// Function : writeTable
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::writeTable(QDomDocument doc, table<double>* t)
{
    QDomElement rootElem = doc.documentElement();
    QDomElement tableElem = doc.createElement("Table");
    rootElem.appendChild(tableElem);

    tableElem.setAttribute("name", QString::fromStdString(t->getName()));
    tableElem.setAttribute("num_rows", QString::number(t->getNumRows()));
    tableElem.setAttribute("num_cols", QString::number(t->getNumColumns()));

    std::vector<std::string> colLabels = t->getColumnLabels();
    std::vector<std::string> rowLabels = t->getRowLabels();

    //ublas::matrix<double> &tableMatrix = t->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &tableMatrix = t->getMatrix();

    for (unsigned int i = 0; i < tableMatrix.rows(); i++) {
      QDomElement rowElem = doc.createElement("Row");
      tableElem.appendChild(rowElem);
      rowElem.setAttribute("id", QString::number(i));
      rowElem.setAttribute("name", QString::fromStdString(rowLabels[i].c_str()));

      for (unsigned int j = 0; j < tableMatrix.cols(); j++) {
        QDomElement cellElem = doc.createElement("Cell");
        rowElem.appendChild(cellElem);
        cellElem.setAttribute("col", QString::number(j));
        cellElem.setAttribute("value", QString::number(tableMatrix(i,j)));
      }
    }
    for (unsigned int i = 0; i < tableMatrix.cols(); i++) {
      QDomElement colElem = doc.createElement("Col");
      tableElem.appendChild(colElem);
      colElem.setAttribute("id", QString::number(i));
      colElem.setAttribute("name", QString::fromStdString(colLabels[i]));
    }
}

// ============================================================
// Function : writeTable
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::writeIntTable(QDomDocument doc, table<int>* t)
{
    QDomElement rootElem = doc.documentElement();
    QDomElement tableElem = doc.createElement("Table");
    rootElem.appendChild(tableElem);

    rootElem.setAttribute("name", QString::fromStdString(t->getName()));

    tableElem.setAttribute("name", QString::fromStdString(t->getName()));
    tableElem.setAttribute("num_rows", QString::number(t->getNumRows()));
    tableElem.setAttribute("num_cols", QString::number(t->getNumColumns()));

    std::vector<std::string> colLabels = t->getColumnLabels();
    std::vector<std::string> rowLabels = t->getRowLabels();

    //ublas::matrix<int> &tableMatrix = t->getMatrix();
    Eigen::Matrix<int, Dynamic, Dynamic> &tableMatrix = t->getMatrix();

    for (unsigned int i = 0; i < tableMatrix.rows(); i++) {
      QDomElement rowElem = doc.createElement("Row");
      tableElem.appendChild(rowElem);
      rowElem.setAttribute("id", QString::number(i));
      rowElem.setAttribute("name", QString::fromStdString(rowLabels[i].c_str()));

      for (unsigned int j = 0; j < tableMatrix.cols(); j++) {
        QDomElement cellElem = doc.createElement("Cell");
        rowElem.appendChild(cellElem);
        cellElem.setAttribute("col", QString::number(j));
        cellElem.setAttribute("value", QString::number(tableMatrix(i,j)));
      }
    }
    for (unsigned int i = 0; i < tableMatrix.cols(); i++) {
      QDomElement colElem = doc.createElement("Col");
      tableElem.appendChild(colElem);
      colElem.setAttribute("id", QString::number(i));
      colElem.setAttribute("name", QString::fromStdString(colLabels[i]));
    }
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
int dMParser::write(sheet* s, std::string fileName, bool bComments)
{
    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

    TiXmlElement* root = new TiXmlElement( "Sheet" );
    doc.LinkEndChild( root );
    root->SetAttribute("name", s->getName());

    std::vector<table<double>*> doubleTables = s->getTables();
    std::vector<table<int>*> intTables = s->getIntTables();

    if (doubleTables.size() > 0) {

      if (bComments) {
        TiXmlComment* commentDoubles = new TiXmlComment();
        commentDoubles->SetValue(" Table containing Doubles " );
        root->LinkEndChild(commentDoubles);
      }

      for (unsigned int i = 0; i < doubleTables.size(); i++) {
        this->writeTable(root, doubleTables[i]);
      }
    }

    if (intTables.size() > 0) {
      if (bComments) {
        TiXmlComment* commentInts = new TiXmlComment();
        commentInts->SetValue(" Table containing Integers " );
        root->LinkEndChild(commentInts);
      }

      for (unsigned int i = 0; i < intTables.size(); i++) {
        this->writeIntTable(root, intTables[i]);
      }
    }

    doc.SaveFile(fileName);

/*
    TiXmlElement * msgs = new TiXmlElement( "Messages" );
    root->LinkEndChild( msgs );

    msg = new TiXmlElement( "Welcome" );
    msg->LinkEndChild( new TiXmlText( "Welcome to MyApp" ));
    msgs->LinkEndChild( msg );

    msg = new TiXmlElement( "Farewell" );
    msg->LinkEndChild( new TiXmlText( "Thank you for using MyApp" ));
    msgs->LinkEndChild( msg );

    TiXmlElement * windows = new TiXmlElement( "Windows" );
    root->LinkEndChild( windows );

    TiXmlElement * window;
    window = new TiXmlElement( "Window" );
    windows->LinkEndChild( window );
    window->SetAttribute("name", "MainFrame");
    window->SetAttribute("x", 5);
    window->SetAttribute("y", 15);
    window->SetAttribute("w", 400);
    window->SetAttribute("h", 250);

    TiXmlElement * cxn = new TiXmlElement( "Connection" );
    root->LinkEndChild( cxn );
    cxn->SetAttribute("ip", "192.168.0.1");
    cxn->SetDoubleAttribute("timeout", 123.456); // floating point attrib
*/
/*
<?xml version="1.0" ?>
<MyApp>
    <!-- Settings for MyApp -->
    <Messages>
        <Welcome>Welcome to MyApp</Welcome>
        <Farewell>Thank you for using MyApp</Farewell>
    </Messages>
    <Windows>
        <Window name="MainFrame" x="5" y="15" w="400" h="250" />
    </Windows>
    <Connection ip="192.168.0.1" timeout="123.456000" />
</MyApp>
*/

    //dump_to_stdout(&doc);

    return 0;
}

// ============================================================
// Function : writeTable
// ------------------------------------------------------------
// Write table using tinyxml
// ------------------------------------------------------------
void dMParser::writeTable(TiXmlElement* root, table<double>* t)
{
    TiXmlElement* table = new TiXmlElement("Table");
    root->LinkEndChild(table);

    table->SetAttribute("name", t->getName());
    table->SetAttribute("num_rows", int2String(t->getNumRows()));
    table->SetAttribute("num_cols", int2String(t->getNumColumns()));
    table->SetAttribute("type", "d");

    std::vector<std::string> colLabels = t->getColumnLabels();
    std::vector<std::string> rowLabels = t->getRowLabels();

    //ublas::matrix<double> &tableMatrix = t->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &tableMatrix = t->getMatrix();

    for (unsigned int i = 0; i < tableMatrix.rows(); i++) {

      TiXmlElement* row = new TiXmlElement("Row");
      table->LinkEndChild(row);

      row->SetAttribute("id", int2String(i).c_str());
      row->SetAttribute("name", rowLabels[i].c_str());

      for (unsigned int j = 0; j < tableMatrix.cols(); j++) {
        TiXmlElement* cell = new TiXmlElement("Cell");
        row->LinkEndChild(cell);

        cell->SetAttribute("col", int2String(j).c_str());
        cell->SetAttribute("value", double2String(tableMatrix(i,j)).c_str());
        //cell->SetDoubleAttribute("value", tableMatrix(i,j));
      }
    }

    for (unsigned int i = 0; i < tableMatrix.cols(); i++) {
      TiXmlElement* col = new TiXmlElement("Col");
      table->LinkEndChild(col);

      col->SetAttribute("id", int2String(i).c_str());
      col->SetAttribute("name", colLabels[i].c_str());
    }
}

// ============================================================
// Function : writeIntTable
// ------------------------------------------------------------
// Write table using tinyxml
// ------------------------------------------------------------
void dMParser::writeIntTable(TiXmlElement* root, table<int>* t)
{
    TiXmlElement* table = new TiXmlElement("Table");
    root->LinkEndChild(table);

    table->SetAttribute("name", t->getName());
    table->SetAttribute("num_rows", int2String(t->getNumRows()));
    table->SetAttribute("num_cols", int2String(t->getNumColumns()));
    table->SetAttribute("type", "i");

    std::vector<std::string> colLabels = t->getColumnLabels();
    std::vector<std::string> rowLabels = t->getRowLabels();

    //ublas::matrix<int> &tableMatrix = t->getMatrix();
    Eigen::Matrix<int, Dynamic, Dynamic> &tableMatrix = t->getMatrix();

    for (unsigned int i = 0; i < tableMatrix.rows(); i++) {

      TiXmlElement* row = new TiXmlElement("Row");
      table->LinkEndChild(row);

      row->SetAttribute("id", int2String(i).c_str());
      row->SetAttribute("name", rowLabels[i].c_str());

      for (unsigned int j = 0; j < tableMatrix.cols(); j++) {
        TiXmlElement* cell = new TiXmlElement("Cell");
        row->LinkEndChild(cell);

        cell->SetAttribute("col", int2String(j).c_str());
        cell->SetAttribute("value", int2String(tableMatrix(i,j)).c_str());
      }
    }

    for (unsigned int i = 0; i < tableMatrix.cols(); i++) {
      TiXmlElement* col = new TiXmlElement("Col");
      table->LinkEndChild(col);

      col->SetAttribute("id", int2String(i).c_str());
      col->SetAttribute("name", colLabels[i].c_str());
    }
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
int dMParser::write(sheet* s, std::string fileName, bool bComments)
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
          DOMLSSerializer *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer(XMLPlatformUtils::fgMemoryManager);
#else
          DOMWriter *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
#endif
          //DOMErrorHandler *myErrorHandler = new DOMPrintErrorHandler();
          //theSerializer->setErrorHandler(myErrorHandler);

          XERCES_CPP_NAMESPACE::DOMDocument* doc = impl->createDocument(
            0,           // root element namespace URI.
            X("Sheet"), // root element name
            0);          // document type object (DTD).

          DOMElement* rootElem = doc->getDocumentElement();

          rootElem->setAttribute(X("name"), X(s->getName().c_str()));

          std::vector<table<double>*> doubleTables = s->getTables();
          std::vector<table<int>*> intTables = s->getIntTables();

          for (unsigned int i = 0; i < doubleTables.size(); i++) {
            this->writeTable(doc, doubleTables[i]);
          }

          for (unsigned int i = 0; i < intTables.size(); i++) {
            this->writeIntTable(doc, intTables[i]);
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
    return errorCode;
}

// ============================================================
// Function : writeTable
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::writeTable(XERCES_CPP_NAMESPACE::DOMDocument* doc, table<double>* t)
{
    DOMElement* rootElem = doc->getDocumentElement();
    DOMElement* tableElem = doc->createElement(X("Table"));
    rootElem->appendChild(tableElem);

    tableElem->setAttribute(X("name"), X(t->getName().c_str()));
    tableElem->setAttribute(X("num_rows"), X(int2String(t->getNumRows()).c_str()));
    tableElem->setAttribute(X("num_cols"), X(int2String(t->getNumColumns()).c_str()));

    std::vector<std::string> colLabels = t->getColumnLabels();
    std::vector<std::string> rowLabels = t->getRowLabels();

    //ublas::matrix<double> &tableMatrix = t->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> &tableMatrix = t->getMatrix();

    for (unsigned int i = 0; i < tableMatrix.rows(); i++) {
      DOMElement* rowElem = doc->createElement(X("Row"));
      tableElem->appendChild(rowElem);
      rowElem->setAttribute(X("id"), X(int2String(i).c_str()));
      rowElem->setAttribute(X("name"), X(rowLabels[i].c_str()));
      for (unsigned int j = 0; j < tableMatrix.cols(); j++) {
        DOMElement* cellElem = doc->createElement(X("Cell"));
        rowElem->appendChild(cellElem);
        cellElem->setAttribute(X("col"), X(int2String(j).c_str()));
        cellElem->setAttribute(X("value"), X(int2String(j).c_str()));
        cellElem->setAttribute(X("value"), X(double2String(tableMatrix(i,j)).c_str()));
      }
    }
    for (unsigned int i = 0; i < tableMatrix.cols(); i++) {
      DOMElement* colElem = doc->createElement(X("Col"));
      tableElem->appendChild(colElem);
      colElem->setAttribute(X("id"), X(int2String(i).c_str()));
      colElem->setAttribute(X("name"), X(colLabels[i].c_str()));
    }
}

// ============================================================
// Function : writeIntTable
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void dMParser::writeIntTable(XERCES_CPP_NAMESPACE::DOMDocument* doc, table<int>* t)
{
    DOMElement* rootElem = doc->getDocumentElement();
    DOMElement* tableElem = doc->createElement(X("Table"));
    rootElem->appendChild(tableElem);

    tableElem->setAttribute(X("name"), X(t->getName().c_str()));
    tableElem->setAttribute(X("num_rows"), X(int2String(t->getNumRows()).c_str()));
    tableElem->setAttribute(X("num_cols"), X(int2String(t->getNumColumns()).c_str()));

    std::vector<std::string> colLabels = t->getColumnLabels();
    std::vector<std::string> rowLabels = t->getRowLabels();

    //ublas::matrix<int> &tableMatrix = t->getMatrix();
    Eigen::Matrix<int, Dynamic, Dynamic> &tableMatrix = t->getMatrix();

    for (unsigned int i = 0; i < tableMatrix.rows(); i++) {
      DOMElement* rowElem = doc->createElement(X("Row"));
      tableElem->appendChild(rowElem);
      rowElem->setAttribute(X("id"), X(int2String(i).c_str()));
      rowElem->setAttribute(X("name"), X(rowLabels[i].c_str()));
      for (unsigned int j = 0; j < tableMatrix.cols(); j++) {
        DOMElement* cellElem = doc->createElement(X("Cell"));
        rowElem->appendChild(cellElem);
        cellElem->setAttribute(X("col"), X(int2String(j).c_str()));
        cellElem->setAttribute(X("value"), X(int2String(j).c_str()));
        cellElem->setAttribute(X("value"), X(int2String(tableMatrix(i,j)).c_str()));
      }
    }
    for (unsigned int i = 0; i < tableMatrix.cols(); i++) {
      DOMElement* colElem = doc->createElement(X("Col"));
      tableElem->appendChild(colElem);
      colElem->setAttribute(X("id"), X(int2String(i).c_str()));
      colElem->setAttribute(X("name"), X(colLabels[i].c_str()));
    }
}
#endif // USE_XERCES

} // MTKpp namespace
